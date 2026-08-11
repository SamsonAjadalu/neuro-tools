#!/usr/bin/env python3
"""Create an OPPNI-D input file from a BIDS diffusion dataset."""

import argparse
import fnmatch
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_DWI_PATTERNS = ["*dwi*.nii.gz", "*dwi*.nii"]
DEFAULT_ANAT_PATTERNS = ["*T1w*.nii.gz", "*T1w*.nii"]
DEFAULT_REVERSE_PATTERNS = ["*dir-PA*dwi*", "*dir-rev*dwi*", "*acq-rev*dwi*", "*b0ref*dwi*"]
REV_MODE_CHOICES = {"NONE", "REF"}
PATTERN_LIST_SECTIONS = {
    "anat_patterns",
    "func_patterns",
    "dwi_patterns",
    "func_exclude_patterns",
    "reverse_pe_patterns",
    "func_reverse_pe_patterns",
    "dwi_reverse_pe_patterns",
    "fieldmap_magnitude_patterns",
    "fieldmap_phasediff_patterns",
}


def image_stem(path):
    if path.name.endswith(".nii.gz"):
        return path.name[:-7]
    if path.name.endswith(".nii"):
        return path.name[:-4]
    return path.stem


def sidecar_for(path):
    return path.with_name(image_stem(path) + ".json")


def matching_files(folder, patterns):
    if not folder.is_dir():
        return []
    files = []
    for pattern in patterns:
        for path in sorted(folder.glob(pattern)):
            if path.is_file() and not path.name.startswith(".") and path not in files:
                files.append(path)
    return sorted(files)


def matches_any(path, patterns):
    return any(fnmatch.fnmatchcase(path.name, pattern) for pattern in patterns)


def is_metadata_section(section):
    return (
        section.lower() == "defaults"
        or re.match(r"^(?:bold|dwi)\s+manufacturer:", section, re.IGNORECASE)
        or section.lower().startswith("manufacturer:")
    )


def load_patterns(pattern_file):
    patterns = {
        "anat_patterns": list(DEFAULT_ANAT_PATTERNS),
        "dwi_patterns": list(DEFAULT_DWI_PATTERNS),
        "dwi_reverse_pe_patterns": list(DEFAULT_REVERSE_PATTERNS),
        "metadata_rules": {},
    }
    if pattern_file is None:
        return patterns
    section = None
    seen = set()
    with Path(pattern_file).expanduser().open(encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip()
                if section not in PATTERN_LIST_SECTIONS and not is_metadata_section(section):
                    valid = ", ".join(sorted(PATTERN_LIST_SECTIONS)) + ", [dwi manufacturer:NAME]"
                    raise ValueError(f"unknown pattern section [{section}] on line {line_number}; valid sections: {valid}")
                if section not in seen:
                    if section in PATTERN_LIST_SECTIONS:
                        patterns[section] = []
                    else:
                        patterns.setdefault("metadata_rules", {})[section] = {}
                    seen.add(section)
                continue
            if section is None:
                raise ValueError(f"pattern line outside a section on line {line_number}: {line}")
            if is_metadata_section(section):
                if section.lower().startswith("bold "):
                    continue
                if "=" not in line:
                    raise ValueError(f"expected key=value in [{section}] on line {line_number}")
                key, value = line.split("=", 1)
                key = key.strip().lower()
                if key not in {
                    "pe_fwd",
                    "tro_msec",
                    "rev_mode",
                    "echo_spacing_msec",
                    "accel_factor",
                    "recon_matrix_pe",
                }:
                    raise ValueError(f"unknown DWI rule field {key!r} on line {line_number}")
                patterns["metadata_rules"][section][key] = value.strip()
            elif section in PATTERN_LIST_SECTIONS:
                patterns[section].append(line)
    return patterns


def normalize_manufacturer(value):
    value = str(value or "").strip().lower()
    if "siemens" in value:
        return "Siemens"
    if "philips" in value:
        return "Philips"
    if value in {"ge", "ge medical systems"} or "general electric" in value:
        return "GE"
    return str(value).strip()


def normalize_text(value):
    return " ".join(str(value or "").split()).casefold()


def parse_dwi_rule_section(section):
    body = section.strip()
    if body.lower().startswith("dwi "):
        body = body[4:].strip()
    if body.lower() == "defaults":
        return {"specificity": 0, "manufacturer": "", "model": "", "protocol": "", "tr": None, "te": None}

    match = re.fullmatch(
        r"manufacturer:(?P<manufacturer>.*?)\s+model:(?P<model>.*?)\s+protocol:(?P<protocol>.*?)"
        r"(?:\s+tr:(?P<tr>\S+))?(?:\s+te:(?P<te>\S+))?",
        body,
        re.IGNORECASE,
    )
    if match:
        return {
            "specificity": 3 + bool(match.group("tr")) + bool(match.group("te")),
            "manufacturer": normalize_manufacturer(match.group("manufacturer")),
            "model": normalize_text(match.group("model")),
            "protocol": normalize_text(match.group("protocol")),
            "tr": match.group("tr"),
            "te": match.group("te"),
        }

    match = re.fullmatch(r"manufacturer:(.+)", body, re.IGNORECASE)
    if match:
        return {
            "specificity": 1,
            "manufacturer": normalize_manufacturer(match.group(1)),
            "model": "",
            "protocol": "",
            "tr": None,
            "te": None,
        }
    return None


def numeric_match(rule_value, metadata_value):
    if rule_value is None:
        return True
    try:
        return math.isclose(float(rule_value), float(metadata_value), rel_tol=0, abs_tol=1e-6)
    except (TypeError, ValueError):
        return False


def dwi_rule_values(rules, metadata):
    manufacturer = normalize_manufacturer(metadata.get("Manufacturer"))
    model = normalize_text(metadata.get("ManufacturersModelName"))
    protocol = normalize_text(metadata.get("ProtocolName"))
    values = {}
    sources = {}
    for section, fields in rules.items():
        rule = parse_dwi_rule_section(section)
        if rule is None:
            continue
        specificity = rule["specificity"]
        matches = (
            (not rule["manufacturer"] or manufacturer.casefold() == rule["manufacturer"].casefold())
            and (not rule["model"] or model == rule["model"])
            and (not rule["protocol"] or protocol == rule["protocol"])
            and numeric_match(rule["tr"], metadata.get("RepetitionTime"))
            and numeric_match(rule["te"], metadata.get("EchoTime"))
        )
        source = section
        if not matches:
            continue
        for key, value in fields.items():
            if key in values and values[key] != value and sources[key][0] == specificity:
                raise ValueError(f"conflicting DWI rules for {key}: {sources[key][1]} and {source}")
            if key not in values or specificity >= sources[key][0]:
                values[key] = value
                sources[key] = (specificity, source)
    return values


def find_sessions(subject_folder):
    sessions = sorted(path for path in subject_folder.glob("ses-*") if path.is_dir())
    if sessions:
        return sessions
    if (subject_folder / "anat").is_dir() or (subject_folder / "dwi").is_dir():
        return [subject_folder]
    return []


def label_for(subject_folder, session_folder):
    if session_folder == subject_folder:
        return subject_folder.name
    return f"{subject_folder.name}_{session_folder.name}"


def load_json(path):
    try:
        with path.open(encoding="utf-8") as file:
            data = json.load(file)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read JSON sidecar {path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"JSON sidecar must contain an object: {path}")
    return data


def load_nifti(path):
    try:
        import nibabel as nib
    except ImportError as error:
        raise ValueError(f"cannot validate NIfTI files without nibabel: {path}") from error
    try:
        return nib.load(str(path))
    except Exception as error:
        raise ValueError(f"cannot read DWI NIfTI {path}: {error}") from error


def read_numeric_file(path, label):
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise ValueError(f"cannot read {label} file {path}: {error}") from error
    try:
        values = [float(value) for value in text.split()]
    except ValueError as error:
        raise ValueError(f"{label} contains nonnumeric values: {path}") from error
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError(f"{label} is empty or contains non-finite values: {path}")
    return values


def validate_dwi(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty DWI: {path}")
    image = load_nifti(path)
    if len(image.shape) != 4:
        raise ValueError(f"DWI must be 4D: {path} has shape {image.shape}")
    bval = path.with_name(image_stem(path) + ".bval")
    bvec = path.with_name(image_stem(path) + ".bvec")
    if not bval.is_file() or bval.stat().st_size == 0:
        raise ValueError(f"missing or empty bval for {path}: {bval}")
    if not bvec.is_file() or bvec.stat().st_size == 0:
        raise ValueError(f"missing or empty bvec for {path}: {bvec}")
    bvals = read_numeric_file(bval, "bval")
    if len(bvals) != image.shape[3]:
        raise ValueError(f"bval count {len(bvals)} does not match DWI volumes {image.shape[3]}: {path}")
    if any(value < 0 for value in bvals):
        raise ValueError(f"bval contains negative values: {bval}")
    try:
        raw_rows = [line.split() for line in bvec.read_text(encoding="utf-8").splitlines() if line.strip()]
        rows = [[float(value) for value in row] for row in raw_rows]
    except (OSError, ValueError) as error:
        raise ValueError(f"bvec contains invalid values: {bvec}") from error
    if len(rows) != 3 or any(len(row) != image.shape[3] for row in rows):
        raise ValueError(f"bvec must have shape 3 x {image.shape[3]}: {bvec}")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise ValueError(f"bvec contains non-finite values: {bvec}")
    return image.shape, bval, bvec


def check_dwi_sidecars(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty DWI: {path}")
    bval = path.with_name(image_stem(path) + ".bval")
    bvec = path.with_name(image_stem(path) + ".bvec")
    if not bval.is_file() or bval.stat().st_size == 0:
        raise ValueError(f"missing or empty bval for {path}: {bval}")
    if not bvec.is_file() or bvec.stat().st_size == 0:
        raise ValueError(f"missing or empty bvec for {path}: {bvec}")
    return None, bval, bvec


def axis_direction(image, phase_encoding):
    """Convert a BIDS phase-encoding axis to an OPPNI world-axis direction."""
    import nibabel as nib

    direction = str(phase_encoding or "")
    if not direction or direction[0].lower() not in "ijk" or direction[1:] not in {"", "-"}:
        raise ValueError(f"invalid PhaseEncodingDirection {phase_encoding!r}; expected i, i-, j, j-, k, or k-")
    axis = "ijk".index(direction[0].lower())
    axcodes = nib.aff2axcodes(image.affine)
    positive = axcodes[axis]
    if positive not in {"R", "L", "A", "P", "S", "I"}:
        raise ValueError(f"cannot map NIfTI phase-encoding axis {axis} to PE direction")
    pairs = {"R": ("L", "R"), "L": ("R", "L"), "A": ("P", "A"), "P": ("A", "P"), "S": ("I", "S"), "I": ("S", "I")}
    source, target = pairs[positive]
    if direction.endswith("-"):
        source, target = target, source
    return f"{source}>>{target}"


def resolve_pe(image, metadata, explicit, fallback, path):
    value = explicit or metadata.get("PhaseEncodingDirection") or fallback
    if not value:
        raise ValueError(f"missing PE_FWD/PhaseEncodingDirection for {path}")
    if str(value).upper() in {"A>>P", "P>>A", "R>>L", "L>>R", "I>>S", "S>>I"}:
        return str(value).upper()
    if image is None:
        image = load_nifti(path)
    return axis_direction(image, value)


def resolve_tro(metadata, explicit, fallback, path):
    if explicit is not None:
        try:
            value = float(explicit)
        except ValueError as error:
            raise ValueError(f"invalid TRO_MSEC {explicit!r} for {path}") from error
        if value <= 0:
            raise ValueError(f"TRO_MSEC must be positive for {path}")
        return f"{value:g}"
    value = metadata.get("TotalReadoutTime")
    if value is None and fallback is not None:
        if fallback.get("tro_msec") is not None:
            value = fallback["tro_msec"]
        else:
            try:
                echo_spacing = float(fallback["echo_spacing_msec"])
                accel_factor = float(fallback["accel_factor"])
                recon_matrix_pe = int(fallback["recon_matrix_pe"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"rule must provide tro_msec or echo_spacing_msec, accel_factor, and recon_matrix_pe for {path}"
                ) from error
            if echo_spacing <= 0 or accel_factor <= 0 or recon_matrix_pe <= 1:
                raise ValueError(f"invalid DWI timing rule for {path}")
            # Convert protocol echo spacing to OPPNI total readout time.
            value = echo_spacing / accel_factor * (recon_matrix_pe - 1)
        try:
            value = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid rule TRO_MSEC {fallback!r} for {path}") from error
        if value <= 0:
            raise ValueError(f"rule TRO_MSEC must be positive for {path}")
        return f"{value:g}"
    if value is None:
        raise ValueError(f"missing TRO_MSEC: JSON lacks TotalReadoutTime for {path}")
    try:
        seconds = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid TotalReadoutTime for {path}") from error
    if seconds <= 0:
        raise ValueError(f"TotalReadoutTime must be positive for {path}")
    return f"{seconds * 1000:g}"


def find_anat(subject_folder, session_folder, patterns):
    anat = matching_files(session_folder / "anat", patterns)
    if anat:
        return anat[0]
    if session_folder != subject_folder:
        anat = matching_files(subject_folder / "anat", patterns)
        if anat:
            return anat[0]
    return None


def record_for(subject_folder, session_folder, args, patterns):
    dwi_patterns = args.dwi_pattern or patterns["dwi_patterns"]
    reverse_patterns = args.reverse_dwi_pattern or patterns["dwi_reverse_pe_patterns"]
    anat_patterns = args.anat_pattern or patterns["anat_patterns"]
    dwi_files = matching_files(session_folder / "dwi", dwi_patterns)
    if not dwi_files:
        raise SkipRecord("missing DWI NIfTI")
    reverse_files = [path for path in dwi_files if matches_any(path, reverse_patterns)]
    forward_files = [path for path in dwi_files if path not in reverse_files]
    if not forward_files:
        forward_files = dwi_files
        reverse_files = []
    anat = find_anat(subject_folder, session_folder, anat_patterns)
    if anat is None:
        raise SkipRecord("missing T1 anatomy")
    if not anat.is_file() or anat.stat().st_size == 0:
        raise SkipRecord(f"missing or empty T1 anatomy: {anat}")
    if args.validate:
        anat_image = load_nifti(anat)
        if len(anat_image.shape) != 3:
            raise ValueError(f"T1 anatomy must be 3D: {anat} has shape {anat_image.shape}")

    validated = []
    rule_values = None
    for path in forward_files:
        if args.validate:
            shape, _, _ = validate_dwi(path)
        else:
            shape, _, _ = check_dwi_sidecars(path)
        metadata = load_json(sidecar_for(path)) if sidecar_for(path).is_file() else None
        if metadata is None:
            raise ValueError(f"missing JSON sidecar for DWI: {path}")
        current_rules = dwi_rule_values(patterns["metadata_rules"], metadata)
        if rule_values is None:
            rule_values = current_rules
        pe = resolve_pe(None, metadata, args.pe_fwd, current_rules.get("pe_fwd"), path)
        tro = resolve_tro(metadata, args.tro_msec, current_rules, path)
        validated.append((path, pe, tro, shape))

    pe_values = {(item[1], item[2]) for item in validated}
    if len(pe_values) != 1:
        raise ValueError(f"forward DWI runs disagree on PE_FWD/TRO_MSEC in {label_for(subject_folder, session_folder)}")
    pe_fwd, tro_msec = next(iter(pe_values))

    reverse = None
    pe_rev = None
    rev_mode = str(args.rev_mode or (rule_values or {}).get("rev_mode") or "NONE").upper()
    if rev_mode not in REV_MODE_CHOICES:
        raise ValueError(f"invalid REV_MODE {rev_mode!r}; expected NONE or REF")
    if rev_mode == "REF":
        if not reverse_files:
            raise SkipRecord("REV_MODE=REF but no reverse-PE DWI was found")
        reverse = reverse_files[0]
        if args.validate:
            validate_dwi(reverse)
        else:
            check_dwi_sidecars(reverse)
        reverse_json = sidecar_for(reverse)
        if not reverse_json.is_file():
            raise ValueError(f"missing JSON sidecar for reverse DWI: {reverse}")
        pe_rev = resolve_pe(None, load_json(reverse_json), args.pe_rev, None, reverse)
    return {
        "prefix": label_for(subject_folder, session_folder),
        "anat": anat,
        "forward": [item[0] for item in validated],
        "reverse": reverse,
        "pe_fwd": pe_fwd,
        "pe_rev": pe_rev,
        "tro_msec": tro_msec,
        "rev_mode": rev_mode,
        "shapes": [item[3] for item in validated],
    }


class SkipRecord(Exception):
    """A missing acquisition that can be reported and skipped."""


def row_for_record(record, args):
    fields = [
        f"PREFIX={record['prefix']}",
        "DIFF_FWD=" + ",".join(str(path) for path in record["forward"]),
        f"ANAT={record['anat']}",
        f"ZCLIP={args.zclip}",
        f"TRO_MSEC={record['tro_msec']}",
        f"PE_FWD={record['pe_fwd']}",
        f"REV_MODE={record['rev_mode']}",
    ]
    if record["reverse"] is not None:
        fields.insert(2, f"DIFF_REV={record['reverse']}")
        fields.append(f"PE_REV={record['pe_rev']}")
    return " ".join(fields)


def inspect_record(record):
    lines = [
        f"{record['prefix']}:",
        f"  DIFF_FWD: {len(record['forward'])} file(s)",
    ]
    for path, shape in zip(record["forward"], record["shapes"]):
        lines.append(f"    - {path} shape={shape}")
    lines.extend([
        f"  DIFF_REV: {record['reverse'] or 'none'}",
        f"  ANAT: {record['anat']}",
        f"  TRO_MSEC: {record['tro_msec']}",
        f"  PE_FWD: {record['pe_fwd']}",
        f"  PE_REV: {record['pe_rev'] or 'none'}",
        f"  REV_MODE: {record['rev_mode']}",
    ])
    return "\n".join(lines)


def build_records(args, patterns):
    root = Path(args.bids_root).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"ERROR: missing BIDS root: {root}")
    subjects = [path for path in sorted(root.glob("sub-*")) if path.is_dir() and (not args.subject or args.subject in path.name)]
    records = []
    skipped = []
    errors = []
    sessions = 0
    for index, subject in enumerate(subjects, start=1):
        print(f"[{index}/{len(subjects)}] Processing {subject.name}", flush=True)
        for session in find_sessions(subject):
            if args.session and args.session not in ("" if session == subject else session.name):
                continue
            sessions += 1
            prefix = label_for(subject, session)
            try:
                records.append(record_for(subject, session, args, patterns))
            except SkipRecord as error:
                skipped.append((prefix, str(error)))
                print(f"[{index}/{len(subjects)}] SKIPPING {prefix}: {error}", flush=True)
            except ValueError as error:
                errors.append(f"{prefix}: {error}")
                print(f"[{index}/{len(subjects)}] ERROR for {prefix}: {error}", flush=True)
                print(f"[{index}/{len(subjects)}] SKIPPING {prefix}: invalid DWI inputs", flush=True)
    print(
        f"Discovery summary: subjects={len(subjects)} sessions={sessions} rows={len(records)} "
        f"skipped={len(skipped)} errors={len(errors)}",
        flush=True,
    )
    return records, skipped, errors


def parse_args():
    parser = argparse.ArgumentParser(description="Generate OPPNI-D input rows from a BIDS diffusion dataset")
    required = parser.add_argument_group("required")
    required.add_argument("--bids-root", required=True, metavar="PATH", help="BIDS dataset root")
    required.add_argument("--output", required=True, metavar="FILE", help="Output OPPNI-D input file")
    optional = parser.add_argument_group("optional")
    optional.add_argument("--patterns", metavar="FILE", help="Shared dataset filename and DWI metadata rules")
    optional.add_argument("--dwi-pattern", action="append", help="DWI filename pattern; repeatable")
    optional.add_argument("--anat-pattern", action="append", help="T1 filename pattern; repeatable")
    optional.add_argument("--reverse-dwi-pattern", action="append", help="Reverse-PE DWI pattern; repeatable")
    optional.add_argument("--pe-fwd", metavar="VALUE", help="Override PE_FWD, for example A>>P")
    optional.add_argument("--pe-rev", metavar="VALUE", help="Override PE_REV, for example P>>A")
    optional.add_argument("--tro-msec", metavar="VALUE", help="Override TRO_MSEC in milliseconds")
    optional.add_argument("--rev-mode", choices=sorted(REV_MODE_CHOICES), default=None, help="Reverse-PE mode")
    optional.add_argument("--zclip", default="AUTO", help="OPPNI-D ZCLIP value")
    optional.add_argument("--subject", help="Only include subjects whose folder name contains this text")
    optional.add_argument("--session", help="Only include sessions whose folder name contains this text")
    optional.add_argument("--inspect", action="store_true", help="Print resolved files and metadata")
    optional.add_argument(
        "--validate",
        action="store_true",
        help="Perform full NIfTI, bval, and bvec consistency validation",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    try:
        patterns = load_patterns(args.patterns)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    records, skipped, errors = build_records(args, patterns)
    if errors:
        print("ERROR: one or more acquisitions failed validation; no output was written", file=sys.stderr)
        sys.exit(1)
    if args.inspect:
        for record in records:
            print(inspect_record(record))
    output = Path(args.output).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(row_for_record(record, args) for record in records) + ("\n" if records else ""), encoding="utf-8")
    print(f"Wrote {len(records)} row(s) to {output}")
    if skipped:
        print(f"Skipped {len(skipped)} record(s)")


if __name__ == "__main__":
    main()
