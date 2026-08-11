#!/usr/bin/env python3
"""Create an OPPNI-B input file from a BIDS-like dataset."""

import argparse
import fnmatch
import json
import re
import sys
from pathlib import Path


DEFAULT_PATTERNS = {
    "anat_patterns": ["*T1w*.nii.gz", "*T1w*.nii"],
    "func_patterns": ["*bold*.nii.gz", "*bold*.nii"],
    "dwi_patterns": ["*dwi*.nii.gz", "*dwi*.nii"],
    "func_reverse_pe_patterns": [],
    "dwi_reverse_pe_patterns": [],
    "func_exclude_patterns": ["*sbref*", "*acq-rev*", "*dir-PA*", "*dir-flipped*"],
    "reverse_pe_patterns": [
        "*acq-rev*bold*.nii.gz",
        "*acq-rev*bold*.nii",
        "*dir-PA*bold*.nii.gz",
        "*dir-PA*bold*.nii",
        "*dir-flipped*bold*.nii.gz",
        "*dir-flipped*bold*.nii",
    ],
    "fieldmap_magnitude_patterns": [
        "*magnitude1*.nii.gz",
        "*magnitude1*.nii",
        "*magnitude*.nii.gz",
        "*magnitude*.nii",
    ],
    "fieldmap_phasediff_patterns": ["*phasediff*.nii.gz", "*phasediff*.nii"],
}

UNDIST_CHOICES = {"none", "auto", "blip", "fieldmap"}
TASK_MODE_CHOICES = {"none", "auto", "require"}

DEFAULT_ARGS = {
    "undist": "none",
    "tr_msec": None,
    "drop": "[3,0]",
    "task_mode": "none",
    "inspect": False,
}


def apply_defaults(args):
    for key, value in DEFAULT_ARGS.items():
        if getattr(args, key) is None:
            setattr(args, key, value)

    if args.undist not in UNDIST_CHOICES:
        valid = ", ".join(sorted(UNDIST_CHOICES))
        sys.exit(f"ERROR: undist must be one of: {valid}")
    if args.task_mode not in TASK_MODE_CHOICES:
        valid = ", ".join(sorted(TASK_MODE_CHOICES))
        sys.exit(f"ERROR: task_mode must be one of: {valid}")

    return args


def load_text_patterns(pattern_file, patterns):
    section = None
    seen_sections = set()

    with open(pattern_file, encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip()
                if section not in patterns and not is_metadata_section(section):
                    valid = ", ".join(sorted(patterns)) + ", [defaults], [manufacturer:NAME], [site:NAME], [site:NAME manufacturer:NAME]"
                    sys.exit(f"ERROR: unknown pattern section [{section}] on line {line_number}. Valid sections: {valid}")
                if section not in seen_sections:
                    patterns[section] = {} if is_metadata_section(section) else []
                    seen_sections.add(section)
                continue

            if section is None:
                sys.exit(f"ERROR: pattern line outside a section on line {line_number}: {line}")

            if is_metadata_section(section):
                if section.lower().startswith("dwi "):
                    continue
                if "=" not in line:
                    sys.exit(f"ERROR: expected key=value in [{section}] on line {line_number}")
                key, value = line.split("=", 1)
                key = key.strip().lower()
                if key not in {"tpattern", "slice_axis", "tr_msec"}:
                    sys.exit(f"ERROR: unknown metadata rule field {key!r} on line {line_number}")
                patterns[section][key] = value.strip()
            else:
                patterns[section].append(line)


def is_metadata_section(section):
    return (
        section.lower() == "defaults"
        or re.fullmatch(r"(?:bold|dwi) manufacturer:.+", section, re.IGNORECASE) is not None
        or re.fullmatch(r"(?:bold|dwi) site:.+", section, re.IGNORECASE) is not None
        or re.fullmatch(r"(?:bold|dwi) site:.+ manufacturer:.+", section, re.IGNORECASE) is not None
        or re.fullmatch(r"manufacturer:.+", section, re.IGNORECASE) is not None
        or re.fullmatch(r"site:.+", section, re.IGNORECASE) is not None
        or re.fullmatch(r"site:.+ manufacturer:.+", section, re.IGNORECASE) is not None
    )


def load_patterns(pattern_file):
    patterns = {key: list(value) for key, value in DEFAULT_PATTERNS.items()}
    patterns["metadata_rules"] = {}

    if pattern_file is None:
        return patterns

    load_text_patterns(pattern_file, patterns)

    for section in list(patterns):
        if is_metadata_section(section):
            patterns["metadata_rules"][section] = patterns.pop(section)

    return patterns


def extend_patterns(patterns, key, values):
    if values:
        patterns[key] = values


def matches_any(path, patterns):
    name = path.name
    return any(fnmatch.fnmatchcase(name, pattern) for pattern in patterns)


def matching_files(folder, include_patterns, exclude_patterns=None):
    exclude_patterns = exclude_patterns or []
    files = []

    if not folder.is_dir():
        return files

    for pattern in include_patterns:
        for path in sorted(folder.glob(pattern)):
            if path.name.startswith(".") or path in files:
                continue
            if matches_any(path, exclude_patterns):
                continue
            files.append(path)

    return files


def first_matching_file(folder, include_patterns):
    files = matching_files(folder, include_patterns)
    return files[0] if files else None


def find_sessions(subject_folder):
    sessions = sorted(path for path in subject_folder.glob("ses-*") if path.is_dir())

    if sessions:
        return sessions

    if (subject_folder / "anat").is_dir() or (subject_folder / "func").is_dir():
        return [subject_folder]

    return []


def label_for(subject_folder, session_folder):
    if session_folder == subject_folder:
        return subject_folder.name

    return f"{subject_folder.name}_{session_folder.name}"


def session_label(subject_folder, session_folder):
    if session_folder == subject_folder:
        return ""

    return session_folder.name


def find_anat(subject_folder, session_folder, patterns):
    anat_folder = session_folder / "anat"
    anat = first_matching_file(anat_folder, patterns["anat_patterns"])

    if anat is not None:
        return anat

    if session_folder != subject_folder:
        return first_matching_file(subject_folder / "anat", patterns["anat_patterns"])

    return None


def find_func(session_folder, patterns):
    return matching_files(
        session_folder / "func",
        patterns["func_patterns"],
        patterns["func_exclude_patterns"],
    )


def find_reverse_pe(session_folder, patterns):
    reverse_patterns = patterns.get("func_reverse_pe_patterns") or patterns["reverse_pe_patterns"]
    return first_matching_file(session_folder / "func", reverse_patterns)


def find_fieldmap(session_folder, patterns):
    fmap_folder = session_folder / "fmap"
    magnitude = first_matching_file(fmap_folder, patterns["fieldmap_magnitude_patterns"])
    phasediff = first_matching_file(fmap_folder, patterns["fieldmap_phasediff_patterns"])

    if magnitude and phasediff:
        return [magnitude, phasediff]

    return []


def event_file_for(func_file):
    name = func_file.name

    if name.endswith(".nii.gz"):
        stem = name[:-7]
    elif name.endswith(".nii"):
        stem = name[:-4]
    else:
        return None

    if not stem.endswith("_bold"):
        return None

    event_file = func_file.with_name(stem[:-5] + "_events.tsv")
    return event_file if event_file.is_file() else None


def image_stem(image_file):
    name = image_file.name

    if name.endswith(".nii.gz"):
        return name[:-7]
    if name.endswith(".nii"):
        return name[:-4]

    return image_file.stem


def timing_file_name(func_file):
    stem = image_stem(func_file)

    if stem.endswith("_bold"):
        stem = stem[:-5]

    return stem + "_slicetime.txt"


def find_tpattern(tpattern_path, func_file):
    if not tpattern_path.is_dir():
        return None

    tpattern_file = tpattern_path / timing_file_name(func_file)
    return tpattern_file if tpattern_file.is_file() else None


def json_file_for(func_file):
    if func_file.name.endswith(".nii.gz"):
        return func_file.with_name(func_file.name[:-7] + ".json")
    return func_file.with_suffix(".json")


def load_sidecar(func_file):
    json_path = json_file_for(func_file)
    if not json_path.is_file():
        return json_path, {}
    try:
        with json_path.open(encoding="utf-8") as file:
            data = json.load(file)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read JSON sidecar {json_path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"JSON sidecar must contain an object: {json_path}")
    return json_path, data


def nifti_shape(func_file):
    try:
        import nibabel as nib
    except ImportError as error:
        raise ValueError(f"cannot determine NIfTI shape for {func_file}; install nibabel") from error
    try:
        return tuple(nib.load(str(func_file)).shape)
    except Exception as error:
        raise ValueError(f"cannot read NIfTI shape for {func_file}: {error}") from error


def normalize_manufacturer(value):
    value = str(value or "").strip()
    lowered = value.lower()
    if "siemens" in lowered:
        return "Siemens"
    if "philips" in lowered:
        return "Philips"
    if lowered in {"ge", "ge medical systems"} or "general electric" in lowered:
        return "GE"
    return value


def metadata_value(data, *keys):
    for key in keys:
        value = data.get(key)
        if value not in (None, ""):
            return str(value).strip()
    return ""


def participants_sites(bids_root):
    path = bids_root / "participants.tsv"
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as file:
        rows = [line.rstrip("\n").split("\t") for line in file]
    if not rows:
        return {}
    header = {name: index for index, name in enumerate(rows[0])}
    if "participant_id" not in header or "site" not in header:
        return {}
    result = {}
    for row in rows[1:]:
        if len(row) > max(header["participant_id"], header["site"]):
            result[row[header["participant_id"]].strip()] = row[header["site"]].strip()
    return result


def slice_axis(direction):
    match = re.fullmatch(r"([ijk])(?:[-+])?", str(direction or ""), re.IGNORECASE)
    if not match:
        return None
    return {"i": 0, "j": 1, "k": 2}[match.group(1).lower()]


def validate_slice_timing(func_file, json_path, data):
    values = data.get("SliceTiming")
    if values is None:
        return None, None
    if not isinstance(values, list) or not values:
        raise ValueError(f"invalid SliceTiming in {json_path}: expected a nonempty numeric list")
    try:
        values = [float(value) for value in values]
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid SliceTiming in {json_path}: expected a nonempty numeric list") from error
    if any(value < 0 for value in values):
        raise ValueError(f"invalid SliceTiming in {json_path}: values cannot be negative")
    direction = data.get("SliceEncodingDirection")
    shape = None
    if direction is not None:
        axis = slice_axis(direction)
        if axis is None:
            raise ValueError(f"invalid SliceEncodingDirection {direction!r} in {json_path}")
        shape = nifti_shape(func_file)
        if axis >= len(shape) or len(values) != shape[axis]:
            raise ValueError(
                f"SliceTiming length mismatch: BOLD {func_file}; JSON {json_path}; timing length {len(values)}; "
                f"selected axis {axis} ({direction}); NIfTI shape {shape}"
            )
    return values, shape


def write_timing_file(func_file, values, metadata_dir):
    metadata_dir.mkdir(parents=True, exist_ok=True)
    path = metadata_dir / timing_file_name(func_file)
    content = "\n".join(str(value).rstrip("0").rstrip(".") if value % 1 else str(int(value)) for value in values) + "\n"
    if not path.is_file() or path.read_text(encoding="utf-8") != content:
        path.write_text(content, encoding="utf-8")
    return path


def parse_rule_section(section):
    scope = "bold"
    scoped = re.fullmatch(r"(bold|dwi)\s+(.+)", section, re.IGNORECASE)
    if scoped:
        scope = scoped.group(1).lower()
        section = scoped.group(2).strip()
    if section.lower() == "defaults":
        return "defaults", "", "", scope
    combined = re.fullmatch(r"site:(.+) manufacturer:(.+)", section, re.IGNORECASE)
    if combined:
        return "combined", combined.group(1).strip(), normalize_manufacturer(combined.group(2)), scope
    manufacturer = re.fullmatch(r"manufacturer:(.+)", section, re.IGNORECASE)
    if manufacturer:
        return "manufacturer", "", normalize_manufacturer(manufacturer.group(1)), scope
    site = re.fullmatch(r"site:(.+)", section, re.IGNORECASE)
    if site:
        return "site", site.group(1).strip(), "", scope
    return None, "", "", scope


def rule_matches(section, site, manufacturer):
    kind, rule_site, rule_manufacturer, _ = parse_rule_section(section)
    if kind == "defaults":
        return True
    if kind == "manufacturer":
        return bool(manufacturer) and manufacturer.casefold() == rule_manufacturer.casefold()
    if kind == "site":
        return bool(site) and site.casefold() == rule_site.casefold()
    if kind == "combined":
        return bool(site and manufacturer) and site.casefold() == rule_site.casefold() and manufacturer.casefold() == rule_manufacturer.casefold()
    return False


def resolve_rule_values(rules, site, manufacturer, scope="bold"):
    resolved = {}
    sources = {}
    specificity = {"defaults": 0, "manufacturer": 1, "site": 2, "combined": 3}
    for section, values in rules.items():
        kind, _, _, rule_scope = parse_rule_section(section)
        if kind is None or (rule_scope not in {scope, "common"}) or not rule_matches(section, site, manufacturer):
            continue
        for field, value in values.items():
            if field in resolved and specificity[kind] == specificity[sources[field][0]] and resolved[field] != value:
                raise ValueError(f"conflicting equal-specificity rules for {field}: {sources[field][1]} and {section}")
            if field not in resolved or specificity[kind] >= specificity[sources[field][0]]:
                resolved[field] = value
                sources[field] = (kind, section)
    return resolved, {field: f"rule:{section}" for field, (_, section) in sources.items()}


def parse_tr_msec(value, source):
    if value is None:
        return None
    match = re.fullmatch(r"\[?\s*([0-9]+(?:\.[0-9]+)?)\s*\]?", str(value).strip())
    if not match or float(match.group(1)) <= 0:
        raise ValueError(f"invalid TR_MSEC value {value!r} from {source}")
    number = float(match.group(1))
    return f"[{int(number) if number.is_integer() else number:g}]"


def parse_slice_axis(value, source):
    axis = slice_axis(value)
    if axis is None:
        raise ValueError(f"invalid slice_axis value {value!r} from {source}; use i, j, or k")
    return axis


def resolve_run_metadata(func_file, args, patterns, bids_root, metadata_dir, participants):
    json_path, data = load_sidecar(func_file)
    manufacturer = normalize_manufacturer(data.get("Manufacturer"))
    participant = next((part for part in func_file.parts if part.startswith("sub-")), "")
    site = participants.get(participant, "") or metadata_value(data, "SiteName", "InstitutionName", "StationName")
    rules, rule_sources = resolve_rule_values(patterns["metadata_rules"], site, manufacturer, scope="bold")
    tpattern = args.tpattern
    tr_msec = args.tr_msec
    tpattern_source = "explicit" if tpattern is not None else None
    tr_source = "explicit" if tr_msec is not None else None
    slice_timing, shape = validate_slice_timing(func_file, json_path, data)
    shape = shape or nifti_shape(func_file)

    if tpattern is None and slice_timing is not None:
        tpattern = write_timing_file(func_file, slice_timing, metadata_dir)
        tpattern_source = "json:SliceTiming"
    if tr_msec is None and data.get("RepetitionTime") is not None:
        try:
            seconds = float(data["RepetitionTime"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid RepetitionTime in {json_path}: expected a positive number") from error
        if seconds <= 0:
            raise ValueError(f"invalid RepetitionTime in {json_path}: expected a positive number")
        tr_msec = f"[{seconds * 1000:g}]"
        tr_source = "json:RepetitionTime"
    if tpattern is None:
        tpattern = rules.get("tpattern")
        tpattern_source = rule_sources.get("tpattern")
    if tr_msec is None:
        tr_msec = rules.get("tr_msec")
        tr_source = rule_sources.get("tr_msec")
    if tpattern is None:
        raise ValueError(f"unresolved TPATTERN for BOLD {func_file} (JSON {json_path})")
    if tr_msec is None:
        raise ValueError(f"unresolved TR_MSEC for BOLD {func_file} (JSON {json_path})")
    tr_msec = parse_tr_msec(tr_msec, tr_source)

    tpattern_path = Path(str(tpattern)).expanduser()
    if tpattern_source == "explicit" and tpattern_path.is_dir():
        tpattern_file = find_tpattern(tpattern_path, func_file)
        if tpattern_file is None:
            raise ValueError(f"missing per-run TPATTERN file for {func_file} in {tpattern_path}")
        tpattern = tpattern_file
    elif str(tpattern).casefold() == "parity_alt+z":
        slice_count = len(slice_timing) if slice_timing is not None else None
        if slice_count is None:
            axis = slice_axis(data.get("SliceEncodingDirection"))
            if axis is None and rules.get("slice_axis") is not None:
                axis = parse_slice_axis(rules["slice_axis"], rule_sources.get("slice_axis", "rule"))
            if axis is None:
                raise ValueError(f"cannot resolve parity_alt+z for {func_file}: missing valid SliceEncodingDirection and SliceTiming")
            shape = shape or nifti_shape(func_file)
            if axis >= len(shape):
                raise ValueError(f"slice_axis {axis} is outside NIfTI shape {shape} for {func_file}")
            slice_count = shape[axis]
        tpattern = "alt+z" if slice_count % 2 else "alt+z2"
    return {
        "tpattern": str(tpattern),
        "tpattern_source": tpattern_source or rule_sources.get("tpattern", "defaults"),
        "tr_msec": tr_msec,
        "tr_source": tr_source or rule_sources.get("tr_msec", "defaults"),
        "json": json_path,
        "manufacturer": manufacturer or "unknown",
        "site": site or "unknown",
        "shape": shape,
        "slice_count": len(slice_timing) if slice_timing is not None else None,
    }


def prefix_for_func(func_file):
    stem = image_stem(func_file)

    if stem.endswith("_bold"):
        stem = stem[:-5]

    return stem


def bracket(path):
    return f"[{path}]" if path else "[]"


def bracket_list(paths):
    return "[" + ",".join(str(path) for path in paths) + "]" if paths else "[]"


def row_for_record(record, args):
    funcs = record["func"]
    drop_text = ",".join([args.drop] * len(funcs))
    fields = [
        f"PREFIX={record['prefix']}",
        f"ANAT={record['anat']}",
        "FUNC=" + ",".join(str(path) for path in funcs),
        f"TPATTERN={record['tpattern']}",
        f"TR_MSEC={record['tr_msec']}",
        f"DROP={drop_text}",
    ]

    if args.seed:
        fields.append(f"SEED={args.seed}")

    if args.task_mode != "none":
        task_files = [event_file_for(path) for path in funcs]
        if args.task_mode == "require" and any(path is None for path in task_files):
            missing = [str(funcs[i]) for i, path in enumerate(task_files) if path is None]
            sys.exit("ERROR: missing events files for:\n  " + "\n  ".join(missing))
        if all(path is not None for path in task_files):
            fields.append("TASK=" + ",".join(str(path) for path in task_files))

    if args.undist != "none":
        fields.append(f"DIST={bracket(args.dist_file)}")
        fields.append(f"PE_rev={bracket(record['reverse_pe'])}")
        fields.append(f"FIELD={bracket_list(record['fieldmap'])}")

    return " ".join(fields)


def inspect_record(record):
    lines = [
        f"{record['prefix']}:",
        f"  ANAT: {record['anat'] or 'missing'}",
        f"  FUNC: {len(record['func'])}",
    ]
    for index, path in enumerate(record["func"]):
        event = event_file_for(path)
        suffix = f" events={event}" if event else ""
        lines.append(f"    - {path}{suffix}")
        metadata = record["metadata"][index]
        lines.extend([
            f"      BOLD: {path}",
            f"      JSON: {metadata['json']}",
            f"      Manufacturer: {metadata['manufacturer']}",
            f"      Site: {metadata['site']}",
            f"      NIfTI shape: {metadata['shape'] or 'unknown'}",
            f"      Slice count: {metadata['slice_count'] or 'unknown'}",
            f"      TPATTERN: {metadata['tpattern']}",
            f"      TPATTERN source: {metadata['tpattern_source']}",
            f"      TR_MSEC: {metadata['tr_msec']}",
            f"      TR_MSEC source: {metadata['tr_source']}",
        ])
    lines.append(f"  PE_rev: {record['reverse_pe'] or 'missing'}")
    fieldmap = ", ".join(str(path) for path in record["fieldmap"]) or "missing"
    lines.append(f"  FIELD: {fieldmap}")
    lines.append(f"  TPATTERN: {record['tpattern']}")
    lines.append(f"  TR_MSEC: {record['tr_msec']}")
    return "\n".join(lines)


def build_records(args, patterns):
    bids_root = Path(args.bids_root).expanduser().resolve()
    metadata_dir = Path(args.output).expanduser().with_name(Path(args.output).stem + "_metadata")
    participants = participants_sites(bids_root)
    records = []
    skipped = []
    errors = []

    if not bids_root.is_dir():
        sys.exit(f"ERROR: missing BIDS root: {bids_root}")

    subject_folders = [
        path for path in sorted(bids_root.glob("sub-*"))
        if path.is_dir() and (not args.subject or args.subject in path.name)
    ]
    total_subjects = len(subject_folders)
    processed_subjects = 0
    processed_sessions = 0

    for subject_index, subject_folder in enumerate(subject_folders, start=1):
        processed_subjects += 1
        print(f"[{subject_index}/{total_subjects}] Processing {subject_folder.name}", flush=True)

        for session_folder in find_sessions(subject_folder):
            ses_label = session_label(subject_folder, session_folder)
            if args.session and args.session not in ses_label:
                continue
            processed_sessions += 1

            anat = find_anat(subject_folder, session_folder, patterns)
            funcs = find_func(session_folder, patterns)
            reverse_pe = find_reverse_pe(session_folder, patterns)
            fieldmap = find_fieldmap(session_folder, patterns)
            prefix = label_for(subject_folder, session_folder)

            reasons = []
            if anat is None:
                reasons.append("missing anat")
            if not funcs:
                reasons.append("missing func")
            if args.undist == "blip" and reverse_pe is None:
                reasons.append("missing reverse PE")
            if args.undist == "fieldmap" and not fieldmap:
                reasons.append("missing fieldmap")
            if args.undist in {"auto", "blip", "fieldmap"} and not args.dist_file:
                reasons.append("missing --dist-file")
            if args.undist == "auto" and reverse_pe is None and not fieldmap:
                reasons.append("missing reverse PE or fieldmap")

            record = {
                "prefix": prefix,
                "anat": anat,
                "func": funcs,
                "reverse_pe": reverse_pe,
                "fieldmap": fieldmap,
            }

            if reasons:
                skipped.append((record, reasons))
                print(
                    f"[{subject_index}/{total_subjects}] SKIPPING {prefix}: {', '.join(reasons)}",
                    flush=True,
                )
            else:
                metadata = []
                for func in funcs:
                    try:
                        metadata.append(resolve_run_metadata(func, args, patterns, bids_root, metadata_dir, participants))
                    except ValueError as error:
                        message = str(error)
                        errors.append(message)
                        print(
                            f"[{subject_index}/{total_subjects}] ERROR for {prefix}: {message}",
                            flush=True,
                        )
                if len(metadata) != len(funcs):
                    skipped.append((record, ["metadata resolution failed"]))
                    print(
                        f"[{subject_index}/{total_subjects}] SKIPPING {prefix}: "
                        "one or more functional runs failed metadata resolution",
                        flush=True,
                    )
                    continue
                record["metadata"] = metadata
                first = metadata[0]
                same_metadata = all(
                    item["tpattern"] == first["tpattern"] and item["tr_msec"] == first["tr_msec"]
                    for item in metadata
                )
                if same_metadata:
                    record.update({
                        "tpattern": first["tpattern"],
                        "tr_msec": first["tr_msec"],
                    })
                    records.append(record)
                else:
                    for func, item in zip(funcs, metadata):
                        records.append({
                            **record,
                            "prefix": prefix_for_func(func),
                            "func": [func],
                            "metadata": [item],
                            "tpattern": item["tpattern"],
                            "tr_msec": item["tr_msec"],
                        })

    print(
        f"Discovery summary: subjects={processed_subjects} sessions={processed_sessions} "
        f"rows={len(records)} skipped={len(skipped)} errors={len(errors)}",
        flush=True,
    )
    return records, skipped, errors


def parse_args():
    parser = argparse.ArgumentParser(description="Generate OPPNI-B input rows from BIDS-like data")
    required = parser.add_argument_group("required")
    required.add_argument("--bids-root", required=True, metavar="PATH", help="BIDS dataset root")
    required.add_argument("--output", required=True, metavar="FILE", help="Output OPPNI input file")
    optional = parser.add_argument_group("optional")
    optional.add_argument("--patterns", metavar="FILE", help="Dataset-specific filename and metadata rules")
    optional.add_argument("--func-pattern", action="append", help="Override FUNC include pattern; repeatable")
    optional.add_argument("--exclude-func-pattern", action="append", help="Override FUNC exclude pattern; repeatable")
    optional.add_argument("--reverse-pe-pattern", action="append", help="Override reverse-PE pattern; repeatable")
    optional.add_argument("--subject", help="Only include subjects whose folder name contains this text")
    optional.add_argument("--session", help="Only include sessions whose folder name contains this text")
    optional.add_argument("--undist", choices=sorted(UNDIST_CHOICES), default=None, metavar="MODE", help="Distortion-correction mode")
    optional.add_argument("--dist-file", help="OPPNI distortion parameter file")
    optional.add_argument("--tpattern", metavar="VALUE", help="Override TPATTERN with an AFNI code, shared timing file, or timing folder")
    optional.add_argument("--tr-msec", metavar="VALUE", help="Override TR_MSEC in milliseconds")
    optional.add_argument("--drop", help="OPPNI DROP value per functional run")
    optional.add_argument("--seed", help="Optional OPPNI SEED file")
    optional.add_argument("--task-mode", choices=sorted(TASK_MODE_CHOICES), default=None, metavar="MODE", help="Task-event handling mode")
    optional.add_argument("--inspect", action="store_true", help="Print resolved files, metadata, and value sources")
    return parser.parse_args()


def main():
    args = apply_defaults(parse_args())
    patterns = load_patterns(args.patterns)
    extend_patterns(patterns, "func_patterns", args.func_pattern)
    extend_patterns(patterns, "func_exclude_patterns", args.exclude_func_pattern)
    extend_patterns(patterns, "reverse_pe_patterns", args.reverse_pe_pattern)

    records, skipped, errors = build_records(args, patterns)

    if errors:
        sys.exit(1)

    if args.inspect:
        for record in records:
            print(inspect_record(record))
        for record, reasons in skipped:
            print(f"{record['prefix']}: skipped ({', '.join(reasons)})")

    output = Path(args.output).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(row_for_record(record, args) for record in records) + ("\n" if records else ""))

    print(f"Wrote {len(records)} row(s) to {output}")
    if skipped:
        print(f"Skipped {len(skipped)} record(s)")


if __name__ == "__main__":
    main()
