#!/usr/bin/env python3
"""Convert MINC2/HDF5 source files into a configured BIDS dataset."""

import argparse
import configparser
import csv
import json
import re
import shutil
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path


class ConversionError(Exception):
    """A source can be skipped without stopping the remaining batch."""


@dataclass
class MincHeader:
    shape: tuple
    affine: object


class Tee:
    """Write the complete run stream to the terminal and a log file."""

    def __init__(self, console, log_file):
        self.console = console
        self.log_file = log_file

    def write(self, text):
        self.console.write(text)
        self.log_file.write(text)
        self.log_file.flush()

    def flush(self):
        self.console.flush()
        self.log_file.flush()


def text_value(value):
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip(" _")
    if isinstance(value, list):
        return [text_value(item) for item in value]
    return value.item() if hasattr(value, "item") else value


def numeric(value):
    value = text_value(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return float(value) if isinstance(value, (int, float)) else None


def numeric_array(value):
    value = text_value(value)
    if value is None:
        return None
    if not isinstance(value, list):
        value = [value]
    result = [numeric(item) for item in value]
    return result if all(item is not None for item in result) else None


def _read_minc_source(source, load_data=True):
    """Read MINC2 data without relying on nibabel's regular-time check."""
    try:
        import h5py
        import numpy as np
        import nibabel as nib
    except ModuleNotFoundError as error:
        raise ConversionError(
            "MINC2 conversion requires h5py, numpy, and nibabel in the active environment"
        ) from error

    with h5py.File(source, "r") as file:
        root = file.get("/minc-2.0")
        if root is None or "/minc-2.0/image/0/image" not in file:
            raise ConversionError(f"not a supported MINC2/HDF5 file: {source}")

        dataset = file["/minc-2.0/image/0/image"]
        source_shape = tuple(dataset.shape)
        dimorder = text_value(dataset.attrs.get("dimorder"))
        if isinstance(dimorder, str):
            dimorder = [item.strip() for item in dimorder.split(",") if item.strip()]
        if not dimorder or len(dimorder) != len(dataset.shape):
            raise ConversionError(f"missing or invalid MINC dimension order: {source}")

        allowed = {"xspace", "yspace", "zspace", "time"}
        if any(axis not in allowed for axis in dimorder):
            raise ConversionError(f"unsupported MINC dimensions {dimorder}: {source}")
        required = {"xspace", "yspace", "zspace"}
        if not required.issubset(dimorder):
            raise ConversionError(f"MINC file lacks xspace, yspace, or zspace: {source}")

        raw = np.asarray(dataset[...]) if load_data else None
        target_axes = ["xspace", "yspace", "zspace"]
        if "time" in dimorder:
            target_axes.append("time")
        transpose_order = [dimorder.index(axis) for axis in target_axes]
        data = None
        if load_data:
            data = np.transpose(raw, transpose_order) if transpose_order != list(range(raw.ndim)) else raw

        image_group = file["/minc-2.0/image/0"]
        valid_range = numeric_array(dataset.attrs.get("valid_range")) or [0.0, 4095.0]
        image_min = image_group.get("image-min")
        image_max = image_group.get("image-max")
        if load_data and image_min is not None and image_max is not None and np.issubdtype(raw.dtype, np.integer):
            # MINC stores integer voxels plus per-slice/volume real-value bounds.
            mins = np.asarray(image_min[...])
            maxs = np.asarray(image_max[...])
            raw_float = raw.astype(np.float32)
            scale = (maxs - mins) / (valid_range[1] - valid_range[0])
            offset = mins - valid_range[0] * scale
            scale_dimorder = text_value(image_min.attrs.get("dimorder"))
            if isinstance(scale_dimorder, str):
                scale_dimorder = [item.strip() for item in scale_dimorder.split(",") if item.strip()]
            scale_shape = [1] * raw.ndim
            for axis in scale_dimorder or []:
                if axis in dimorder:
                    scale_shape[dimorder.index(axis)] = dataset.shape[dimorder.index(axis)]
            scale = scale.reshape(scale_shape)
            offset = offset.reshape(scale_shape)
            raw_float = raw_float * scale + offset
            data = np.transpose(raw_float, transpose_order) if transpose_order != list(range(raw.ndim)) else raw_float

        dimensions = file["/minc-2.0/dimensions"]
        affine = np.eye(4, dtype=float)
        origin = np.zeros(3, dtype=float)
        for column, axis in enumerate(("xspace", "yspace", "zspace")):
            attrs = dimensions[axis].attrs
            step = numeric(attrs.get("step"))
            start = numeric(attrs.get("start"))
            direction = numeric_array(attrs.get("direction_cosines"))
            if step is None or start is None or direction is None or len(direction) != 3:
                raise ConversionError(f"incomplete MINC spatial metadata for {axis}: {source}")
            affine[:3, column] = np.asarray(direction) * step
            origin += np.asarray(direction) * start
        affine[:3, 3] = origin

        metadata, diffusion = read_minc_metadata(file)
        if "time" in dimorder:
            time_step, time_source = read_trusted_time_step(file, dimensions, source)
            metadata["RepetitionTime"] = float(time_step)
        else:
            time_step = None
            time_source = None

    shape = tuple(source_shape[index] for index in transpose_order)
    zooms = tuple(float(np.linalg.norm(affine[:3, i])) for i in range(3))
    zooms = zooms + ((float(time_step),) if time_step is not None and len(shape) == 4 else ())
    if not load_data:
        return MincHeader(shape, affine), metadata, diffusion, time_source
    image = nib.Nifti1Image(data, affine)
    image.header.set_zooms(zooms)
    image.header.set_xyzt_units("mm", "sec" if len(shape) == 4 else None)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    return image, metadata, diffusion, time_source


def read_minc_source(source, load_data=True):
    try:
        return _read_minc_source(source, load_data=load_data)
    except ConversionError:
        raise
    except (KeyError, OSError, RuntimeError, TypeError, ValueError) as error:
        raise ConversionError(f"malformed or incomplete MINC metadata: {source}: {error}") from error


def read_nifti_source(source):
    try:
        import nibabel as nib
    except ModuleNotFoundError as error:
        raise ConversionError("NIfTI organization requires nibabel in the active environment") from error
    try:
        image = nib.load(str(source))
        return MincHeader(tuple(image.shape), image.affine), {}, {}, None
    except (OSError, ValueError, RuntimeError) as error:
        raise ConversionError(f"could not read NIfTI header: {source}: {error}") from error


def read_trusted_time_step(file, dimensions, source):
    """Return a positive time step only when MINC units are defensible."""
    time_attrs = dimensions["time"].attrs
    dimension_step = numeric(time_attrs.get("step"))
    dimension_units = str(text_value(time_attrs.get("units", ""))).strip().lower()
    unit_factors = {"s": 1.0, "sec": 1.0, "secs": 1.0, "second": 1.0, "seconds": 1.0, "ms": 0.001, "msec": 0.001, "millisecond": 0.001, "milliseconds": 0.001}
    info = file["/minc-2.0/info"]
    acquisition_group = info.get("acquisition")
    acquisition = {} if acquisition_group is None else {key: text_value(value) for key, value in acquisition_group.attrs.items()}
    repetition_time = numeric(acquisition.get("repetition_time"))
    dicom_group = info.get("dicom_0x0018")
    dicom_tr = None if dicom_group is None else numeric(text_value(dicom_group.attrs.get("el_0x0080")))
    # COMPASS stores acquisition repetition_time in seconds; confirm it against
    # DICOM 0018,0080, which is stored in milliseconds.
    if repetition_time is not None and repetition_time > 0 and dicom_tr is not None and abs(repetition_time - dicom_tr / 1000.0) < 1e-3:
        return repetition_time, "MINC acquisition repetition_time confirmed against DICOM milliseconds"
    if dimension_step is not None and dimension_step > 0 and dimension_units in unit_factors:
        return dimension_step * unit_factors[dimension_units], f"MINC time dimension step ({dimension_units})"
    raise ConversionError(f"no trustworthy positive TR/time step with known units: {source}")


def read_minc_metadata(file):
    info = file["/minc-2.0/info"]

    def attrs(group_name):
        group = info.get(group_name)
        return {} if group is None else {key: text_value(value) for key, value in group.attrs.items()}

    acquisition = attrs("acquisition")

    def dicom(tag):
        group = attrs(f"dicom_0x{tag[:4]}")
        return group.get(f"el_0x{tag[4:]}")

    metadata = {}

    def add(name, value):
        if value not in (None, ""):
            metadata[name] = value

    add("Modality", "MR")
    add("MagneticFieldStrength", numeric(dicom("00180087")))
    add("Manufacturer", dicom("00080070"))
    add("ManufacturersModelName", dicom("00081090"))
    add("InstitutionName", dicom("00080080"))
    add("StationName", dicom("00081010"))
    add("SoftwareVersions", dicom("00181020"))
    add("SeriesDescription", dicom("0008103E"))
    add("ProtocolName", dicom("00181030"))
    add("ScanningSequence", dicom("00180020"))
    add("SequenceVariant", dicom("00180021"))
    add("ScanOptions", dicom("00180022"))
    add("MRAcquisitionType", dicom("00180023"))
    add("PulseSequenceName", dicom("00180024"))
    add("SeriesNumber", numeric(dicom("00200011")))
    acquisition_time = dicom("00080032") or acquisition.get("acquisition_time")
    if acquisition_time not in (None, ""):
        parsed_time = format_acquisition_time(acquisition_time)
        if parsed_time is not None:
            add("AcquisitionTime", parsed_time)
    add("SliceThickness", numeric(dicom("00180050")))
    add("SpacingBetweenSlices", numeric(dicom("00180088")))
    add("FlipAngle", numeric(dicom("00181314")))
    add("NumberOfAverages", numeric(dicom("00180083")))
    add("EchoTime", numeric(acquisition.get("echo_time")))
    add("RepetitionTime", numeric(acquisition.get("repetition_time")))
    add("InversionTime", numeric(acquisition.get("inversion_time")))

    diffusion = {key: numeric_array(acquisition.get(key)) for key in ("bvalues", "direction_x", "direction_y", "direction_z")}
    return metadata, diffusion


def format_acquisition_time(value):
    """Convert DICOM HHMMSS.frac to the BIDS HH:MM:SS.frac form."""
    value = str(text_value(value)).strip()
    match = re.fullmatch(r"(\d{2})(\d{2})(\d{2})(?:\.(\d+))?", value)
    if not match:
        warnings.warn(f"omitting unsafe AcquisitionTime value: {value!r}", RuntimeWarning)
        return None
    hour, minute, second = (int(match.group(index)) for index in (1, 2, 3))
    if hour > 23 or minute > 59 or second > 59:
        warnings.warn(f"omitting unsafe AcquisitionTime value: {value!r}", RuntimeWarning)
        return None
    fraction = match.group(4)
    return f"{hour:02d}:{minute:02d}:{second:02d}" + (f".{fraction}" if fraction else "")


def validate_dwi(image, diffusion, source):
    import numpy as np

    volumes = image.shape[3] if len(image.shape) == 4 else 1
    bvalues = diffusion.get("bvalues")
    directions = [diffusion.get(f"direction_{axis}") for axis in ("x", "y", "z")]
    if bvalues is None or any(direction is None for direction in directions):
        raise ConversionError(f"invalid DWI metadata: missing bvalues or diffusion directions: {source.name}")
    if len(bvalues) != volumes or any(len(direction) != volumes for direction in directions):
        raise ConversionError(
            f"invalid DWI metadata: {volumes} volumes but bvalues/directions have lengths "
            f"{len(bvalues)}, {[len(direction) for direction in directions]}: {source.name}"
        )

    # COMPASS direction_x/y/z are already final FSL/NIfTI image-axis vectors.
    # This is verified against the trusted existing 100659 conversion test.
    bvalues = np.asarray(bvalues, dtype=float)
    image_bvecs = np.asarray(directions, dtype=float)
    diffusion_weighted = bvalues > 0
    norms = np.linalg.norm(image_bvecs, axis=0)
    invalid = diffusion_weighted & (norms <= 1e-12)
    if np.any(invalid):
        indices = np.flatnonzero(invalid).tolist()
        raise ConversionError(
            f"positive b-values have zero diffusion directions "
            f"at indices {indices}: {source.name}"
        )
    bvecs = np.zeros_like(image_bvecs)
    bvecs[:, diffusion_weighted] = image_bvecs[:, diffusion_weighted] / norms[diffusion_weighted]
    return bvalues, bvecs


def load_config(path):
    parser = configparser.ConfigParser()
    parser.optionxform = str
    if not path.is_file():
        raise ConversionError(f"missing config file: {path}")
    parser.read(path, encoding="utf-8")
    if "dataset" not in parser or "scans" not in parser:
        raise ConversionError("config must contain [dataset] and [scans] sections")
    return parser


def parse_scan_types(parser):
    scan_types = {}
    if "bids_scan_types" not in parser:
        return scan_types
    for target, value in parser.items("bids_scan_types"):
        fields = {}
        for item in value.split(","):
            key, separator, field_value = item.partition("=")
            if separator:
                fields[key.strip()] = field_value.strip()
        scan_types[target] = fields
    return scan_types


def filename_regex(template, scan_names, extension, allow_trailing=False):
    parts = re.findall(r"\{([^}]+)\}|([^{}]+)", template)
    pattern = []
    for placeholder, literal in parts:
        if literal:
            pattern.append(re.escape(literal))
        elif placeholder == "scan":
            values = "|".join(re.escape(value) for value in sorted(scan_names, key=len, reverse=True))
            pattern.append(f"(?P<scan>{values})")
        elif placeholder == "session":
            pattern.append(r"(?P<session>.+?)")
        else:
            pattern.append(fr"(?P<{placeholder}>[^_]+)")
    trailing = r"(?:_[^_]+(?:-[^_]+)*)?" if allow_trailing else ""
    return re.compile("^" + "".join(pattern) + trailing + re.escape(extension) + "$", re.IGNORECASE)


def parse_entities(source, parser, extension, source_format="minc"):
    dataset = parser["dataset"]
    scan_names = list(parser["scans"].keys())
    match = filename_regex(dataset["filename_template"], scan_names, extension, source_format == "nifti").match(source.name)
    if not match:
        raise ConversionError(f"filename does not match configured template: {source.name}")
    entities = match.groupdict()
    configured_scan = next((key for key in parser["scans"] if key.casefold() == entities["scan"].casefold()), None)
    if configured_scan is None:
        raise ConversionError(f"no configured scan rule for: {entities['scan']}")
    scan_target = parser["scans"][configured_scan]
    scan_type = parse_scan_types(parser).get(scan_target, {})
    if not scan_type:
        raise ConversionError(f"no [bids_scan_types] rule for scan target: {scan_target}")
    subject = dataset.get("subject_template", "{subject}").format(**entities)
    session = dataset.get("session_template", "{session}").format(**entities)
    subject = re.sub(r"[^A-Za-z0-9_+-]+", "", subject)
    session = re.sub(r"[^A-Za-z0-9_+-]+", "", session)
    if not subject or not session:
        raise ConversionError(f"could not construct BIDS subject/session from: {source.name}")
    return {"subject": f"sub-{subject}", "session": f"ses-{session}", "scan": entities["scan"], "run": entities.get("run"), **scan_type}


def bids_destination(output_root, entities):
    name = f"{entities['subject']}_{entities['session']}"
    if entities.get("task"):
        name += f"_task-{entities['task']}"
    if entities.get("acq"):
        name += f"_acq-{entities['acq']}"
    if entities.get("run"):
        name += f"_run-{entities['run']}"
    name += f"_{entities['suffix']}.nii.gz"
    return output_root / entities["subject"] / entities["session"] / entities["datatype"] / name


def convert_file(source, destination, entities, overwrite, source_format):
    json_path = destination.with_suffix("").with_suffix(".json")
    bval_path = destination.with_suffix("").with_suffix(".bval")
    bvec_path = destination.with_suffix("").with_suffix(".bvec")
    if destination.exists():
        required = [destination, json_path]
        if entities.get("datatype") == "dwi":
            required.extend((bval_path, bvec_path))
        if not all(path.exists() for path in required):
            missing = ", ".join(str(path) for path in required if not path.exists())
            raise ConversionError(f"incomplete existing output for {source.name}; missing: {missing}")
        if not overwrite:
            print(f"skip existing complete output: {destination}")
            return "skipped"
    if source_format == "nifti":
        image, metadata, diffusion, time_source = read_nifti_source(source)
    else:
        image, metadata, diffusion, time_source = read_minc_source(source)
    bvalues = bvecs = None
    if entities.get("datatype") == "dwi":
        bvalues, bvecs = validate_dwi(image, diffusion, source)

    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata.update({"ConversionSoftware": "source_to_bids", "ConversionSoftwareVersion": "experimental", "SourceFile": source.name})
    if time_source:
        metadata["ConversionTimeStepSource"] = time_source
    if bvecs is not None:
        metadata["BvecCoordinateSystem"] = "FSL radiological voxel convention; COMPASS MINC image-axis directions"
    if entities.get("task"):
        metadata["TaskName"] = entities["task"]
    if source_format == "nifti":
        shutil.copy2(source, destination)
    else:
        import nibabel as nib
        nib.save(image, str(destination))
    json_path.write_text(json.dumps(metadata, indent=2) + "\n")
    if bvalues is not None:
        bval_path.write_text(" ".join(f"{value:g}" for value in bvalues) + "\n")
        bvec_path.write_text("\n".join(" ".join(f"{value:.8g}" for value in row) for row in bvecs) + "\n")
    print(f"converted: {source} -> {destination}")
    return "converted"


def parse_extensions(value):
    extensions = [item.strip().lower() for item in value.split(",") if item.strip()]
    extensions = [item if item.startswith(".") else f".{item}" for item in extensions]
    supported = {".mnc", ".nii", ".nii.gz"}
    unsupported = sorted(set(extensions) - supported)
    if unsupported:
        raise ConversionError(
            f"supported source extensions are .mnc, .nii, and .nii.gz; unsupported: {', '.join(unsupported)}"
        )
    return sorted(set(extensions), key=len, reverse=True)


def source_files(input_path, extensions):
    if input_path.is_file():
        return [input_path]
    files = {
        path
        for extension in extensions
        for path in input_path.rglob(f"*{extension}")
        if path.is_file()
    }
    return sorted(files)


def source_extension(source, extensions):
    source_name = source.name.lower()
    for extension in extensions:
        if source_name.endswith(extension):
            return extension
    raise ConversionError(f"unsupported source extension: {source.name}")


def source_format_for(source, extensions):
    extension = source_extension(source, extensions)
    return "nifti" if extension in {".nii", ".nii.gz"} else "minc"


def read_source_header(source, source_format):
    if source_format == "nifti":
        return read_nifti_source(source)
    return read_minc_source(source, load_data=False)


def write_dataset_description(output_root, config):
    description = {
        "Name": config["dataset"].get("name", "Converted imaging dataset"),
        "BIDSVersion": "1.10.0",
        "DatasetType": "raw",
        "GeneratedBy": [{"Name": "source_to_bids", "Version": "experimental"}],
    }
    (output_root / "dataset_description.json").write_text(json.dumps(description, indent=2) + "\n")


def write_bidsignore(output_root):
    (output_root / ".bidsignore").write_text("/*conversion.log\n/*conversion_manifest.tsv\n\n")


def mark_collisions(plans):
    destinations = {}
    for plan in plans:
        if plan["destination"] is not None:
            destinations.setdefault(plan["destination"], []).append(plan)
    collision_count = 0
    for destination, matching in destinations.items():
        if len(matching) > 1:
            collision_count += 1
            sources = ", ".join(str(plan["source"]) for plan in matching)
            print(f"collision: {destination}\n  sources: {sources}", file=sys.stderr)
            for plan in matching:
                plan["error"] = f"destination collision: {destination}"
    return collision_count


def preflight(files, config, extensions, output_root):
    """Parse, inspect, and validate every source before writing any image."""
    plans = []
    total = len(files)
    print(f"starting preflight: inspecting {total} source file(s)")
    for index, source in enumerate(files, start=1):
        print(f"preflight [{index}/{total}]: inspecting {source}")
        try:
            file_extension = source_extension(source, extensions)
            file_format = source_format_for(source, extensions)
            entities = parse_entities(source, config, file_extension, file_format)
            destination = bids_destination(output_root, entities)
            image, _, diffusion, _ = read_source_header(source, file_format)
            if entities.get("datatype") == "dwi":
                validate_dwi(image, diffusion, source)
            plans.append({"source": source, "entities": entities, "destination": destination, "source_format": file_format, "error": None})
        except (ConversionError, OSError, ValueError) as error:
            plans.append({"source": source, "entities": None, "destination": None, "error": str(error)})

    return plans, mark_collisions(plans)


def write_manifest(path, plans):
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(("source_file", "status", "output_file", "reason"))
        for plan in plans:
            destination = plan.get("manifest_destination", "")
            writer.writerow((str(plan["source"]), plan.get("status", "failed"), destination, plan.get("reason", "")))


def main():
    parser = argparse.ArgumentParser(description="Convert MINC2/HDF5 files or organize existing NIfTI files into BIDS.")
    parser.add_argument("input", help="Source MINC file or folder")
    parser.add_argument("--config", help="Configuration file")
    parser.add_argument("--output", help="BIDS output root; default is input folder/bids")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Parse filenames and inspect MINC headers without writing outputs")
    parser.add_argument("--log", help="Run log path; default is OUTPUT/conversion.log")
    parser.add_argument("--manifest", help="Manifest path; default is OUTPUT/conversion_manifest.tsv")
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        sys.exit(f"ERROR: missing input: {input_path}")
    script_folder = Path(__file__).resolve().parent
    config_path = Path(args.config).expanduser().resolve() if args.config else script_folder / "source_to_bids_config.txt"
    try:
        config = load_config(config_path)
        extensions = parse_extensions(config["dataset"].get("source_extension", ".mnc"))
        if input_path.is_file() and input_path.name.lower().endswith((".nii", ".nii.gz")):
            extensions = [".nii.gz" if input_path.name.lower().endswith(".nii.gz") else ".nii"]
        files = source_files(input_path, extensions)
        if input_path.is_dir() and not files and extensions == [".mnc"]:
            nifti_files = source_files(input_path, [".nii.gz", ".nii"])
            if nifti_files:
                extensions = [".nii.gz", ".nii"]
                files = nifti_files
        if not files:
            raise ConversionError(f"no source files with extension(s) {', '.join(extensions)} found in {input_path}")
    except ConversionError as error:
        sys.exit(f"ERROR: {error}")

    output_root = Path(args.output).expanduser().resolve() if args.output else (input_path / "bids" if input_path.is_dir() else input_path.parent / "bids")
    log_file = None
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)
        log_path = Path(args.log).expanduser().resolve() if args.log else output_root / "conversion.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_file = log_path.open("w", encoding="utf-8")
        sys.stdout = Tee(original_stdout, log_file)
        sys.stderr = Tee(original_stderr, log_file)

    exit_code = 0
    try:
        plans, collision_count = preflight(files, config, extensions, output_root)
        counts = {"converted": 0, "skipped": 0, "failed": 0}
        if not args.dry_run:
            write_dataset_description(output_root, config)
            write_bidsignore(output_root)

        total = len(plans)
        for index, plan in enumerate(plans, start=1):
            source = plan["source"]
            if plan["error"]:
                counts["failed"] += 1
                plan.update(status="failed", reason=plan["error"], manifest_destination="")
                print(f"failed [{index}/{total}]: {source}: {plan['error']}", file=sys.stderr)
                continue
            if args.dry_run:
                print(f"would convert [{index}/{total}]: {source} -> {plan['destination']}")
                continue
            print(f"starting conversion [{index}/{total}]: {source}")
            try:
                status = convert_file(source, plan["destination"], plan["entities"], args.overwrite, plan["source_format"])
                counts[status] += 1
                plan.update(status=status, reason="", manifest_destination=str(plan["destination"]))
                print(f"completed [{index}/{total}]: {source} status={status}")
            except (ConversionError, OSError, ValueError) as error:
                counts["failed"] += 1
                plan.update(status="failed", reason=str(error), manifest_destination="")
                print(f"failed [{index}/{total}]: {source}: {error}", file=sys.stderr)

        if args.dry_run:
            dry_failures = sum(plan["error"] is not None for plan in plans)
            print(f"dry-run summary: parsed={len(files)} failed={dry_failures} collisions={collision_count}")
            if dry_failures:
                exit_code = 1
        else:
            print(f"summary: converted={counts['converted']} skipped={counts['skipped']} failed={counts['failed']} collisions={collision_count}")
            manifest_path = Path(args.manifest).expanduser().resolve() if args.manifest else output_root / "conversion_manifest.tsv"
            write_manifest(manifest_path, plans)
            if counts["failed"]:
                exit_code = 1
    finally:
        if log_file is not None:
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            log_file.close()

    if exit_code:
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
