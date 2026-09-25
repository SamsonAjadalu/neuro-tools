#!/usr/bin/env python3
"""Render a generic 3D NIfTI underlay/overlay mosaic."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
VENV_ROOT = ROOT / ".venv"
VENV_PYTHON = VENV_ROOT / "bin" / "python"

if not VENV_PYTHON.is_file():
    raise SystemExit(
        "Environment not installed.\n"
        "Run: bash setup.sh"
    )

if Path(sys.prefix).resolve() != VENV_ROOT.resolve():
    try:
        os.execv(str(VENV_PYTHON), [str(VENV_PYTHON), *sys.argv])
    except OSError as exc:
        raise SystemExit(
            "Required dependencies are missing or the environment is incomplete.\n"
            "Run: bash setup.sh"
        ) from exc


try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import nibabel as nib
    import numpy as np
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
except (ImportError, ModuleNotFoundError) as exc:
    raise SystemExit(
        "Required dependencies are missing or the environment is incomplete.\n"
        "Run: bash setup.sh"
    ) from exc


AFFINE_ATOL = 1e-5
AFFINE_RTOL = 1e-5
SUPPORTED_CMAPS = ("gray", "jet", "RdBu_r")
ORIENTATIONS = {"S": ("x", 0), "C": ("y", 1), "A": ("z", 2)}
ORIENTATION_NAMES = {"S": "sagittal", "C": "coronal", "A": "axial"}
DEFAULTS = {
    "n_slices": 10,
    "min_active_percent": 1.0,
    "min_active_voxels": 10,
    "underlay_cmap": "gray",
    "overlay_cmap": "jet",
    "overlay_alpha": 0.5,
    "dpi": 300,
}
PATTERN_PATH_KEYS = {"underlay", "overlay", "mask", "output"}
PATTERN_VALUE_KEYS = {
    "sagittal_mm",
    "coronal_mm",
    "axial_mm",
    "n_slices",
    "min_active_percent",
    "min_active_voxels",
    "underlay_cmap",
    "overlay_cmap",
    "overlay_alpha",
    "vmin",
    "vmax",
    "title",
    "dpi",
}
PATTERN_KEYS = PATTERN_PATH_KEYS | PATTERN_VALUE_KEYS


@dataclass(frozen=True)
class View:
    orientation: str
    axis: str
    mm_position: float
    voxel_index: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a generic 3D NIfTI underlay/overlay mosaic."
    )
    parser.add_argument("--underlay", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--patterns", type=Path)
    parser.add_argument("--mask", type=Path, help="Optional binary overlay-visibility mask.")
    parser.add_argument("--sagittal-mm", type=float, nargs="+")
    parser.add_argument("--coronal-mm", type=float, nargs="+")
    parser.add_argument("--axial-mm", type=float, nargs="+")
    parser.add_argument("--n-slices", type=int)
    parser.add_argument("--min-active-percent", type=float)
    parser.add_argument("--min-active-voxels", type=int)
    parser.add_argument("--underlay-cmap", choices=SUPPORTED_CMAPS)
    parser.add_argument("--overlay-cmap", choices=SUPPORTED_CMAPS)
    parser.add_argument("--overlay-alpha", type=float)
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--title")
    parser.add_argument("--dpi", type=int)
    return parser.parse_args()


def is_nifti_path(path):
    value = str(path).lower()
    return value.endswith(".nii") or value.endswith(".nii.gz")


def parse_pattern_value(key, value, section):
    value = value.strip()
    if key in {"sagittal_mm", "coronal_mm", "axial_mm"}:
        try:
            return [float(item) for item in value.split()]
        except ValueError as exc:
            raise ValueError(f"section [{section}] has invalid {key} values") from exc
    if key in {"n_slices", "min_active_voxels", "dpi"}:
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"section [{section}] has invalid integer {key}") from exc
    if key in {"min_active_percent", "overlay_alpha", "vmin", "vmax"}:
        try:
            return float(value)
        except ValueError as exc:
            raise ValueError(f"section [{section}] has invalid numeric {key}") from exc
    return value


def parse_patterns(path):
    if not path.is_file():
        raise FileNotFoundError(f"patterns file does not exist: {path}")
    sections = {}
    current = None
    try:
        lines = path.read_text().splitlines()
    except OSError as exc:
        raise ValueError(f"could not read patterns file: {path}") from exc
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            if not current:
                raise ValueError(f"patterns file line {line_number}: empty section name")
            if current in sections:
                raise ValueError(f"patterns file has duplicate section [{current}]")
            sections[current] = {}
            continue
        if current is None or "=" not in line:
            raise ValueError(f"patterns file line {line_number}: expected [section] or key=value")
        key, value = (part.strip() for part in line.split("=", 1))
        if key not in PATTERN_KEYS:
            raise ValueError(f"section [{current}] has invalid key: {key}")
        if not value:
            raise ValueError(f"section [{current}] has an empty value for {key}")
        sections[current][key] = parse_pattern_value(key, value, current)
    if not sections:
        raise ValueError(f"patterns file contains no sections: {path}")
    for name, section in sections.items():
        for required in ("underlay", "overlay", "output"):
            if required not in section:
                raise ValueError(f"section [{name}] is missing required key: {required}")
        for key in ("underlay", "overlay", "mask"):
            if key in section and not is_nifti_path(section[key]):
                raise ValueError(
                    f"{key} pattern must resolve to a .nii or .nii.gz file"
                )
        if not str(section["output"]).lower().endswith(".png"):
            raise ValueError(f"section [{name}] output must be a PNG path")
    return sections


def discover_section_cases(input_root, section_name, section):
    directories = [
        path for path in input_root.iterdir()
        if path.is_dir() and not path.name.startswith((".", "_"))
    ]
    subject_directories = [path for path in directories if path.name.startswith("sub-")]
    subjects = sorted(path.name for path in (subject_directories or directories))
    if not any("{subject}" in str(section.get(key, "")) for key in ("underlay", "overlay", "mask", "output")):
        subjects = ["."]

    def resolve_pattern(pattern, subject):
        has_subject = "{subject}" in pattern
        if has_subject:
            try:
                pattern = pattern.format(subject=subject)
            except KeyError as exc:
                raise ValueError(f"section [{section_name}] has an unsupported pattern field: {exc.args[0]}") from exc
        direct = input_root / pattern
        if direct.is_file():
            return [direct]
        if has_subject:
            return []
        return sorted(path for path in input_root.rglob(pattern) if path.is_file())

    cases = []
    skipped = 0
    for subject in subjects:
        matches = {
            key: resolve_pattern(section[key], subject)
            for key in ("underlay", "overlay")
        }
        if "mask" in section:
            matches["mask"] = resolve_pattern(section["mask"], subject)
        for key in ("underlay", "overlay", "mask"):
            if key in section and "{subject}" not in section[key] and len(matches[key]) > 1:
                raise ValueError(
                    f"{key} pattern must include {{subject}} when it matches multiple subject files"
                )
        if len(matches["underlay"]) != 1:
            print(f"[{section_name}] skipped {subject}: expected one underlay match")
            skipped += 1
            continue
        if len(matches["overlay"]) != 1:
            print(f"[{section_name}] skipped {subject}: expected one overlay match")
            skipped += 1
            continue
        if "mask" in section and len(matches["mask"]) != 1:
            print(f"[{section_name}] skipped {subject}: expected one mask match")
            skipped += 1
            continue
        cases.append(
            {
                "case": subject,
                "underlay": matches["underlay"][0],
                "overlay": matches["overlay"][0],
                "mask": matches.get("mask", [None])[0],
            }
        )
    return cases, skipped


def resolve_case_output(output_root, case, output_pattern):
    case_root = output_root if case == "." else output_root / case
    if "{subject}" in output_pattern:
        return output_root / output_pattern.format(subject=case)
    return case_root / output_pattern


def apply_pattern_settings(args, section):
    """Apply section values only where the CLI did not provide a value."""
    option_names = {
        "sagittal_mm": "sagittal_mm",
        "coronal_mm": "coronal_mm",
        "axial_mm": "axial_mm",
        "n_slices": "n_slices",
        "min_active_percent": "min_active_percent",
        "min_active_voxels": "min_active_voxels",
        "underlay_cmap": "underlay_cmap",
        "overlay_cmap": "overlay_cmap",
        "overlay_alpha": "overlay_alpha",
        "vmin": "vmin",
        "vmax": "vmax",
        "title": "title",
        "dpi": "dpi",
    }
    for pattern_key, arg_name in option_names.items():
        if getattr(args, arg_name) is None and pattern_key in section:
            setattr(args, arg_name, section[pattern_key])


def build_case_args(base_args, section, match, output_root):
    args = argparse.Namespace(**vars(base_args))
    args.underlay = match["underlay"]
    args.overlay = match["overlay"]
    args.mask = base_args.mask if base_args.mask is not None else match["mask"]
    args.output = resolve_case_output(output_root, match["case"], section["output"])
    apply_pattern_settings(args, section)
    return args


def load_and_validate(args: argparse.Namespace):
    for label, path in (("underlay", args.underlay), ("overlay", args.overlay)):
        if not path.is_file():
            raise FileNotFoundError(f"{label} NIfTI does not exist: {path}")
    if args.mask is not None and not args.mask.is_file():
        raise FileNotFoundError(f"mask NIfTI does not exist: {args.mask}")

    underlay_img = nib.load(str(args.underlay))
    overlay_img = nib.load(str(args.overlay))
    mask_img = nib.load(str(args.mask)) if args.mask is not None else None
    images = [("underlay", underlay_img), ("overlay", overlay_img)]
    if mask_img is not None:
        images.append(("mask", mask_img))

    for label, image in images:
        if len(image.shape) != 3:
            raise ValueError(f"{label} must be 3D; found shape {image.shape}")
        if not np.isfinite(image.affine).all():
            raise ValueError(f"{label} has a non-finite affine")
    if underlay_img.shape != overlay_img.shape:
        raise ValueError(
            f"underlay and overlay shapes differ: {underlay_img.shape} vs {overlay_img.shape}"
        )
    if mask_img is not None and mask_img.shape != overlay_img.shape:
        raise ValueError(
            f"mask shape {mask_img.shape} does not match overlay shape {overlay_img.shape}"
        )
    for label, image in images[1:]:
        if not np.allclose(
            underlay_img.affine, image.affine, atol=AFFINE_ATOL, rtol=AFFINE_RTOL
        ):
            raise ValueError(f"underlay and {label} affines are incompatible")

    underlay = np.asanyarray(underlay_img.dataobj, dtype=np.float64)
    overlay = np.asanyarray(overlay_img.dataobj, dtype=np.float64)
    mask = np.asanyarray(mask_img.dataobj, dtype=np.float64) if mask_img is not None else None
    for label, data in (("underlay", underlay), ("overlay", overlay), ("mask", mask)):
        if data is not None and not np.isfinite(data).all():
            raise ValueError(f"{label} contains NaN or infinite values")
    if mask is not None and not np.all(np.isin(np.unique(mask), [0.0, 1.0])):
        raise ValueError("mask must be strictly binary and contain only 0 and 1")

    underlay_img = nib.as_closest_canonical(underlay_img)
    overlay_img = nib.as_closest_canonical(overlay_img)
    if mask_img is not None:
        mask_img = nib.as_closest_canonical(mask_img)
    underlay = np.asanyarray(underlay_img.dataobj, dtype=np.float64)
    overlay = np.asanyarray(overlay_img.dataobj, dtype=np.float64)
    mask = np.asanyarray(mask_img.dataobj, dtype=np.float64) if mask_img is not None else None
    return underlay_img, overlay_img, mask_img, underlay, overlay, mask


def is_binary(data: np.ndarray) -> bool:
    return bool(np.all(np.isin(np.unique(data), [0.0, 1.0])))


def visible_overlay(overlay: np.ndarray, mask: np.ndarray | None) -> np.ndarray:
    return overlay if mask is None else overlay[mask == 1]


def compute_display_limits(overlay, mask, requested_vmin, requested_vmax):
    if (requested_vmin is None) != (requested_vmax is None):
        raise ValueError("--vmin and --vmax must be supplied together")
    if requested_vmin is not None:
        if not np.isfinite([requested_vmin, requested_vmax]).all() or requested_vmin >= requested_vmax:
            raise ValueError("display range requires finite --vmin < --vmax")
        return float(requested_vmin), float(requested_vmax)
    if is_binary(overlay):
        return 0.0, 1.0
    values = visible_overlay(overlay, mask)
    values = values[np.isfinite(values) & (values != 0)]
    if values.size == 0:
        raise ValueError("overlay has no visible nonzero values for display normalization")
    low, high = np.percentile(values, [2.5, 97.5])
    if low == high:
        low, high = float(values.min()), float(values.max())
    if low == high:
        delta = max(abs(float(low)) * 0.01, 1.0)
        low, high = float(low) - delta, float(high) + delta
    return float(low), float(high)


def compute_underlay_limits(underlay):
    values = underlay[np.isfinite(underlay) & (underlay != 0)]
    if values.size == 0:
        values = underlay[np.isfinite(underlay)]
    if values.size == 0:
        raise ValueError("underlay contains no finite values")
    low, high = np.percentile(values, [2.5, 97.5])
    if low == high:
        low, high = float(values.min()), float(values.max())
    if low == high:
        delta = max(abs(float(low)) * 0.01, 1.0)
        low, high = float(low) - delta, float(high) + delta
    return float(low), float(high)


def world_coordinate_bounds(affine, shape, world_axis):
    corners = np.array(np.meshgrid(*[(0, size - 1) for size in shape], indexing="ij")).reshape(3, -1).T
    world = nib.affines.apply_affine(affine, corners)
    values = world[:, world_axis]
    return float(values.min()), float(values.max())


def world_mm_to_voxel_index(affine, shape, orientation, mm_position):
    world_axis = ORIENTATIONS[orientation][1]
    low, high = world_coordinate_bounds(affine, shape, world_axis)
    if mm_position < low - AFFINE_ATOL or mm_position > high + AFFINE_ATOL:
        axis = ORIENTATIONS[orientation][0]
        raise ValueError(
            f"{orientation} coordinate {mm_position:g} mm is outside the {axis}-range [{low:g}, {high:g}] mm"
        )
    center_voxel = (np.asarray(shape, dtype=float) - 1.0) / 2.0
    target_world = nib.affines.apply_affine(affine, center_voxel)
    target_world[world_axis] = mm_position
    voxel = nib.affines.apply_affine(np.linalg.inv(affine), target_world)
    index = int(np.rint(voxel[world_axis]))
    if index < 0 or index >= shape[world_axis]:
        raise ValueError(f"{orientation} coordinate {mm_position:g} mm maps outside the image")
    return index


def view_mm_position(affine, orientation, voxel_index, shape):
    voxel = (np.asarray(shape, dtype=float) - 1.0) / 2.0
    axis = ORIENTATIONS[orientation][1]
    voxel[axis] = voxel_index
    return float(nib.affines.apply_affine(affine, voxel)[axis])


def manual_views(args, affine, shape):
    views = []
    for orientation, positions in (("S", args.sagittal_mm), ("C", args.coronal_mm), ("A", args.axial_mm)):
        if positions is None:
            continue
        for position in positions:
            index = world_mm_to_voxel_index(affine, shape, orientation, float(position))
            actual_mm_position = view_mm_position(affine, orientation, index, shape)
            views.append(
                View(
                    orientation,
                    ORIENTATIONS[orientation][0],
                    actual_mm_position,
                    index,
                )
            )
    return views


def automatic_axial_views(overlay, mask, affine, n_slices, min_active_percent, min_active_voxels):
    if mask is not None:
        activity = np.count_nonzero(mask > 0, axis=(0, 1))
        maximum = int(activity.max())
        if maximum == 0:
            raise ValueError("mask contains no active voxels")
        percentage_threshold = (min_active_percent / 100.0) * maximum
        threshold = max(percentage_threshold, min_active_voxels)
        eligible = np.flatnonzero(activity >= threshold)
        if eligible.size == 0:
            raise ValueError(
                "no axial slice qualifies: "
                f"maximum activity={maximum}, percentage threshold={percentage_threshold:g}, "
                f"min-active-voxels={min_active_voxels}, final threshold={threshold:g}"
            )
    else:
        activity = np.sum(np.abs(overlay), axis=(0, 1))
        maximum = float(activity.max())
        if maximum == 0:
            raise ValueError("overlay is all zero; automatic axial selection is impossible")
        threshold = (min_active_percent / 100.0) * maximum
        eligible = np.flatnonzero(activity >= threshold)
        if eligible.size == 0:
            raise ValueError(
                f"no axial slice qualifies: maximum activity={maximum:g}, percentage threshold={threshold:g}"
            )
    first, last = int(eligible[0]), int(eligible[-1])
    selected = list(dict.fromkeys(np.rint(np.linspace(first, last, n_slices)).astype(int).tolist()))
    return [
        View("A", "z", view_mm_position(affine, "A", index, overlay.shape[0:3]), index)
        for index in selected
    ]


def extract_plane(data, orientation, index):
    axis = ORIENTATIONS[orientation][1]
    plane = data[index, :, :] if axis == 0 else data[:, index, :] if axis == 1 else data[:, :, index]
    return np.rot90(plane)


def balanced_grid(n_views, max_columns=5):
    rows = math.ceil(n_views / max_columns)
    columns = math.ceil(n_views / rows)
    return rows, columns


def render(args, underlay, overlay, mask, views, vmin, vmax):
    max_columns = 5
    grouped_views = {
        orientation: [view for view in views if view.orientation == orientation]
        for orientation in ORIENTATIONS
    }
    present = [(orientation, row_views) for orientation, row_views in grouped_views.items() if row_views]
    equal_counts = len(present) > 1 and len({len(row_views) for _, row_views in present}) == 1
    underlay_low, underlay_high = compute_underlay_limits(underlay)
    overlay_norm = Normalize(vmin=vmin, vmax=vmax)
    overlay_cmap = plt.get_cmap(args.overlay_cmap).copy()
    overlay_cmap.set_bad(alpha=0)
    underlay_cmap = plt.get_cmap(args.underlay_cmap)
    if mask is None:
        visible = np.where(overlay == 0, np.nan, overlay)
    else:
        visible = np.where(mask == 1, overlay, np.nan)
    image_axes = []

    def draw_panel(panel, view, title):
        image_axes.append(panel)
        panel.imshow(
            extract_plane(underlay, view.orientation, view.voxel_index),
            cmap=underlay_cmap,
            vmin=underlay_low,
            vmax=underlay_high,
            interpolation="nearest",
        )
        panel.imshow(
            extract_plane(visible, view.orientation, view.voxel_index),
            cmap=overlay_cmap,
            norm=overlay_norm,
            alpha=args.overlay_alpha,
            interpolation="nearest",
        )
        panel.set_title(title)
        panel.axis("off")

    figure_title = args.title
    if len(present) == 1:
        orientation_name = ORIENTATION_NAMES[present[0][0]].capitalize()
        figure_title = f"{args.title} | {orientation_name}" if args.title else orientation_name

    if len(present) > 1 and not equal_counts:
        flat_views = [view for _, row_views in present for view in row_views]
        n_rows, columns = balanced_grid(len(flat_views), max_columns)
        layout_rows = [flat_views[start : start + columns] for start in range(0, len(flat_views), columns)]
        fig = plt.figure(figsize=(3.2 * columns, 3.2 * len(layout_rows)), constrained_layout=True)
        grid = fig.add_gridspec(len(layout_rows), columns)
        for row, row_views in enumerate(layout_rows):
            for column, view in enumerate(row_views):
                panel = fig.add_subplot(grid[row, column])
                axis_name = ORIENTATIONS[view.orientation][0]
                draw_panel(panel, view, f"{view.orientation} {axis_name}={view.mm_position:+.1f} mm")
            for column in range(len(row_views), columns):
                empty_axis = fig.add_subplot(grid[row, column])
                empty_axis.axis("off")
    elif len(present) == 1:
        orientation, orientation_views = present[0]
        n_rows, columns = balanced_grid(len(orientation_views), max_columns)
        layout_rows = [
            orientation_views[start : start + columns]
            for start in range(0, len(orientation_views), columns)
        ]
        fig = plt.figure(figsize=(3.2 * columns, 3.2 * n_rows), constrained_layout=True)
        grid = fig.add_gridspec(n_rows, columns)
        for row, row_views in enumerate(layout_rows):
            for column, view in enumerate(row_views):
                panel = fig.add_subplot(grid[row, column])
                draw_panel(panel, view, f"{view.mm_position:+.1f} mm")
            for column in range(len(row_views), columns):
                empty_axis = fig.add_subplot(grid[row, column])
                empty_axis.axis("off")
    else:
        layout_rows = []
        for orientation, orientation_views in present:
            for start in range(0, len(orientation_views), max_columns):
                layout_rows.append((orientation, orientation_views[start : start + max_columns], start == 0))
        columns = max(len(row_views) for _, row_views, _ in layout_rows)
        label_width = 0.65
        fig = plt.figure(figsize=(3.2 * columns + label_width, 3.2 * len(layout_rows)), constrained_layout=True)
        grid = fig.add_gridspec(
            len(layout_rows),
            columns + 1,
            width_ratios=[label_width] + [1.0] * columns,
        )
        for row, (orientation, row_views, show_label) in enumerate(layout_rows):
            label_axis = fig.add_subplot(grid[row, 0])
            label_axis.axis("off")
            if show_label:
                label = (
                    ORIENTATION_NAMES[orientation]
                    if len(present) == 1
                    else f"{orientation} ({ORIENTATIONS[orientation][0]}, mm)"
                )
                label_axis.text(1.0, 0.5, label, ha="right", va="center", transform=label_axis.transAxes)
            for column, view in enumerate(row_views):
                panel = fig.add_subplot(grid[row, column + 1])
                title = f"{view.mm_position:+.1f} mm"
                draw_panel(panel, view, title)
            for column in range(len(row_views), columns):
                empty_axis = fig.add_subplot(grid[row, column + 1])
                empty_axis.axis("off")

    scalar_map = ScalarMappable(norm=overlay_norm, cmap=overlay_cmap)
    scalar_map.set_array([])
    fig.colorbar(scalar_map, ax=image_axes, shrink=0.82, pad=0.02)
    if figure_title:
        fig.suptitle(figure_title, fontsize=14)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def validate_args(args):
    if args.n_slices <= 0:
        raise ValueError("--n-slices must be positive")
    if args.min_active_percent < 0:
        raise ValueError("--min-active-percent must be >= 0")
    if args.min_active_voxels < 0:
        raise ValueError("--min-active-voxels must be >= 0")
    if not 0 <= args.overlay_alpha <= 1:
        raise ValueError("--overlay-alpha must be between 0 and 1")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")


def process_single(args, print_result=True):
    _, overlay_img, _, underlay, overlay, mask = load_and_validate(args)
    shape = tuple(overlay.shape)
    affine = overlay_img.affine
    for name, default in DEFAULTS.items():
        if getattr(args, name) is None:
            setattr(args, name, default)
    validate_args(args)
    manual_supplied = any(x is not None for x in (args.sagittal_mm, args.coronal_mm, args.axial_mm))
    if manual_supplied:
        views, selection_method = manual_views(args, affine, shape), "manual"
    else:
        views = automatic_axial_views(overlay, mask, affine, args.n_slices,
                                      args.min_active_percent, args.min_active_voxels)
        selection_method = "automatic"
    if not views:
        raise ValueError("no views were selected")
    vmin, vmax = compute_display_limits(overlay, mask, args.vmin, args.vmax)
    render(args, underlay, overlay, mask, views, vmin, vmax)
    if print_result:
        print(json.dumps({"output": str(args.output), "views": [view.__dict__ for view in views],
                          "vmin": vmin, "vmax": vmax}, indent=2))


def validate_cli_paths(args):
    if args.input is not None:
        if not args.input.is_dir():
            raise ValueError(f"--input must be a directory: {args.input}")
        if args.patterns is None:
            raise ValueError("folder input requires --patterns FILE")
        if args.underlay is not None or args.overlay is not None:
            raise ValueError("folder input cannot be combined with --underlay or --overlay")
        return "folder"
    if args.patterns is not None:
        raise ValueError("--patterns requires --input DIRECTORY")
    if args.underlay is None or args.overlay is None:
        raise ValueError("single-file usage requires --underlay FILE and --overlay FILE")
    for label, path in (("underlay", args.underlay), ("overlay", args.overlay), ("mask", args.mask)):
        if path is not None and not is_nifti_path(path):
            raise ValueError(f"{label} must end in .nii or .nii.gz")
    return "single"


def process_folder(args):
    sections = parse_patterns(args.patterns)
    processed = 0
    skipped = 0
    errors = 0
    for section_name, section in sections.items():
        cases, discovery_skipped = discover_section_cases(args.input, section_name, section)
        skipped += discovery_skipped
        if not cases:
            print(f"[{section_name}] no matches found; skipping")
            if discovery_skipped == 0:
                skipped += 1
            continue
        for match in cases:
            case_args = build_case_args(args, section, match, args.output)
            try:
                process_single(case_args, print_result=False)
                print(f"[{section_name}] processed {match['case']}: {case_args.output}")
                processed += 1
            except (FileNotFoundError, ValueError, OSError) as exc:
                print(f"[{section_name}] ERROR {match['case']}: {exc}")
                errors += 1
    print(f"processed={processed} skipped={skipped} errors={errors}")
    return 0 if errors == 0 else 1


def main() -> int:
    args = parse_args()
    workflow = validate_cli_paths(args)
    if workflow == "folder":
        return process_folder(args)
    for name, default in DEFAULTS.items():
        if getattr(args, name) is None:
            setattr(args, name, default)
    validate_args(args)
    process_single(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
