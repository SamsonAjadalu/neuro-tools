#!/usr/bin/env python3
"""Create positive and negative frequency maps from NIfTI t-statistic maps."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VENV_ROOT = ROOT / ".venv"
VENV_PYTHON = VENV_ROOT / "bin" / "python"

if not VENV_PYTHON.is_file():
    raise SystemExit("Environment not installed.\nRun: bash setup.sh")

if Path(sys.prefix).resolve() != VENV_ROOT.resolve():
    try:
        os.execv(str(VENV_PYTHON), [str(VENV_PYTHON), *sys.argv])
    except OSError as exc:
        raise SystemExit(
            "Required dependencies are missing or the environment is incomplete.\n"
            "Run: bash setup.sh"
        ) from exc


try:
    import nibabel as nib
    import numpy as np
except (ImportError, ModuleNotFoundError) as exc:
    raise SystemExit(
        "Required dependencies are missing or the environment is incomplete.\n"
        "Run: bash setup.sh"
    ) from exc


AFFINE_ATOL = 1e-5
AFFINE_RTOL = 1e-5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create positive and negative voxelwise seed-frequency maps."
    )
    parser.add_argument(
        "--maps", required=True, help="Recursive glob pattern for input NIfTI maps"
    )
    parser.add_argument(
        "--contrast",
        type=int,
        default=0,
        help="Zero-based contrast volume for 4D maps (default: 0)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=2.0,
        help="Positive t-statistic threshold (default: 2)",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Output directory"
    )
    args = parser.parse_args()
    if args.contrast < 0:
        parser.error("--contrast must be zero or greater")
    if not np.isfinite(args.threshold):
        parser.error("--threshold must be finite")
    return args


def matched_paths(pattern: str) -> list[Path]:
    paths = sorted(Path().glob(pattern) if not os.path.isabs(pattern) else Path("/").glob(pattern.lstrip("/")))
    resolved = [path.resolve() for path in paths]
    if len(set(resolved)) != len(resolved):
        raise ValueError("input pattern matched duplicate files")
    if not paths:
        raise ValueError(f"input pattern matched no files: {pattern}")
    return paths


def selected_volume(path: Path, contrast: int) -> tuple[np.ndarray, nib.Nifti1Image]:
    try:
        image = nib.load(str(path))
    except Exception as exc:
        raise ValueError(f"could not read NIfTI file: {path}: {exc}") from exc

    if len(image.shape) == 3:
        if contrast != 0:
            raise ValueError(
                f"contrast {contrast} is invalid for 3D image: {path}"
            )
        data = np.asanyarray(image.dataobj)
    elif len(image.shape) == 4:
        if contrast >= image.shape[3]:
            raise ValueError(
                f"contrast {contrast} is out of range for {path} with "
                f"{image.shape[3]} volumes"
            )
        data = np.asanyarray(image.dataobj[..., contrast])
    else:
        raise ValueError(
            f"input must be 3D or 4D, got shape {image.shape}: {path}"
        )

    return np.asarray(data, dtype=np.float32), image


def create_frequency_maps(
    paths: list[Path], contrast: int, threshold: float
) -> tuple[np.ndarray, np.ndarray, nib.Nifti1Image]:
    positive = negative = None
    reference = None

    for path in paths:
        data, image = selected_volume(path, contrast)
        if reference is None:
            reference = image
            positive = np.zeros(data.shape, dtype=np.int32)
            negative = np.zeros(data.shape, dtype=np.int32)
        else:
            if data.shape != positive.shape:
                raise ValueError(
                    f"incompatible dimensions: {path} has {data.shape}, "
                    f"expected {positive.shape}"
                )
            if not np.allclose(
                image.affine,
                reference.affine,
                rtol=AFFINE_RTOL,
                atol=AFFINE_ATOL,
            ):
                raise ValueError(f"incompatible affine: {path}")

        positive += data > threshold
        negative += data < -threshold

    denominator = float(len(paths))
    return positive / denominator, negative / denominator, reference


def save_map(data: np.ndarray, reference: nib.Nifti1Image, path: Path) -> None:
    header = reference.header.copy()
    header.set_data_dtype(np.float32)
    image = nib.Nifti1Image(data.astype(np.float32), reference.affine, header=header)
    image.set_qform(reference.get_qform(), code=reference.get_qform(coded=True)[1])
    image.set_sform(reference.get_sform(), code=reference.get_sform(coded=True)[1])
    nib.save(image, str(path))


def main() -> None:
    args = parse_args()
    paths = matched_paths(args.maps)
    positive, negative, reference = create_frequency_maps(
        paths, args.contrast, args.threshold
    )

    args.output.mkdir(parents=True, exist_ok=True)
    positive_path = args.output / "positive_frequency.nii.gz"
    negative_path = args.output / "negative_frequency.nii.gz"
    save_map(positive, reference, positive_path)
    save_map(negative, reference, negative_path)

    print("Matched input files:")
    for path in paths:
        print(f"  {path}")
    print(f"Number of maps: {len(paths)}")
    print(f"Contrast: {args.contrast}")
    print(f"Threshold: {args.threshold}")
    print(f"Positive output: {positive_path}")
    print(f"Negative output: {negative_path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
