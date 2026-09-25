# make-nifti-mosaic

Render compatible 3D NIfTI underlay/overlay pairs as PNG mosaics.

## Installation

```bash
bash setup.sh
```

Run the tool from this directory:

```bash
python3 make_nifti_mosaic.py --help
```

## Single-file workflow

```bash
python3 make_nifti_mosaic.py \
  --underlay anat.nii.gz \
  --overlay func.nii.gz \
  --output qc.png
```

Required options: `--underlay`, `--overlay`, and `--output`.

## Folder workflow

```bash
python3 make_nifti_mosaic.py \
  --input /path/to/fmri_proc \
  --patterns qc_patterns.txt \
  --output /path/to/QC
```

Each section specifies its underlay, overlay, and output. Use `{subject}` when
the matching files are subject-specific:

```ini
underlay = {subject}/anat_proc/anat1_2std.nii.gz
overlay = {subject}/anat_proc/subpipe_001/anat_procss.nii.gz
output = QC/{subject}.png
```

A single matching file can be used for every subject. A pattern that matches
multiple subject files includes `{subject}` to pair them explicitly.

The default patterns example is [qc_patterns.txt](qc_patterns.txt).

## Options

```text
--mask PATH
--sagittal-mm MM [MM ...]
--coronal-mm MM [MM ...]
--axial-mm MM [MM ...]
--n-slices INT                    default: 10
--min-active-percent FLOAT        default: 1
--min-active-voxels INT           default: 10
--underlay-cmap {gray,jet,RdBu_r} default: gray
--overlay-cmap {gray,jet,RdBu_r}  default: jet
--overlay-alpha FLOAT             default: 0.5
--vmin FLOAT --vmax FLOAT
--title TEXT
--dpi INT                         default: 300
```

`S`, `C`, and `A` represent sagittal/world x, coronal/world y, and axial/world
z. Manual coordinates use millimetres.

Automatic selection produces axial views. Masked selection uses active mask
voxels. Unmasked selection uses axial sums of absolute overlay values.

Binary overlays use `[0, 1]` automatically. Continuous overlays use
percentile-based display
limits when `--vmin` and `--vmax` are omitted; provide them together to set
explicit limits.

## Patterns file

Each section defines one output configuration. Underlay, overlay, and mask
entries resolve to `.nii`
or `.nii.gz` files. The output entry resolves to a `.png` file.

```ini
[func_anat]
underlay = {subject}/anat_proc/subpipe_001/anat_warped.nii.gz
overlay = {subject}/func_proc_p1/subpipe_001/postwarp/func1_warped_tav.nii.gz
output = {subject}/func_anat_alignment.png
overlay_alpha = 0.33
sagittal_mm = -50 -8 30
coronal_mm = -65 -20 54
axial_mm = -6 13 58
```

Supported section settings:

```text
mask
sagittal_mm
coronal_mm
axial_mm
n_slices
min_active_percent
min_active_voxels
underlay_cmap
overlay_cmap
overlay_alpha
vmin
vmax
title
dpi
```

## Precedence

```text
CLI > pattern section > built-in default
```

Coordinate lists follow the same order. A higher-priority list replaces the
lower-priority list.
