# Source To BIDS

Convert configured MINC files, or organize existing NIfTI files, into a BIDS dataset.

## Run

Use `compass_to_bids_config.txt` for COMPASS-ND. Use `source_to_bids_config.txt` as a template for another dataset.

`source_extension` accepts one extension or a comma-separated list, for example:

```ini
source_extension=.mnc,.nii,.nii.gz
```

The converter selects MINC or NIfTI separately for each file from its extension.

```bash
python source_to_bids.py INPUT_FOLDER \
  --config compass_to_bids_config.txt \
  --output BIDS_OUTPUT
```

The input may be one file or a folder. A normal run performs preflight checks before conversion and writes `conversion.log` and `conversion_manifest.tsv` in the output folder.

To check filenames, headers, dimensions, DWI metadata, and destination collisions without writing outputs:

```bash
python source_to_bids.py INPUT_FOLDER \
  --config compass_to_bids_config.txt \
  --output BIDS_OUTPUT \
  --dry-run
```

Use `--overwrite` only when replacing an existing complete output is intentional. A collision between two source files in the same batch always fails.
