# BIDS To OPPNI

Create OPPNI-B or OPPNI-D input files from a BIDS dataset.

## Run

Use `bids_patterns.txt` for general BIDS filename patterns. Use
`compass_bids_patterns.txt` for shared COMPASS-ND BOLD and DWI rules.

```bash
python bids_to_oppni_b_input.py \
  --bids-root BIDS_INPUT \
  --patterns compass_bids_patterns.txt \
  --output input_auto.txt
```

The generator resolves `TR_MSEC` and `TPATTERN` separately for each functional
run from explicit settings, JSON metadata, pattern rules, or defaults. It
creates slice-timing files beside the output when `SliceTiming` is available.

To print the selected files, metadata, resolved values, and their sources:

```bash
python bids_to_oppni_b_input.py \
  --bids-root BIDS_INPUT \
  --patterns compass_bids_patterns.txt \
  --output input_auto.txt \
  --inspect
```

Inputs that cannot produce a valid OPPNI row are skipped with their reasons,
such as missing T1 anatomy, missing functional data, or missing distortion
inputs. Metadata errors are reported immediately and stop output creation.

Run the following for the complete list of CLI options:

```bash
python bids_to_oppni_b_input.py --help
```

## Diffusion

```bash
python bids_to_oppni_d_input.py \
  --bids-root BIDS_INPUT \
  --patterns compass_bids_patterns.txt \
  --output input_dwi.txt \
  --inspect
```

The shared pattern file contains common anatomy and modality-specific sections.
The BOLD generator reads the BOLD sections; the DWI generator reads the DWI
sections. `REV_MODE=NONE` is supported for acquisitions without reverse PE.

The DWI generator is a fast mapper by default. It checks that required files
and metadata are present, but does not open NIfTI files or parse bval/bvec
contents. Use `--validate` when full NIfTI, bval, and bvec consistency checks
are needed.

Use the following for the complete DWI CLI options:

```bash
python bids_to_oppni_d_input.py --help
```

## Pattern files

Pattern files describe dataset-specific filename patterns and metadata
fallbacks. Common sections include:

- `[anat_patterns]` - T1 anatomy used by both generators
- `[func_patterns]` - BOLD files used by OPPNI-B
- `[dwi_patterns]` - diffusion files used by OPPNI-D
- `[func_reverse_pe_patterns]` - reverse-PE BOLD files
- `[dwi_reverse_pe_patterns]` - reverse-PE diffusion files

Scanner-specific metadata can provide values missing from BIDS JSON sidecars.
For example:

```ini
[bold manufacturer:Siemens]
tpattern=parity_alt+z
slice_axis=k
tr_msec=[2130]

[dwi manufacturer:Siemens]
pe_fwd=A>>P
tro_msec=39.2393
rev_mode=NONE
```

These are dataset-specific examples, not universal Siemens settings. See
`compass_bids_patterns.txt` for a complete COMPASS-ND example.

`parity_alt+z` resolves to `alt+z` for an odd number of slices and `alt+z2`
for an even number. When JSON provides `SliceTiming`, the generator uses the
exact generated timing file instead.
