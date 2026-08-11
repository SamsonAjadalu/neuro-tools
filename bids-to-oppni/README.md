# BIDS To OPPNI

Create an OPPNI-B input file from a BIDS dataset.

## Run

Use `bids_patterns.txt` for general BIDS filename patterns. Use
`compass_bids_patterns.txt` for COMPASS-ND scanner rules.

```bash
python bids_to_oppni_input.py \
  --bids-root BIDS_INPUT \
  --patterns compass_bids_patterns.txt \
  --output input_auto.txt
```

The generator resolves `TR_MSEC` and `TPATTERN` separately for each functional
run from explicit settings, JSON metadata, pattern rules, or defaults. It
creates slice-timing files beside the output when `SliceTiming` is available.

To print the selected files, metadata, resolved values, and their sources:

```bash
python bids_to_oppni_input.py \
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
python bids_to_oppni_input.py --help
```
