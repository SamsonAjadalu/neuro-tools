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

Use `--inspect` to print the selected files, metadata, resolved values, and
their sources while the dataset is processed.
