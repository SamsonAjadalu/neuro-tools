# seed-frequency

Create positive and negative voxelwise frequency maps from P4 t-statistic
NIfTI files. The input maps can be anywhere and may use any directory or seed
naming convention.

From the repository root:

```bash
python3 seed-frequency/seed_frequency.py \
  --maps '/path/to/results/**/tstat_unthresh.nii' \
  --contrast 0 \
  --threshold 2 \
  --output ./frequency
```

The tool writes `positive_frequency.nii.gz` and `negative_frequency.nii.gz`.
Run `bash setup.sh` once first to install the shared `.venv` dependencies.
