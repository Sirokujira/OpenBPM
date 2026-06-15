# Numerical regression tests

`check_regression.py` guards against silent numerical breakage in the OpenBPM
solver. For each sample input it reads the HDF5 result, computes stable scalar
metrics from the final field, and compares them against the committed golden
values in `reference.json` within a tolerance. For a curated subset it also
asserts a physical relation against the analytic theory documented in
`ReadMe.md`.

## Metrics (all samples)

Computed from the final complex field `Efinal` (intensity `I = |E|^2`) using
the physical grid coordinates `Xc`/`Yc` stored in the result:

| metric | meaning |
|---|---|
| `power` | total output power `sum(I)` |
| `peak`  | peak intensity `max(I)` |
| `cx`,`cy` | intensity centroid (m) — detects beam steering (tilt/bend/offset) |
| `wx`,`wy` | 1/e² beam width `2·sqrt(second moment)` (m) — detects focusing/diffraction |

Tolerance (see `_meta.tol` in `reference.json`): relative `1e-2`, with an
absolute floor of `1e-8 m` for the metre-valued metrics so that near-zero
centroids do not trip the relative check.

## Physics checks (curated subset)

These assert the solver reproduces the analytic result, not merely that it is
unchanged:

- **freespace** — Gaussian diffraction `w(z) = w0·sqrt(1+(z/zR)²)`; measured
  output width agrees to ~0.5% (tol 3%).
- **tilt_wideangle** — wide-angle BPM; centroid shift converges to `z·tan(20°)`
  (21.8 µm), distinct from the paraxial `z·sin(20°)` (tol 5%).

## Reproducibility

Results are bit-identical only with a single OpenMP thread. Always run the
solver with `OMP_NUM_THREADS=1` when generating or checking, as the CI does.

## Usage

```sh
# check a result against the reference
OMP_NUM_THREADS=1 ./bin/obpm data/freespace.ofd
python3 tests/regression/check_regression.py --name freespace --h5 time_series_data.h5

# regenerate / add golden values (run with a known-good build, then commit)
python3 tests/regression/check_regression.py --name freespace --h5 time_series_data.h5 --update
```

Requires `h5py` and `numpy`. Note: Debian's `python3-h5py` package can fail to
import (serial/mpi auto-selection bug); install via `pip` in a venv instead, as
the CI workflow does.
