#!/usr/bin/env python3
"""OpenBPM numerical regression check.

Reads an OpenBPM HDF5 result (``time_series_data.h5``), computes a set of
stable scalar metrics from the final field, and compares them against a
committed golden reference (``reference.json``) within a tolerance. For a
curated subset of samples it additionally asserts a physical relation
against the analytic theory documented in ReadMe.md.

Usage
-----
Check a result against the reference::

    check_regression.py --name freespace --h5 time_series_data.h5

Regenerate (or add) the golden values for a sample (run locally with a
known-good build, then commit reference.json)::

    check_regression.py --name freespace --h5 time_series_data.h5 --update

Results are reproducible only with a single OpenMP thread, so always run the
solver with ``OMP_NUM_THREADS=1`` when generating or checking.
"""
import argparse
import json
import math
import os
import sys

import h5py
import numpy as np

REF_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference.json")

# Metrics compared against the golden reference for every sample.
# Position-like metrics (centroids, in metres) are compared with an absolute
# floor so that values near zero do not trip the relative tolerance.
METRIC_FIELDS = ("power", "peak", "cx", "cy", "wx", "wy")
DEFAULT_RTOL = 0.01           # 1% relative tolerance
DEFAULT_ATOL_M = 1.0e-8       # 0.01 um absolute floor for metres-valued metrics
METRES_METRICS = ("cx", "cy", "wx", "wy")


def load(h5path):
    """Return (intensity[Ny,Nx], metadata dict, Ixz or None) from an h5 file."""
    with h5py.File(h5path, "r") as f:
        er = f["field/Efinal_r"][()].astype(np.float64)
        ei = f["field/Efinal_i"][()].astype(np.float64)
        intensity = er * er + ei * ei
        md = {k: np.array(f["metadata"][k]).ravel() for k in f["metadata"]}
        ixz = f["field/Ixz"][()].astype(np.float64) if "field/Ixz" in f else None
    return intensity, md, ixz


def metrics(intensity, md):
    """Compute stable scalar metrics from the final field.

    Efinal is stored as [iy][ix]; Xc (len Nx) is the x axis, Yc (len Ny) the y
    axis. Widths are 1/e^2 radii (= 2*sqrt(second moment)) for a Gaussian.
    """
    xc = md["Xc"]
    yc = md["Yc"]
    total = intensity.sum()
    cx = float((intensity * xc[None, :]).sum() / total)
    cy = float((intensity * yc[:, None]).sum() / total)
    wx = 2.0 * math.sqrt(float((intensity * (xc[None, :] - cx) ** 2).sum() / total))
    wy = 2.0 * math.sqrt(float((intensity * (yc[:, None] - cy) ** 2).sum() / total))
    return {
        "power": float(total),
        "peak": float(intensity.max()),
        "cx": cx,
        "cy": cy,
        "wx": wx,
        "wy": wy,
    }


# --- Physical theory checks (curated subset) ---------------------------------
# Each returns a list of (label, measured, theory, rtol) tuples. These assert
# that the solver reproduces the analytic result documented in ReadMe.md, not
# merely that it is unchanged.

def _physics_freespace(intensity, md, m):
    # Gaussian beam diffraction: w(z) = w0*sqrt(1 + (z/zR)^2), zR = pi*w0^2*n0/lambda
    w0 = md["beam_w0"][0]
    lam = md["lambda"][0]
    n0 = md["n_0"][0]
    length = md["Zn"][-1] - md["Zn"][0]
    z_r = math.pi * w0 ** 2 * n0 / lam
    w_theory = w0 * math.sqrt(1.0 + (length / z_r) ** 2)
    return [("output beam width w(L) [um]", m["wx"] * 1e6, w_theory * 1e6, 0.03)]


def _physics_tilt(intensity, md, m):
    # Wide-angle BPM: 20 deg tilt -> transverse shift converges to z*tan(20),
    # vs paraxial z*sin(20). ReadMe: z*tan(20)=21.8um.
    tilt_deg = 20.0
    length = md["Zn"][-1] - md["Zn"][0]
    x0 = md["beam_x0"][0]
    disp_theory = length * math.tan(math.radians(tilt_deg))
    disp_meas = m["cx"] - x0
    return [("centroid shift = z*tan(20) [um]", disp_meas * 1e6, disp_theory * 1e6, 0.05)]


PHYSICS_CHECKS = {
    "freespace": _physics_freespace,
    "tilt_wideangle": _physics_tilt,
}


def compare(name, m, ref):
    """Return (ok, lines) comparing metrics m against reference entry ref."""
    tol = ref.get("_meta", {}).get("tol", {})
    rtol = tol.get("rtol", DEFAULT_RTOL)
    atol_m = tol.get("atol_m", DEFAULT_ATOL_M)
    golden = ref["samples"].get(name, {}).get("metrics")
    lines = []
    ok = True
    if golden is None:
        return False, [f"  no golden reference for sample '{name}'"]
    for key in METRIC_FIELDS:
        got = m[key]
        want = golden[key]
        atol = atol_m if key in METRES_METRICS else 0.0
        diff = abs(got - want)
        limit = atol + rtol * abs(want)
        passed = diff <= limit
        ok = ok and passed
        lines.append(
            f"  [{'OK' if passed else 'FAIL'}] {key:6s} got={got:+.6e} ref={want:+.6e} "
            f"diff={diff:.2e} limit={limit:.2e}"
        )
    return ok, lines


def run_physics(name, intensity, md, m):
    """Return (ok, lines) for the physics theory checks of a sample, if any."""
    check = PHYSICS_CHECKS.get(name)
    if check is None:
        return True, []
    lines = ["  physics (vs analytic theory):"]
    ok = True
    for label, meas, theory, rtol in check(intensity, md, m):
        rel = abs(meas - theory) / abs(theory) if theory else 0.0
        passed = rel <= rtol
        ok = ok and passed
        lines.append(
            f"    [{'OK' if passed else 'FAIL'}] {label}: meas={meas:.4f} "
            f"theory={theory:.4f} rel={rel:.2%} tol={rtol:.0%}"
        )
    return ok, lines


def main():
    ap = argparse.ArgumentParser(description="OpenBPM numerical regression check")
    ap.add_argument("--name", required=True, help="sample name (basename without .ofd)")
    ap.add_argument("--h5", required=True, help="path to time_series_data.h5")
    ap.add_argument("--ref", default=REF_DEFAULT, help="reference.json path")
    ap.add_argument("--update", action="store_true", help="write metrics into reference.json")
    args = ap.parse_args()

    intensity, md, _ = load(args.h5)
    m = metrics(intensity, md)

    if args.update:
        if os.path.exists(args.ref):
            with open(args.ref) as fh:
                ref = json.load(fh)
        else:
            ref = {"_meta": {"tol": {"rtol": DEFAULT_RTOL, "atol_m": DEFAULT_ATOL_M},
                             "omp_threads": 1}, "samples": {}}
        ref.setdefault("samples", {})[args.name] = {"metrics": m}
        with open(args.ref, "w") as fh:
            json.dump(ref, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"updated reference for '{args.name}' in {args.ref}")
        return 0

    with open(args.ref) as fh:
        ref = json.load(fh)

    print(f"== {args.name} ==")
    ok_metrics, lines = compare(args.name, m, ref)
    for line in lines:
        print(line)
    ok_physics, plines = run_physics(args.name, intensity, md, m)
    for line in plines:
        print(line)

    if ok_metrics and ok_physics:
        print(f"  PASS")
        return 0
    print(f"  REGRESSION DETECTED for '{args.name}'", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
