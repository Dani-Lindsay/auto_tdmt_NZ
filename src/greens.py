"""Green's functions — attribution and provenance.

This module was compiled with Claude (Anthropic) assistance. The CPS
computation chain (hprep96 | hspec96 | hpulse96 -D -i | f96tosac -B) is
taken from mttime example notebook 02 by Andrea Chiang (LLNL):
https://github.com/LLNL/mttime/tree/master/examples/notebooks, using
Computer Programs in Seismology 3.30 by Robert Herrmann:
https://rbherrmann.github.io/ComputerProgramsSeismology/ (see also the CPS
Green's-function tutorial at
https://www.eas.slu.edu/eqc/eqc_cps/TUTORIAL/GREEN/index.html).

Deviations from the notebook: a precomputed distance x depth library per
velocity model (built once, reused by nearest-distance lookup) instead of
per-event GF computation; unfiltered storage with per-event band-pass;
manifest with model SHA-256 for provenance.

Green's function library: one-time CPS build (local) + runtime lookup.

Build (local only, needs CPS binaries — never runs in CI):
    pixi run python greens.py --build
    pixi run python greens.py --build --depths 4 6 8   # subset first

Library layout:
    <GF_LIBRARY_DIR>/<model>/<version>/d<depth km, %06.2f>/<dist km, %04d>.<GRN>
with the ten Herrmann fundamental sources GRN in
    ZDD RDD ZDS RDS TDS ZSS RSS TSS ZEX REX
stored UNFILTERED (raw hpulse96 -D -i output, 1e20 dyne-cm source, cm).
Per-event code filters copies to the event band — the library itself is
band-agnostic.

The CPS chain per depth (mttime example notebook 02, verbatim):
    hprep96 -M <model> -d dfile -HS <depth> -HR 0 -EQEX
    hspec96
    hpulse96 -D -i > file96
    f96tosac -B file96
which yields B%03d%02d<GRN>.sac numbered by (distance line, source index).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from obspy import read

import config

# The ten Herrmann & Wang (1985) fundamental-source time series. This set,
# its names, and its ORDER are a fixed convention shared by CPS (f96tosac -B
# writes B<dist><j><name>.sac with j=01..10 in exactly this sequence) and by
# mttime's green="herrmann" reader — not a tunable parameter. Z/R/T =
# vertical/radial/transverse; DD = 45-deg dip-slip, DS = 90-deg dip-slip,
# SS = vertical strike-slip, EX = explosion. TDD/TEX do not exist because
# those sources are azimuthally symmetric (no transverse motion). If CPS
# ever changed the convention, build_depth's file-existence assert fails
# loudly rather than mis-assigning a component.
GREENS = ("ZDD", "RDD", "ZDS", "RDS", "TDS", "ZSS", "RSS", "TSS", "ZEX", "REX")


def model_path(model: str) -> Path:
    p = config.REPO_DIR / "models" / f"{model}.d"
    assert p.exists(), f"velocity model missing: {p}"
    return p


def library_root(model: str) -> Path:
    return config.GF_LIBRARY_DIR / model / config.GF_VERSION


def depth_dir(model: str, depth_km: float) -> Path:
    return library_root(model) / f"d{depth_km:06.2f}"


def gf_file(model: str, depth_km: float, dist_km: int, grn: str) -> Path:
    return depth_dir(model, depth_km) / f"{dist_km:04d}.{grn}"


def nearest_grid_distance(dist_km: float) -> int:
    grid = config.GF_DIST_KM
    assert grid[0] <= dist_km <= grid[-1], (
        f"distance {dist_km:.1f} km outside GF grid {grid[0]}-{grid[-1]} km"
    )
    return min(grid, key=lambda d: abs(d - dist_km))


def available_depths(model: str) -> list[float]:
    root = library_root(model)
    if not root.exists():
        return []
    return sorted(
        float(p.name[1:]) for p in root.iterdir() if p.name.startswith("d")
    )


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def _run(cmd: str, cwd: Path, env: dict) -> None:
    r = subprocess.run(
        cmd, shell=True, cwd=cwd, env=env,
        capture_output=True, text=True, timeout=3600,
    )
    assert r.returncode == 0, f"'{cmd}' failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"


def build_depth(model: str, depth_km: float, workdir: Path, env: dict) -> None:
    """Run the CPS chain for one source depth over the full distance grid."""
    workdir.mkdir(parents=True, exist_ok=True)
    with open(workdir / "dfile", "w") as f:
        for d in config.GF_DIST_KM:
            # dist dt npts t0 vred ; vred=0 + t0=0 -> traces start at origin
            f.write(f"{d:.0f} {config.DT:.2f} {config.GF_NPTS:d} 0 0.0\n")
    # CPS Fortran arg parsing rejects long absolute paths -> relative only
    shutil.copy(model_path(model), workdir / "model.d")

    _run(
        f"hprep96 -M model.d -d dfile -HS {depth_km:.4f} -HR 0 -EQEX",
        workdir, env,
    )
    _run("hspec96", workdir, env)
    _run("hpulse96 -D -i > file96", workdir, env)
    _run("f96tosac -B file96", workdir, env)

    out = depth_dir(model, depth_km)
    out.mkdir(parents=True, exist_ok=True)
    for i, dist in enumerate(config.GF_DIST_KM):
        for j, grn in enumerate(GREENS):
            src = workdir / f"B{i + 1:03d}{j + 1:02d}{grn}.sac"
            assert src.exists(), f"CPS did not produce {src.name}"
            tr = read(str(src), format="SAC")[0]
            assert tr.stats.npts == config.GF_NPTS, (
                f"{src.name}: npts {tr.stats.npts} != {config.GF_NPTS}"
            )
            assert abs(tr.stats.delta - config.DT) < 1e-6, (
                f"{src.name}: dt {tr.stats.delta} != {config.DT}"
            )
            shutil.move(str(src), str(gf_file(model, depth_km, dist, grn)))
    for leftover in workdir.glob("B*.sac"):
        leftover.unlink()


def build_library(model: str, depths: list[float]) -> None:
    assert config.CPS_BIN.exists(), (
        f"CPS bin dir not found: {config.CPS_BIN} (build CPS first)"
    )
    env = dict(os.environ, PATH=f"{config.CPS_BIN}:{os.environ['PATH']}")
    workdir = library_root(model) / "_build_tmp"

    t_start = time.time()
    for depth in depths:
        t0 = time.time()
        build_depth(model, depth, workdir, env)
        print(f"depth {depth:6.2f} km done in {time.time() - t0:6.1f} s", flush=True)

    manifest = {
        "model": model,
        "model_sha256": hashlib.sha256(model_path(model).read_bytes()).hexdigest(),
        "version": config.GF_VERSION,
        "distances_km": config.GF_DIST_KM,
        "depths_km": available_depths(model),
        "dt_s": config.DT,
        "npts": config.GF_NPTS,
        "cps_chain": "hprep96 -EQEX | hspec96 | hpulse96 -D -i | f96tosac -B",
        "units": "cm displacement for a 1e20 dyne-cm source (mttime convention)",
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "build_seconds": round(time.time() - t_start, 1),
    }
    with open(library_root(model) / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    shutil.rmtree(workdir, ignore_errors=True)
    print(f"library at {library_root(model)} with "
          f"{len(available_depths(model))} depths")


# ---------------------------------------------------------------------------
# Runtime: stage filtered GFs for one event
# ---------------------------------------------------------------------------

def stage_event_greens(
    model: str, stations: list[dict], depths: list[float],
    band_hz: tuple[float, float], green_dir: Path,
) -> None:
    """Copy + filter library GFs into <green_dir> using mttime naming:
    <NET.STA.LOC>.<depth %.4f>.<GRN>, filtered like the data."""
    have = set(available_depths(model))
    missing = [d for d in depths if d not in have]
    assert not missing, f"depths {missing} not in GF library (have {sorted(have)})"

    green_dir.mkdir(parents=True, exist_ok=True)
    fmin, fmax = band_hz
    for row in stations:
        sid = f"{row['network']}.{row['station']}.{row['location']}"
        dist = nearest_grid_distance(row["distance_km"])
        row["gf_distance_km"] = dist
        for depth in depths:
            for grn in GREENS:
                st = read(str(gf_file(model, depth, dist, grn)), format="SAC")
                st.filter(
                    "bandpass", freqmin=fmin, freqmax=fmax,
                    corners=config.FILTER_CORNERS, zerophase=True,
                )
                st.write(
                    str(green_dir / f"{sid}.{depth:.4f}.{grn}"), format="SAC"
                )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build", action="store_true")
    ap.add_argument(
        "--models", nargs="*", default=list(config.GF_MODELS),
        help="velocity model name(s); default = all in config.GF_MODELS",
    )
    ap.add_argument(
        "--depths", nargs="*", type=float, default=None,
        help="subset of depths (km); default = full config grid",
    )
    args = ap.parse_args()
    if args.build:
        for m in args.models:
            build_library(m, [float(d) for d in (args.depths or config.GF_DEPTHS_KM)])
    else:
        ap.print_help()
