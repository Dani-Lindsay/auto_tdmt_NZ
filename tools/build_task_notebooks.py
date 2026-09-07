"""Generate docs/task_1_geonet.ipynb ... task_5_publish.ipynb.

Each notebook walks one task by IMPORTING the task module and calling
its functions on a real event — no logic is copied, so the notebooks
cannot drift from the code. Rebuild after any change:

    pixi run python tools/build_task_notebooks.py
    pixi run python tools/build_task_notebooks.py --execute   # run them too

Notebooks are committed with outputs stripped (the --execute run is a
smoke test, not a deliverable).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
EVENT = "2026p669681"   # Mw 3.7 Opunake: small, fast, a B under v5


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {},
            "source": text.strip("\n").splitlines(keepends=True)}


def code(text: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(keepends=True)}


PREAMBLE = f'''
import os, sys, json
from pathlib import Path
ROOT = Path.cwd().parent if Path.cwd().name == "docs" else Path.cwd()
sys.path.insert(0, str(ROOT / "src"))
os.chdir(ROOT)
# work in a scratch archive so the real events/ are untouched
os.environ.setdefault("AUTO_TDMT_EVENTS", str(Path.home() / "work" / "proj_tdmt_NZ" / "notebook_runs"))
import config
from config import P            # every tunable, from auto_tdmt.cfg
EVENT = "{EVENT}"
print("parameters from", P.source)
'''

NOTEBOOKS = {
    "task_1_geonet": [
        md("""
# Task 1 — GeoNet search

What the pipeline does first: ask GeoNet for recent events, apply the
processing floor from `auto_tdmt.cfg` §1, and fetch one event's origin.
Every function here is the one the pipeline itself runs
(`geonet.py`, `trigger.py`, `watch.py`).
"""),
        code(PREAMBLE),
        code("""
print(json.dumps(P.as_dict()["geonet"], indent=2))
"""),
        md("## 1.1 Recent events above the floor"),
        code("""
import trigger
from geonet import recent_quakes
events = recent_quakes(mmi=3)                 # one request to the quake API
print(f"{len(events)} recent events from the quake API")
for ev in events[:15]:
    ok, why = trigger.passes_processing_floor(ev)
    print(f"{ev.public_id} M{ev.prelim_mag:.1f} {ev.depth_km:5.1f} km "
          f"{ev.locality[:32]:32s} -> {'PROCESS' if ok else why}")
"""),
        md("## 1.2 One event's origin (the input to task 2)"),
        code("""
from geonet import get_event
ev = get_event(EVENT)
print(ev)
print(json.dumps(ev.to_dict(), indent=2))
"""),
        md("""
## What to check
- Is the floor doing what you expect (`geonet.minPrelimMag`, `maxDepthKm`)?
- Depths reported as 5 / 12 / 33 km are fixed default values rather than
  fitted depths — those events are always attempted and the depth search
  decides (task 3).
"""),
    ],
    "task_2_stations": [
        md("""
# Task 2 — station selection and waveform processing

`waveforms.fetch_and_process` does everything in this task. The rules and
their sources are in the module docstring and `auto_tdmt.cfg` §2:
distance window by magnitude; unusable data rejected; per-component SNR;
amplitude outliers; azimuth-balanced round-robin selection.
"""),
        code(PREAMBLE),
        code("""
print(json.dumps(P.as_dict()["station"], indent=2))
"""),
        md("## 2.1 Run the task on one event and one band"),
        code("""
import waveforms
from geonet import get_event
ev = get_event(EVENT)
band = config.band_candidates(ev.prelim_mag)[0]      # first band of the menu
wd = config.EVENTS_DIR / ev.public_id / "task2"
pool, dropped = waveforms.fetch_and_process(ev, wd, band)
print(f"band {config.band_tag(band)}: {len(pool)} selected, {len(dropped)} not")
"""),
        md("## 2.2 What was selected, and why the rest were not"),
        code("""
print(f"{'station':14s} {'dist':>5s} {'az':>4s} sec {'SNR Z/R/T':>17s} {'med':>5s} {'win':>4s} {'amp':>5s}")
for r in pool:
    print(f"{waveforms.station_id(r):14s} {r['distance_km']:5.0f} {r['azimuth']:4.0f} "
          f"{r['sector']:3d} {r['snr']['Z']:5.1f}/{r['snr']['R']:5.1f}/{r['snr']['T']:5.1f} "
          f"{r['snr_med']:5.1f} {r['window_end_s']:4d} {r.get('amp_ratio', float('nan')):5.2f}")
print()
for d in dropped:
    print(f"{d['station']:14s} {d['reason']}")
"""),
        md("## 2.3 The traces the inversion will see"),
        code("""
from obspy import read
import matplotlib.pyplot as plt
fig, axes = plt.subplots(len(pool), 1, figsize=(10, 1.2 * len(pool)), sharex=True)
for ax, r in zip(axes, pool):
    for comp, c in zip("ZRT", ("k", "#0072B2", "#E69F00")):
        tr = read(str(wd / f"{waveforms.station_id(r)}.{comp}.dat"), format="SAC")[0]
        t = tr.times() + tr.stats.sac.b
        ax.plot(t, tr.data, c, lw=0.6, label=comp)
    ax.axvline(r["window_end_s"], color="0.5", ls="--", lw=0.8)
    ax.set_ylabel(r["station"], rotation=0, ha="right", fontsize=8)
    ax.set_yticks([])
axes[0].legend(ncol=3, fontsize=7, loc="upper right")
axes[-1].set_xlabel("s after origin  (dashed = end of the inverted window)")
plt.tight_layout(); plt.show()
"""),
        md("""
## What to check
- Does the SNR threshold agree with your eye? (`station.snrMin`)
- Is the inverted window (dashed) long enough for the far stations?
  (`station.maxWindowS`, `groupVelKms`, `windowTailS`)
- Did the round-robin leave out a station you would have kept?
  (`station.maxStations`)
"""),
    ],
    "task_3_invert": [
        md("""
# Task 3 — moment tensor inversion

mttime does the inversion. `invert.py` writes its control file, bounds
the depth grid around the GeoNet hypocentre, runs the Clinton loop
(drop stations the solution cannot explain, re-invert), and grades the
result with the INGV table. Sources in the module docstring and
`auto_tdmt.cfg` §3.
"""),
        code(PREAMBLE),
        code("""
print(json.dumps(P.as_dict()["invert"], indent=2))
"""),
        md("## 3.1 Stations from task 2, Green's functions staged"),
        code("""
import greens, invert, waveforms
from geonet import get_event
ev = get_event(EVENT)
band = config.band_candidates(ev.prelim_mag)[0]
wd = config.EVENTS_DIR / ev.public_id / "task3"
pool, dropped = waveforms.fetch_and_process(ev, wd, band)
model = config.model_for_event(ev.latitude, ev.longitude)
depths = invert.search_depths(ev, model)
print(f"{len(pool)} stations; model {model}; depth grid "
      f"{depths[0]:g}-{depths[-1]:g} km ({len(depths)} depths) around GeoNet {ev.depth_km:g} km")
greens.stage_event_greens(model, pool, depths, band, wd / "greens")
"""),
        md("## 3.2 The Clinton loop"),
        code("""
os.chdir(wd)
inv, used, rejected, rounds = invert.clinton_loop(ev, pool, depths, wd, wd / "greens")
os.chdir(ROOT)
print(json.dumps(rounds, indent=1))
for r in rejected:
    print(f"dropped {r['station']}: {r['reason']}")
"""),
        md("## 3.3 The solution, its depth curve and its grade"),
        code("""
sol = invert.summarize(inv, ev, used, dropped + rejected, model, rounds)
os.chdir(wd)
sol["jackknife"] = invert.jackknife(ev, used, sol["preferred"]["depth_km"], wd, wd / "greens", sol["preferred"]["plane1"])
os.chdir(ROOT)
sol["quality"] = invert.quality_gates(sol)
p, q = sol["preferred"], sol["quality"]
print(f"Mw {p['mw']:.2f}  depth {p['depth_km']:g} km [{q['depth_range_km'][0]:g}-{q['depth_range_km'][1]:g}]  "
      f"VR {p['vr']:.1f}  DC {p['pdc']:.0f}  grade {q['grade']}  ({q['n_stations_used']} stations)")
print("plane1", p["plane1"]); print("plane2", p["plane2"])
print("checks", q["checks"]); print("jackknife", {k: v for k, v in sol["jackknife"].items() if k != "subsets"})
import matplotlib.pyplot as plt
d = sorted(sol["depth_search"], key=lambda r: r["depth_km"])
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot([r["depth_km"] for r in d], [r["vr"] for r in d], "k.-", label="VR")
ax.axvspan(*q["depth_range_km"], color="0.85", label=f"within {P.invert.depthUncPct:g}% of max")
ax.axvline(ev.depth_km, color="#D55E00", ls="--", label="GeoNet")
ax2 = ax.twinx(); ax2.plot([r["depth_km"] for r in d], [r["pdc"] for r in d], "-", color="#0072B2", alpha=0.6, label="%DC")
ax.set_xlabel("depth (km)"); ax.set_ylabel("VR (%)"); ax2.set_ylabel("%DC", color="#0072B2")
ax.legend(loc="lower right", fontsize=8); plt.tight_layout(); plt.show()
"""),
        md("## 3.4 Per-station fit"),
        code("""
print(f"{'station':14s} {'dist':>5s} {'az':>4s} {'own VR':>7s} {'shift s':>8s}")
for r in used:
    print(f"{invert.station_id(r):14s} {r['distance_km']:5.0f} {r['azimuth']:4.0f} {r['station_vr']:7.1f} {r['zcor_s']:+8.1f}")
"""),
        md("""
## What to check
- Do the dropped stations deserve it? (`invert.stationVRDrop`, `stationVRFloor`, `maxTimeShiftS`)
- Is the depth curve a plateau or a spike? Does the shaded range include GeoNet's depth?
  (`invert.depthWindowKm`, `depthUncPct`)
- Would you publish this? Compare your call with the grade (`invert.gradeA_VR` ... `dcMinPublish`).
"""),
    ],
    "task_4_forward": [
        md("""
# Task 4 — forward model of surface displacement

For each nodal plane, a rectangular uniform-slip Okada (1992) dislocation
via okada4py, with fault dimensions from Wells & Coppersmith (1994).
This is the quantity the pipeline exists to produce: is this event
InSAR-relevant? Parameters in `auto_tdmt.cfg` §4.
"""),
        code(PREAMBLE),
        code("""
print(json.dumps(P.as_dict()["forward"], indent=2))
"""),
        md("## 4.1 From an archived solution"),
        code("""
import okada_forward
cands = sorted(Path(config.REPO_DIR, "events").glob(f"{EVENT}*/solution.json")) or \\
        sorted(config.EVENTS_DIR.glob(f"{EVENT}*/solution.json"))
sol = json.loads(cands[0].read_text())
fwd = okada_forward.forward_both_planes(sol)
print(f"peak |u| {fwd['peak_abs_m']*100:.3f} cm  detectable: {fwd['detectable']}")
for pl in ("plane1", "plane2"):
    print(pl, fwd[pl]["fault"])
"""),
        md("## 4.2 The displacement field"),
        code("""
import matplotlib.pyplot as plt
import numpy as np
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax, pl in zip(axes, ("plane1", "plane2")):
    g = fwd[pl]
    im = ax.pcolormesh(g["x_km"], g["y_km"], g["uz_m"] * 100, cmap="RdBu_r",
                       vmin=-abs(g["uz_m"]).max()*100, vmax=abs(g["uz_m"]).max()*100, shading="auto")
    ax.set_title(f"{pl}: strike {g['fault']['strike']:.0f} dip {g['fault']['dip']:.0f} rake {g['fault']['rake']:.0f}")
    ax.set_aspect("equal"); ax.set_xlabel("east (km)"); ax.set_ylabel("north (km)")
    plt.colorbar(im, ax=ax, label="vertical (cm)")
plt.tight_layout(); plt.show()
"""),
        md("""
## What to check
- Rigidity, the Wells & Coppersmith coefficients and the grid are all in the cfg.
- The publish gate uses `peak_abs_m` against `publish.minDisplacementM`.
"""),
    ],
    "task_5_publish": [
        md("""
# Task 5 — archive and publish decision

Every processed event is archived; the email goes out only when the
gates in `auto_tdmt.cfg` §5 pass. The three tables are rebuilt from the
archived `solution.json` files by `catalogue.py`.
"""),
        code(PREAMBLE),
        code("""
print(json.dumps(P.as_dict()["publish"], indent=2))
"""),
        md("## 5.1 The publish decision for an archived solution"),
        code("""
import trigger, okada_forward, publish
cands = sorted(Path(config.REPO_DIR, "events").glob(f"{EVENT}*/solution.json"))
sol = json.loads(cands[0].read_text())
fwd = okada_forward.forward_both_planes(sol)
decision = trigger.publish_decision(sol, fwd, published_history=[])
print(json.dumps(decision, indent=2))
subject, body = publish.draft_text(sol, fwd, sol.get("nisar_passes", []))
print(subject); print(body[:1500])
"""),
        md("## 5.2 The tables"),
        code("""
import pandas as pd
cat = pd.read_csv(config.REPO_DIR / "catalogue.csv")
print(cat["Grade"].value_counts().sort_index().to_dict(), "of", len(cat), "events")
display(cat[cat.PublicID == EVENT].T)
print("publish flags:", cat["publish_flag"].value_counts().head(8).to_dict())
npub = pd.read_csv(config.REPO_DIR / "not_published.csv")
print(f"{len(npub)} events with no solution; by outcome and stage:")
print(npub.groupby(["Outcome", "Stage"]).size())
ledger = pd.read_json(config.REPO_DIR / "station_ledger.jsonl", lines=True)
rec = ledger[ledger.PublicID == EVENT].iloc[0]
display(pd.DataFrame(rec["stations"]).T)   # one row per station, numbers only
"""),
        md("""
## What to check
- `publish_flag` in `catalogue.csv` says why a solved event did not email.
- `not_published.csv` (repo root) lists the events with no solution at all.
- `station_ledger.jsonl` is the year-one learning table (one line per
  event, a per-station dictionary of numbers): after a year,
  `station_performance.csv` (its aggregate) shows which stations are
  consistently picked or dropped, and by which rule.
"""),
    ],
}


def build(execute: bool) -> None:
    DOCS.mkdir(exist_ok=True)
    for name, cells in NOTEBOOKS.items():
        nb = {"cells": cells, "metadata": {
            "kernelspec": {"display_name": "Python (auto_tdmt_NZ pixi)",
                           "language": "python", "name": "python3"},
            "language_info": {"name": "python"}},
              "nbformat": 4, "nbformat_minor": 5}
        out = DOCS / f"{name}.ipynb"
        out.write_text(json.dumps(nb, indent=1))
        print(f"wrote {out.relative_to(ROOT)}")
        if execute:
            r = subprocess.run(
                [sys.executable, "-m", "jupyter", "nbconvert", "--to",
                 "notebook", "--execute", "--inplace",
                 "--ExecutePreprocessor.timeout=1800", str(out)],
                cwd=ROOT, capture_output=True, text=True)
            print("  executed OK" if r.returncode == 0
                  else f"  FAILED:\n{r.stderr[-2000:]}")
            # strip outputs again so the committed file is clean
            nb = json.loads(out.read_text())
            for c in nb["cells"]:
                if c["cell_type"] == "code":
                    c["outputs"], c["execution_count"] = [], None
            out.write_text(json.dumps(nb, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()
    build(args.execute)
