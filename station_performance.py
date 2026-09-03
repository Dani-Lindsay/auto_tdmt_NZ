"""Per-station performance ledger, aggregated from every archived
solution.json — the raw material for "we know this station does well".

Regenerated alongside the catalogue after each event and consumed by the
validation scripts. Collected per station:

- n_seen / n_used / use_rate: how often the station enters the pool and
  survives selection;
- med_station_vr: median individual variance reduction when used (fit
  quality track record);
- mean_dv_pct / med_abs_dv_pct: signed and absolute velocity-model
  deviation from solved zcor — a consistent sign is a path anomaly, a
  large scatter is an unreliable station;
- med_amp_ratio: median distance-corrected amplitude vs network median —
  drift flags response-metadata problems (e.g. NZ.RDHZ at ~100x, or a
  dead channel near 0x);
- n_snr_drop (SNR or peak/noise floor) / n_amp_outlier / n_eliminated: why it gets excluded.

A future selection prior can read this table directly; today it is the
audit trail.
"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

import config

COLUMNS = [
    "station", "n_seen", "n_used", "use_rate", "med_station_vr",
    "mean_dv_pct", "med_abs_dv_pct", "med_amp_ratio",
    "n_snr_drop", "n_amp_outlier", "n_eliminated", "last_event",
]


def build_station_performance(events_dir: Path | None = None) -> Path | None:
    events_dir = events_dir or config.EVENTS_DIR
    solutions = sorted(events_dir.glob("*/solution.json"))
    if not solutions:
        return None

    ledger: dict[str, dict] = defaultdict(
        lambda: {"used": [], "dv": [], "amp": [], "snr_drop": 0,
                 "amp_out": 0, "elim": 0, "seen": 0, "last": ""})

    def _key(sid: str) -> str:
        parts = sid.split(".")
        return ".".join(parts[:2]) if len(parts) >= 2 else sid

    for path in solutions:
        s = json.loads(path.read_text())
        date = s["event"]["origin_time"][:10]
        for r in s.get("stations_used", []):
            k = f"{r['network']}.{r['station']}"
            e = ledger[k]
            e["seen"] += 1
            e["last"] = max(e["last"], date)
            if "station_vr" in r:
                e["used"].append(r["station_vr"])
            if "zcor_s" in r and r.get("distance_km"):
                e["dv"].append(
                    r["zcor_s"] / (r["distance_km"]
                                   / config.GROUP_VELOCITY_KMS) * 100.0)
            if "amp_ratio" in r:
                e["amp"].append(r["amp_ratio"])
        for d in s.get("stations_dropped", []):
            k = _key(d["station"])
            e = ledger[k]
            e["seen"] += 1
            e["last"] = max(e["last"], date)
            reason = d.get("reason", "")
            if "amp_ratio" in d:
                e["amp"].append(d["amp_ratio"])
            if re.search(r"SNR .* <|peak/noise .* <", reason):
                e["snr_drop"] += 1
            elif "amplitude outlier" in reason:
                e["amp_out"] += 1
            elif "eliminated" in reason or "test-drop" in reason:
                e["elim"] += 1

    out = events_dir / "station_performance.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for k in sorted(ledger):
            e = ledger[k]
            n_used = len(e["used"])
            w.writerow({
                "station": k,
                "n_seen": e["seen"],
                "n_used": n_used,
                "use_rate": round(n_used / e["seen"], 2) if e["seen"] else 0,
                "med_station_vr": (round(float(np.median(e["used"])), 1)
                                   if e["used"] else ""),
                "mean_dv_pct": (round(float(np.mean(e["dv"])), 1)
                                if e["dv"] else ""),
                "med_abs_dv_pct": (round(float(np.median(np.abs(e["dv"]))), 1)
                                   if e["dv"] else ""),
                "med_amp_ratio": (round(float(np.median(e["amp"])), 2)
                                  if e["amp"] else ""),
                "n_snr_drop": e["snr_drop"],
                "n_amp_outlier": e["amp_out"],
                "n_eliminated": e["elim"],
                "last_event": e["last"],
            })
    print(f"station performance: {len(ledger)} stations -> {out}")
    return out


if __name__ == "__main__":
    build_station_performance()
