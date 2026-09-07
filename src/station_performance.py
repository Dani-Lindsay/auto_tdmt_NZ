"""Per-station aggregate of station_ledger.jsonl — "which stations
consistently get picked or dropped", the year-one learning table.

Regenerated with the catalogue after each event. Per station:
  n_seen / n_used / use_rate     entered the pool / survived the loop
  med_snr                        median of its median component SNR
  med_station_vr                 median own VR when used
  med_shift_s, med_abs_shift_s   solved time shift (a consistent sign is a
                                 path/velocity-model anomaly)
  med_amp_ratio                  distance-corrected amplitude vs the network
                                 median (drift = response metadata problem)
  n_<reason class>               why it was left out
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import catalogue
import config

CLASSES = ("nodata", "snr", "amp", "not_selected", "fit", "shift", "abort",
           "other")
COLUMNS = (["station", "n_seen", "n_used", "use_rate", "med_snr",
            "med_station_vr", "med_shift_s", "med_abs_shift_s",
            "med_amp_ratio"] + [f"n_{c}" for c in CLASSES] + ["last_event"])


def build_station_performance(events_dir: Path | None = None) -> Path | None:
    events_dir = events_dir or config.EVENTS_DIR
    out_dir = catalogue.tables_dir(events_dir)
    ledger_path = out_dir / "station_ledger.jsonl"
    if not ledger_path.exists():
        catalogue.build_catalogue(events_dir)
    if not ledger_path.exists():
        return None
    acc: dict[str, dict] = defaultdict(
        lambda: {"seen": 0, "used": 0, "snr": [], "vr": [], "shift": [],
                 "amp": [], "last": "", **{f"n_{c}": 0 for c in CLASSES}})

    with open(ledger_path) as f:
        for line in f:
            rec = json.loads(line)
            for sid, st in rec["stations"].items():
                e = acc[".".join(sid.split(".")[:2])]
                e["seen"] += 1
                e["last"] = max(e["last"], rec["Date"])
                for key, col in (("snr", "snr_med"), ("amp", "amp_ratio")):
                    if st.get(col) is not None:
                        e[key].append(st[col])
                if st.get("used"):
                    e["used"] += 1
                    for key, col in (("vr", "own_vr"), ("shift", "shift_s")):
                        if st.get(col) is not None:
                            e[key].append(st[col])
                else:
                    cls = st.get("drop") if st.get("drop") in CLASSES else "other"
                    e[f"n_{cls}"] += 1

    def _med(xs, nd=1):
        return round(float(np.median(xs)), nd) if xs else ""

    out = out_dir / "station_performance.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for k in sorted(acc):
            e = acc[k]
            w.writerow({
                "station": k, "n_seen": e["seen"], "n_used": e["used"],
                "use_rate": round(e["used"] / e["seen"], 2) if e["seen"] else 0,
                "med_snr": _med(e["snr"]), "med_station_vr": _med(e["vr"]),
                "med_shift_s": _med(e["shift"]),
                "med_abs_shift_s": _med([abs(x) for x in e["shift"]]),
                "med_amp_ratio": _med(e["amp"], 2),
                **{f"n_{c}": e[f"n_{c}"] for c in CLASSES},
                "last_event": e["last"],
            })
    print(f"station performance: {len(acc)} stations -> {out}")
    return out


if __name__ == "__main__":
    build_station_performance()
