"""Machine-readable catalogue of all archived solutions.

Regenerates events/catalogue.csv from every events/<publicID>/solution.json
after each processed event — the CSV is always derived from the archived
sidecars (single source of truth), never edited by hand.

Event detection and original hypocentres (time, location, preliminary
magnitude, initial depth) are sourced from GeoNet; the GeoNet_M and
GeoNet_depth columns carry that original solution alongside our inverted
Mw and centroid depth (CD) so the revision is explicit per event.

Column conventions follow the published NZ regional CMT solutions CSV
(GeoNet/data moment-tensor/GeoNet_CMT_solutions.csv) where they overlap so
the two are directly comparable; note their (and our) MT elements are in
units of 1e20 dyne-cm.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import config

_FLAG_NAMES = {
    "min_stations": "few_stations",
    "vr_floor": "low_VR",
    "az_gap_ok": "wide_az_gap",
}


def _quality_flag(quality: dict) -> str:
    if quality.get("passed"):
        return "True"
    failed = [
        _FLAG_NAMES.get(name, name)
        for name, ok in quality.get("checks", {}).items() if not ok
    ]
    grade = quality.get("grade", "?")
    return ";".join(failed) if failed else f"grade_{grade}"


COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "GeoNet_M", "GeoNet_depth", "Mw", "Depth", "Mo", "NS", "DC", "CLVD", "VR",
    "Mxx", "Mxy", "Mxz", "Myy", "Myz", "Mzz",
    "band", "model", "quality_flag", "published",
]


def build_catalogue(events_dir: Path | None = None) -> Path | None:
    """Scan <events_dir>/*/solution.json -> <events_dir>/catalogue.csv.

    Returns the CSV path, or None if there are no solutions yet.
    """
    events_dir = events_dir or config.EVENTS_DIR
    solutions = sorted(events_dir.glob("*/solution.json"))
    if not solutions:
        return None

    published_ids: set[str] = set()
    if config.STATE_FILE.exists():
        state = json.loads(config.STATE_FILE.read_text())
        published_ids = {p["public_id"] for p in state.get("published", [])}

    rows = []
    for path in solutions:
        s = json.loads(path.read_text())
        ev, p = s["event"], s["preferred"]
        mt = s["preferred"]["tensor_dyne_cm"]  # XYZ basis, dyne-cm
        rows.append({
            "PublicID": ev["public_id"],
            "Date": ev["origin_time"],
            "Latitude": round(ev["latitude"], 4),
            "Longitude": round(ev["longitude"], 4),
            "strike1": round(p["plane1"]["strike"], 1),
            "dip1": round(p["plane1"]["dip"], 1),
            "rake1": round(p["plane1"]["rake"], 1),
            "strike2": round(p["plane2"]["strike"], 1),
            "dip2": round(p["plane2"]["dip"], 1),
            "rake2": round(p["plane2"]["rake"], 1),
            "GeoNet_M": round(ev["prelim_mag"], 2),
            "GeoNet_depth": round(ev["depth_km"], 1),
            "Mw": round(p["mw"], 2),
            # our centroid depth beside our Mw: a large Mw revision is
            # often explained by the depth revision visible in the adjacent
            # GeoNet_depth column
            "Depth": p["depth_km"],
            # Mo in dyne-cm; MT elements in 1e20 dyne-cm (GeoNet convention)
            "Mo": f"{p['m0_dyne_cm']:.3e}",
            "NS": s["quality"]["n_stations_used"],
            "DC": round(p["pdc"], 1),
            "CLVD": round(p["pclvd"], 1),
            "VR": round(p["vr"], 1),
            **{
                k.capitalize().replace("m", "M", 1): round(v / 1e20, 4)
                for k, v in mt.items()
            },
            "band": s.get("chosen_band", ""),
            "model": s["provenance"]["velocity_model"],
            # "True" when the solution passed, else the actual exit
            # reason(s) so a reader can see WHY it fell short at a glance
            "quality_flag": _quality_flag(s["quality"]),
            "published": ev["public_id"] in published_ids,
        })

    rows.sort(key=lambda r: r["Date"])
    out = events_dir / "catalogue.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in COLUMNS})
    print(f"catalogue: {len(rows)} solutions -> {out}")
    return out


if __name__ == "__main__":
    build_catalogue()
