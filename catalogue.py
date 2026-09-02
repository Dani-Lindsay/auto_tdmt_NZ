"""Machine-readable catalogue of all archived solutions.

Regenerates events/catalogue.csv from every events/<publicID>/solution.json
after each processed event — the CSV is always derived from the archived
sidecars (single source of truth), never edited by hand.

Column conventions follow GeoNet's manual CMT catalogue
(GeoNet/data moment-tensor/GeoNet_CMT_solutions.csv) where they overlap so
the two are directly comparable; note their (and our) MT elements are in
units of 1e20 dyne-cm.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import config

COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "Mprelim", "Mw", "Mo", "CD", "NS", "DC", "CLVD", "VR",
    "Mxx", "Mxy", "Mxz", "Myy", "Myz", "Mzz",
    "band", "model", "gates_passed", "published",
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
            "Mprelim": round(ev["prelim_mag"], 2),
            "Mw": round(p["mw"], 2),
            # Mo in dyne-cm; MT elements in 1e20 dyne-cm (GeoNet convention)
            "Mo": f"{p['m0_dyne_cm']:.3e}",
            "CD": p["depth_km"],
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
            "gates_passed": s["quality"]["passed"],
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
