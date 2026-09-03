"""Rebuild events_human/catalogue_human.csv from every reviewed
solution.json — the CSV is a BUILD PRODUCT, never hand-edited or
hand-merged, so student Pull Requests can never conflict on it.

Layouts supported:
  events_human/<event_dir>/<reviewer-slug>/solution.json   (current)
  events_human/<event_dir>/solution.json                   (legacy)

Standard library only (runs on a bare GitHub Actions runner):

    python3 catalogue_human.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HUMAN = ROOT / "events_human"

COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "GeoNet_M", "GeoNet_depth", "Mw", "Depth", "NS", "AzGap",
    "Grade", "DC", "VR", "Band", "Model",
    "Reviewer", "ReviewDate", "Decision", "Changes",
    "Auto_Mw", "Auto_Depth", "Auto_Grade",
]


def band_tag(band_hz) -> str:
    try:
        lo, hi = float(band_hz[0]), float(band_hz[1])
        return f"{round(1 / hi)}-{round(1 / lo)}s"
    except Exception:  # noqa: BLE001 - band is optional provenance
        return ""


def row_from(sol: dict) -> dict:
    ev = sol["event"]
    p = sol["preferred"]
    q = sol["quality"]
    hr = sol.get("human_review", {})
    ref = sol.get("automated_reference", {})
    prov = sol.get("provenance", {})
    return {
        "PublicID": ev["public_id"],
        "Date": ev["origin_time"][:10],
        "Latitude": ev["latitude"],
        "Longitude": ev["longitude"],
        "strike1": round(p["plane1"]["strike"]),
        "dip1": round(p["plane1"]["dip"]),
        "rake1": round(p["plane1"]["rake"]),
        "strike2": round(p["plane2"]["strike"]),
        "dip2": round(p["plane2"]["dip"]),
        "rake2": round(p["plane2"]["rake"]),
        "GeoNet_M": round(ev.get("prelim_mag", float("nan")), 2),
        "GeoNet_depth": round(ev.get("depth_km", float("nan")), 1),
        "Mw": round(p["mw"], 2),
        "Depth": p["depth_km"],
        "NS": q["n_stations_used"],
        "AzGap": round(q["azimuthal_gap_deg"]),
        "Grade": q["grade"],
        "DC": round(p["pdc"]),
        "VR": round(p["vr"], 1),
        "Band": band_tag(sol.get("filter_band_hz", ())),
        "Model": prov.get("velocity_model", sol.get("model", "")),
        "Reviewer": hr.get("reviewer", ""),
        "ReviewDate": hr.get("date", ""),
        "Decision": hr.get("decision", ""),
        "Changes": hr.get("changes", ""),
        "Auto_Mw": ref.get("mw", ""),
        "Auto_Depth": ref.get("depth_km", ""),
        "Auto_Grade": ref.get("grade", ""),
    }


def build() -> Path:
    rows = []
    for sol_path in sorted(HUMAN.glob("*/solution.json")) \
            + sorted(HUMAN.glob("*/*/solution.json")):
        try:
            rows.append(row_from(json.loads(sol_path.read_text())))
        except Exception as e:  # noqa: BLE001 - name the bad file, keep going
            print(f"SKIPPED {sol_path}: {e}")
    rows.sort(key=lambda r: (r["PublicID"], r["Reviewer"]))
    out = HUMAN / "catalogue_human.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} reviews)")
    return out


if __name__ == "__main__":
    build()
