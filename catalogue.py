"""Task 5 — the archive tables, regenerated from every solution.json.

    events/catalogue.csv          one row per event (solved or not)
    events/not_published.csv      every processed event that did NOT email, and why
    events/station_ledger.csv     one row per station per event: the year-one
                                  learning table ("which stations get picked")

All three are BUILD PRODUCTS derived from the archived sidecars, never
edited by hand. Column conventions follow the published NZ regional CMT
CSV (GeoNet/data moment-tensor) where they overlap; MT elements are in
1e20 dyne-cm like theirs. The human-review catalogue reuses COLUMNS so
the two catalogues line up column for column.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import config

_FLAG_NAMES = {
    "vr_grade": "VR_below_grade_B",
    "dc_floor": "low_DC",
    "jackknife_stable": "unstable_mechanism",
    "min_stations": "few_stations",
    # older vintages, kept so archived rows still read sensibly
    "vr_floor": "low_VR", "no_passengers": "station_not_fitting",
    "az_pair_90": "no_90deg_azimuth_pair", "depth_interior": "grid_edge_depth",
    "depth_agrees_with_geonet": "depth_far_from_geonet",
    "az_gap_ok": "wide_az_gap",
}


def _quality_flag(quality: dict) -> str:
    if quality.get("passed"):
        return "True"
    failed = [_FLAG_NAMES.get(k, k)
              for k, ok in quality.get("checks", {}).items() if not ok]
    return ";".join(failed) if failed else f"grade_{quality.get('grade', '?')}"


def _publish_flag(decision: dict) -> str:
    if decision.get("publish"):
        return "True"
    tags = decision.get("reason_tags")
    return ";".join(tags) if tags else "unpublished"


COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "GeoNet_M", "GeoNet_depth", "GeoNet_depth_unc", "Mw", "Depth",
    "Depth_lo", "Depth_hi", "dZ_GeoNet", "Depth_DCmax", "Mo",
    "NS", "AzGap", "Grade", "DC", "CLVD", "VR", "MinStaVR",
    "Jk_n", "Jk_Mw_std", "Jk_DC_std", "Jk_rot_deg",
    "PredDisp_cm", "Detectable",
    "Mxx", "Mxy", "Mxz", "Myy", "Myz", "Mzz",
    "Band", "Model", "Status", "Selection", "Code",
    "quality_flag", "publish_flag", "published",
]

LEDGER_COLUMNS = [
    "PublicID", "Date", "Station", "Distance_km", "Azimuth", "Sector",
    "SNR_Z", "SNR_R", "SNR_T", "SNR_med", "Amp_ratio", "Window_s",
    "Used", "Reason_class", "Reason", "Station_VR", "Shift_s",
    "Band", "Depth", "Mw", "VR", "Grade", "Selection",
]


def row_from(s: dict, published_ids: set[str] | None = None) -> dict:
    """Catalogue row for one solution.json (shared with catalogue_human)."""
    ev = s["event"]
    prov = s.get("provenance", {})
    base = {
        "PublicID": ev["public_id"], "Date": ev["origin_time"],
        "Latitude": round(ev["latitude"], 4),
        "Longitude": round(ev["longitude"], 4),
        "GeoNet_M": round(ev["prelim_mag"], 2),
        "GeoNet_depth": round(ev["depth_km"], 1),
        "GeoNet_depth_unc": ev.get("depth_unc_km", ""),
        "Band": s.get("chosen_band", "").replace("band_", ""),
        "Model": prov.get("velocity_model", ""),
        "Status": s.get("status", config.STATUS_SOLVED),
        "Selection": prov.get("selection_version", "v3"),
        "Code": prov.get("code_commit", ""),
        "published": (ev["public_id"] in published_ids
                      if published_ids is not None else ""),
    }
    if not config.is_solved(s):
        return {**base, "NS": 0, "Grade": "X",
                "quality_flag": f"no_coherent_solution:{s['abort']['stage']}",
                "publish_flag": "no_coherent_solution"}
    p, q = s["preferred"], s["quality"]
    mt = p["tensor_dyne_cm"]  # XYZ basis, dyne-cm
    jk = s.get("jackknife", {}) or {}
    lo, hi = (s.get("depth_pick_flags", {}).get("depth_range_km")
              or q.get("depth_range_km") or ("", ""))
    return {
        **base,
        "strike1": round(p["plane1"]["strike"], 1),
        "dip1": round(p["plane1"]["dip"], 1),
        "rake1": round(p["plane1"]["rake"], 1),
        "strike2": round(p["plane2"]["strike"], 1),
        "dip2": round(p["plane2"]["dip"], 1),
        "rake2": round(p["plane2"]["rake"], 1),
        "Mw": round(p["mw"], 2), "Depth": p["depth_km"],
        "Depth_lo": lo, "Depth_hi": hi,
        "dZ_GeoNet": round(p["depth_km"] - ev["depth_km"], 1),
        "Depth_DCmax": s.get("depth_pick_flags", {}).get("dc_max_depth_km", ""),
        "Mo": f"{p['m0_dyne_cm']:.3e}",
        "NS": q["n_stations_used"], "AzGap": q.get("azimuthal_gap_deg", ""),
        "Grade": q.get("grade", ""), "DC": round(p["pdc"], 1),
        "CLVD": round(p["pclvd"], 1), "VR": round(p["vr"], 1),
        "MinStaVR": q.get("min_own_vr", ""),
        "Jk_n": jk.get("n_subsets", ""), "Jk_Mw_std": jk.get("mw_std", ""),
        "Jk_DC_std": jk.get("dc_std", ""),
        "Jk_rot_deg": jk.get("max_tensor_rotation_deg", ""),
        "PredDisp_cm": (round(s["forward_model"]["peak_abs_cm"], 3)
                        if s.get("forward_model") else ""),
        "Detectable": s.get("forward_model", {}).get("detectable", ""),
        **{k.capitalize().replace("m", "M", 1): round(v / 1e20, 4)
           for k, v in mt.items()},
        "quality_flag": _quality_flag(q),
        "publish_flag": _publish_flag(s.get("publish_decision", {})),
    }


def ledger_rows(s: dict) -> list[dict]:
    """One row per station that entered this event, used or not."""
    import invert
    ev = s["event"]
    solved = config.is_solved(s)
    p = s.get("preferred", {})
    common = {
        "PublicID": ev["public_id"], "Date": ev["origin_time"][:10],
        "Band": s.get("chosen_band", "").replace("band_", ""),
        "Depth": p.get("depth_km", ""), "Mw": round(p["mw"], 2) if solved else "",
        "VR": round(p["vr"], 1) if solved else "",
        "Grade": s.get("quality", {}).get("grade", ""),
        "Selection": s.get("provenance", {}).get("selection_version", ""),
    }

    def _row(r: dict, used: bool) -> dict:
        snr = r.get("snr") or {}
        return {
            **common,
            "Station": r.get("station") if "." in str(r.get("station", ""))
            else f"{r.get('network', '')}.{r.get('station', '')}",
            "Distance_km": round(r["distance_km"], 1) if r.get("distance_km") else "",
            "Azimuth": round(r["azimuth"], 1) if r.get("azimuth") is not None else "",
            "Sector": r.get("sector", ""),
            "SNR_Z": snr.get("Z", ""), "SNR_R": snr.get("R", ""),
            "SNR_T": snr.get("T", ""), "SNR_med": r.get("snr_med", ""),
            "Amp_ratio": r.get("amp_ratio", ""),
            "Window_s": r.get("window_end_s", ""),
            "Used": used,
            "Reason_class": "" if used else invert.reason_class(r.get("reason", "")),
            "Reason": "" if used else r.get("reason", ""),
            "Station_VR": r.get("station_vr", ""),
            "Shift_s": r.get("zcor_s", ""),
        }
    return ([_row(r, True) for r in s.get("stations_used", [])]
            + [_row(r, False) for r in s.get("stations_dropped", [])
               if "network" in r or "." in str(r.get("station", ""))])


def build_catalogue(events_dir: Path | None = None) -> Path | None:
    """Scan <events_dir>/*/solution.json -> the three tables."""
    events_dir = events_dir or config.EVENTS_DIR
    solutions = sorted(events_dir.glob("*/solution.json"))
    if not solutions:
        return None
    published_ids: set[str] = set()
    if config.STATE_FILE.exists():
        state = json.loads(config.STATE_FILE.read_text())
        published_ids = {p["public_id"] for p in state.get("published", [])}

    rows, unpublished, ledger = [], [], []
    for path in solutions:
        s = json.loads(path.read_text())
        row = row_from(s, published_ids)
        rows.append(row)
        ledger.extend(ledger_rows(s))
        if not row["published"]:
            decision = s.get("publish_decision", {})
            unpublished.append({
                "PublicID": row["PublicID"], "Date": row["Date"],
                "GeoNet_M": row["GeoNet_M"], "Mw": row.get("Mw", ""),
                "Depth": row.get("Depth", ""), "Grade": row["Grade"],
                "VR": row.get("VR", ""), "DC": row.get("DC", ""),
                "NS": row["NS"], "PredDisp_cm": row.get("PredDisp_cm", ""),
                "Why_not": " | ".join(decision.get("reasons", [])) or
                           row["quality_flag"],
                "Tags": ";".join(decision.get("reason_tags", [])),
                "Selection": row["Selection"],
            })

    rows.sort(key=lambda r: r["Date"])
    out = events_dir / "catalogue.csv"
    _write(out, COLUMNS, rows)
    unpublished.sort(key=lambda r: r["Date"])
    _write(events_dir / "not_published.csv",
           ["PublicID", "Date", "GeoNet_M", "Mw", "Depth", "Grade", "VR",
            "DC", "NS", "PredDisp_cm", "Why_not", "Tags", "Selection"],
           unpublished)
    ledger.sort(key=lambda r: (r["Date"], r["PublicID"], r["Station"]))
    _write(events_dir / "station_ledger.csv", LEDGER_COLUMNS, ledger)
    print(f"catalogue: {len(rows)} events ({len(unpublished)} not published, "
          f"{len(ledger)} station rows) -> {events_dir}")
    return out


def _write(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in columns})


if __name__ == "__main__":
    build_catalogue()
