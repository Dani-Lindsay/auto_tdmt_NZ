"""Processing floor and publication gate.

Processing floor (cheap, on GeoNet's preliminary magnitude — a mixed bag of
M/MLv/mB types, so deliberately loose): decides what gets computed.

Publication gate (on OUR inverted Mw + the Okada forward model): decides
what gets emailed. Everything processed is archived either way.
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone

from obspy.geodetics.base import gps2dist_azimuth

import config
from geonet import Event


def passes_processing_floor(event: Event) -> tuple[bool, str]:
    if event.quality == "deleted":
        return False, "deleted by GeoNet"
    if event.prelim_mag < config.PROCESS_MIN_PRELIM_MAG:
        return False, (
            f"prelim M{event.prelim_mag:.1f} < {config.PROCESS_MIN_PRELIM_MAG}"
        )
    if (event.depth_km > config.MAX_PROCESS_DEPTH_KM
            and event.depth_km not in config.PLACEHOLDER_DEPTHS_KM):
        return False, (
            f"depth {event.depth_km:g} km > {config.MAX_PROCESS_DEPTH_KM:g} "
            f"km (below GF library; no surface displacement possible)"
        )
    b = config.NZ_BBOX
    if not (b["lat_min"] <= event.latitude <= b["lat_max"]
            and b["lon_min"] <= event.longitude % 360.0 <= b["lon_max"]):
        return False, "outside NZ bounding box"
    return True, "ok"


def publish_decision(solution: dict, forward: dict,
                     published_history: list[dict]) -> dict:
    """Gate a finished solution for publication.

    published_history: list of previously PUBLISHED events, each
    {public_id, mw, latitude, longitude, published_utc}.
    """
    pref = solution["preferred"]
    ev = solution["event"]
    reasons = []

    grade = solution["quality"].get("grade", "?")
    if not solution["quality"]["passed"]:
        reasons.append(
            f"quality grade {grade} (email requires A or B): "
            f"{solution['quality']['checks']}")

    mw_ok = pref["mw"] >= config.PUBLISH_MIN_MW
    disp_ok = forward["peak_abs_m"] >= config.PUBLISH_MIN_PRED_DISP_M
    if not (mw_ok or disp_ok):
        reasons.append(
            f"Mw {pref['mw']:.2f} < {config.PUBLISH_MIN_MW} and predicted "
            f"displacement {forward['peak_abs_m'] * 100:.2f} cm < "
            f"{config.PUBLISH_MIN_PRED_DISP_M * 100:g} cm"
        )

    now = datetime.now(timezone.utc)
    today = [
        p for p in published_history
        if now - datetime.fromisoformat(p["published_utc"]) < timedelta(days=1)
    ]
    if len(today) >= config.MAX_POSTS_PER_DAY:
        reasons.append(f"daily cap {config.MAX_POSTS_PER_DAY} reached")

    # aftershock throttle: nearby recent published mainshock dominates
    for p in published_history:
        age = now - datetime.fromisoformat(p["published_utc"])
        if age > timedelta(days=config.AFTERSHOCK_WINDOW_DAYS):
            continue
        dist_km = gps2dist_azimuth(
            ev["latitude"], ev["longitude"], p["latitude"], p["longitude"]
        )[0] / 1000.0
        if (dist_km <= config.AFTERSHOCK_RADIUS_KM
                and pref["mw"] < p["mw"] - config.AFTERSHOCK_MW_MARGIN
                and pref["mw"] < config.PUBLISH_MIN_MW):
            reasons.append(
                f"aftershock throttle: Mw {pref['mw']:.2f} within "
                f"{dist_km:.0f} km of published {p['public_id']} "
                f"(Mw {p['mw']:.2f})"
            )
            break

    return {
        "publish": not reasons,
        "reasons": reasons or ["passed all publication gates"],
        "mw_gate": mw_ok,
        "displacement_gate": disp_ok,
        "predicted_peak_disp_cm": forward["peak_abs_m"] * 100.0,
    }
