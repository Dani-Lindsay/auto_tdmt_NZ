"""NISAR acquisition timing at a point: last real acquisitions per
track/direction from NASA CMR, and the predicted next pass from the 12-day
repeat cycle.

Query pattern reused from NZ_observation_density/s1_nisar_point_timing.py
(real granules only — if CMR has nothing at the point, we say so; no dates
are ever fabricated).
"""

from __future__ import annotations

import re
from datetime import datetime, timedelta, timezone
from urllib.parse import urlencode

import requests

import config

# NISAR L2 GSLC collections on CMR (provisional + beta, as of 2026)
COLLECTIONS = ("NISAR_L2_GSLC_PROVISIONAL_V1", "NISAR_L2_GSLC_BETA_V1")
NI_RE = re.compile(r"NISAR_L2_PR_GSLC_(\d{3})_(\d{3})_([AD])_(\d{3})_")
LOOKBACK_DAYS = 60


def nisar_passes(lon: float, lat: float, when: datetime | None = None) -> list[dict]:
    """Last + next NISAR pass per (track, direction) covering the point.

    Returns a list of {track, direction, last_utc, next_utc} sorted by
    next_utc; empty list if CMR has no NISAR granules at the point.
    """
    now = when or datetime.now(timezone.utc)
    start = now - timedelta(days=LOOKBACK_DAYS)

    latest: dict[tuple[int, str], datetime] = {}
    for sn in COLLECTIONS:
        q = urlencode({
            "short_name": sn,
            "point": f"{lon:.4f},{lat:.4f}",
            "temporal": f"{start:%Y-%m-%dT%H:%M:%SZ},{now:%Y-%m-%dT%H:%M:%SZ}",
            "page_size": 200,
        })
        r = requests.get(
            f"https://cmr.earthdata.nasa.gov/search/granules.json?{q}",
            headers={"User-Agent": config.USER_AGENT}, timeout=90,
        )
        r.raise_for_status()
        for e in r.json()["feed"]["entry"]:
            name = e.get("producer_granule_id") or e["title"]
            m = NI_RE.match(name)
            if not m:
                continue
            track, direction = int(m.group(2)), m.group(3)
            t = datetime.fromisoformat(
                e["time_start"].replace("Z", "+00:00"))
            key = (track, direction)
            if key not in latest or t > latest[key]:
                latest[key] = t

    out = []
    for (track, direction), last in sorted(latest.items()):
        nxt = last
        while nxt <= now:
            nxt += timedelta(days=config.NISAR_REPEAT_DAYS)
        out.append({
            "track": track,
            "direction": direction,
            "last_utc": last.strftime("%Y-%m-%d %H:%MZ"),
            "next_utc": nxt.strftime("%Y-%m-%d %H:%MZ"),
        })
    out.sort(key=lambda d: d["next_utc"])
    return out
