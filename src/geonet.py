"""GeoNet access: quake API (detection), FDSN clients (waveforms/metadata),
and the manually-produced GeoNet CMT catalogue (validation).

All requests carry a descriptive User-Agent and fail loudly — no silent
retries into stale data.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, asdict

import pandas as pd
import requests
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

import config


@dataclass
class Event:
    public_id: str
    origin_time: str  # ISO8601
    longitude: float
    latitude: float
    depth_km: float
    prelim_mag: float
    mag_type: str
    locality: str
    quality: str
    # from GeoNet QuakeML when available (None otherwise)
    depth_unc_km: float | None = None
    eval_mode: str | None = None
    eval_status: str | None = None
    depth_type: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def _get(url: str, **kwargs) -> requests.Response:
    headers = kwargs.pop("headers", {})
    headers["User-Agent"] = config.USER_AGENT
    headers.setdefault("Accept-Encoding", "gzip")
    r = requests.get(url, headers=headers, timeout=60, **kwargs)
    r.raise_for_status()
    return r


def recent_quakes(mmi: int = 3) -> list[Event]:
    """Recent events from the real-time quake API (GeoJSON, last 365 d,
    max 100). MMI is the only server-side filter; magnitude filtering is
    done by the caller."""
    r = _get(
        f"{config.QUAKE_API}/quake",
        params={"MMI": mmi},
        headers={"Accept": "application/vnd.geo+json;version=2"},
    )
    events = []
    for f in r.json()["features"]:
        p = f["properties"]
        lon, lat = f["geometry"]["coordinates"][:2]
        events.append(
            Event(
                public_id=p["publicID"],
                origin_time=p["time"],
                longitude=float(lon),
                latitude=float(lat),
                depth_km=float(p["depth"]),
                prelim_mag=float(p["magnitude"]),
                mag_type="M",  # quake API reports the summary magnitude only
                locality=p["locality"],
                quality=p["quality"],
            )
        )
    return events


def _quakeml_origin_details(public_id: str) -> dict:
    """Depth uncertainty and review status from GeoNet QuakeML (FDSN event
    service). Best-effort: {} on any failure — never blocks processing."""
    try:
        from obspy import read_events

        cat = read_events(
            f"{config.FDSN_ARCHIVE}/fdsnws/event/1/query"
            f"?eventid={public_id}")
        o = cat[0].preferred_origin() or cat[0].origins[0]
        unc = getattr(o.depth_errors, "uncertainty", None) \
            if o.depth_errors else None
        return {
            "depth_unc_km": round(unc / 1000.0, 1) if unc else None,
            "eval_mode": str(o.evaluation_mode) if o.evaluation_mode else None,
            "eval_status": (str(o.evaluation_status)
                            if o.evaluation_status else None),
            "depth_type": str(o.depth_type) if o.depth_type else None,
        }
    except Exception:  # noqa: BLE001
        return {}


def get_event(public_id: str) -> Event:
    """Single event by publicID from the quake API."""
    r = _get(
        f"{config.QUAKE_API}/quake/{public_id}",
        headers={"Accept": "application/vnd.geo+json;version=2"},
    )
    feats = r.json()["features"]
    assert len(feats) == 1, f"{public_id}: expected 1 feature, got {len(feats)}"
    p = feats[0]["properties"]
    lon, lat = feats[0]["geometry"]["coordinates"][:2]
    return Event(
        public_id=p["publicID"],
        origin_time=p["time"],
        longitude=float(lon),
        latitude=float(lat),
        depth_km=float(p["depth"]),
        prelim_mag=float(p["magnitude"]),
        mag_type="M",
        locality=p["locality"],
        quality=p["quality"],
        **_quakeml_origin_details(public_id),
    )


def is_nrt(origin_time: UTCDateTime) -> bool:
    """True if the event is young enough for the near-real-time FDSN
    services (which only hold the last ~8 days)."""
    age_days = (UTCDateTime.now() - origin_time) / 86400.0
    return age_days < config.NRT_WINDOW_DAYS


def fdsn_client(origin_time: UTCDateTime) -> Client:
    """NRT client for fresh events, archive client otherwise. The archive
    lags ~7 days behind real time, so for a young event NRT is the ONLY
    source — never fall back silently.

    GeoNet's NRT server transiently serves an incomplete service document
    ("client does not have a dataselect service"), so retry construction a
    few times before failing loudly."""
    import time

    base = config.FDSN_NRT if is_nrt(origin_time) else config.FDSN_ARCHIVE
    last_err = None
    for attempt in range(3):
        client = Client(base_url=base)
        if "dataselect" in client.services and "station" in client.services:
            return client
        last_err = f"incomplete service document from {base}"
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"FDSN client unusable after retries: {last_err}")


def load_geonet_cmt() -> pd.DataFrame:
    """GeoNet's manual CMT catalogue (validation ground truth).

    NOTE: MT elements in the CSV are in units of 1e20 dyne-cm.
    """
    r = _get(config.GEONET_CMT_CSV)
    df = pd.read_csv(io.StringIO(r.text))
    assert "PublicID" in df.columns or "EVENT_ID" in df.columns, (
        f"unexpected CMT CSV columns: {list(df.columns)[:8]}"
    )
    return df
