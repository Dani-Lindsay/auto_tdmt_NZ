"""Infrastructure configuration: paths, services, the precomputed Green's-
function grid, channel codes, and the names other modules import.

EVERY SCIENTIFIC CHOICE LIVES IN auto_tdmt.cfg (loaded by params.py) —
thresholds, distances, bands, the depth window, the grade table, the
publish gates. This module only re-exports those values under the names
the task modules use, so that the cfg is the single place to change a
number. If you are looking for a threshold, open auto_tdmt.cfg.
"""

import os
from pathlib import Path

import params

P = params.load()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Code repo (this directory) is ~/tools/internal/auto_tdmt_NZ.
# Local outputs live in ~/work/proj_tdmt_NZ; CI overrides OUTPUT_BASE to the
# repo checkout so results commit into events/.
REPO_DIR = Path(__file__).resolve().parent

OUTPUT_BASE = Path(
    os.environ.get("AUTO_TDMT_OUTPUT", Path.home() / "work" / "proj_tdmt_NZ")
)
EVENTS_DIR = Path(os.environ.get("AUTO_TDMT_EVENTS", OUTPUT_BASE / "outputs"))
GF_LIBRARY_DIR = Path(os.environ.get("AUTO_TDMT_GF", OUTPUT_BASE / "gf_library"))
# Raw-waveform download cache: miniSEED + StationXML per event, reused by
# every later run. DISPOSABLE — delete the whole directory at any time.
WF_CACHE_DIR = Path(os.environ.get("AUTO_TDMT_WFCACHE", OUTPUT_BASE / "wfcache"))
CPS_BIN = Path(
    os.environ.get("AUTO_TDMT_CPS_BIN", OUTPUT_BASE / "cps" / "PROGRAMS.330" / "bin")
)
STATE_FILE = REPO_DIR / "events" / "index.json"

# ---------------------------------------------------------------------------
# Code version stamped into every solution and catalogue row. Bump
# SELECTION_VERSION whenever selection or grading rules change.
#   selection-v3  the pre-funnel state       (git tag)
#   selection-v4  the four-pass funnel       (git tag)
#   v5            the published-rules version (Clinton loop, INGV grades)
# ---------------------------------------------------------------------------
SELECTION_VERSION = "v5"


def code_version() -> str:
    """Short git hash of HEAD (+ '-dirty'), or 'unknown' outside a checkout."""
    import subprocess
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_DIR,
            capture_output=True, text=True, timeout=10)
        if out.returncode != 0:
            return "unknown"
        commit = out.stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=REPO_DIR, capture_output=True, text=True, timeout=10)
        if dirty.returncode == 0 and dirty.stdout.strip():
            commit += "-dirty"
        return commit
    except Exception:  # noqa: BLE001 - provenance must never break a run
        return "unknown"


# Archive status: an event either has a solution or is recorded as having
# no coherent one. Every consumer branches on is_solved().
STATUS_SOLVED = "solved"
STATUS_NO_SOLUTION = "no_coherent_solution"


def is_solved(solution: dict) -> bool:
    return solution.get("status", STATUS_SOLVED) == STATUS_SOLVED


def slugify(text: str) -> str:
    """Locality -> filesystem-safe slug for event directory names."""
    import re
    return re.sub(r"[^A-Za-z0-9-]+", "-", text).strip("-")


def event_dir_name(public_id: str, mw: float, depth_km: float,
                   locality: str, origin_time: str = "") -> str:
    """<publicID>_<YYYY-MM-DD>_Mw<mw>_<depth>km_<locality>."""
    date = f"_{origin_time[:10]}" if origin_time else ""
    return (f"{public_id}{date}_Mw{mw:.1f}_{depth_km:g}km_"
            f"{slugify(locality)[:40]}")


def no_solution_dir_name(public_id: str, locality: str,
                         origin_time: str = "") -> str:
    date = f"_{origin_time[:10]}" if origin_time else ""
    return f"{public_id}{date}_NOSOL_{slugify(locality)[:40]}"


def find_event_dir(public_id: str, events_dir: Path | None = None):
    """Locate an event's directory whether plain or canonically named."""
    events_dir = events_dir or EVENTS_DIR
    matches = sorted(events_dir.glob(f"{public_id}*"))
    return matches[0] if matches else None

# ---------------------------------------------------------------------------
# GeoNet services
# ---------------------------------------------------------------------------
QUAKE_API = "https://api.geonet.org.nz"
FDSN_NRT = "https://service-nrt.geonet.org.nz"  # last ~8 days, <5 min latency
FDSN_ARCHIVE = "https://service.geonet.org.nz"  # complete, ~7 days behind
GEONET_CMT_CSV = (
    "https://raw.githubusercontent.com/GeoNet/data/main/"
    "moment-tensor/GeoNet_CMT_solutions.csv"
)
USER_AGENT = "auto_tdmt_NZ/0.1 (research MT pipeline; danielle.lindsay@earthsciences.nz)"
NRT_WINDOW_DAYS = 8  # older events must use the archive FDSN services

# ---------------------------------------------------------------------------
# Task 1 — processing floor (auto_tdmt.cfg section 1)
# ---------------------------------------------------------------------------
PROCESS_MIN_PRELIM_MAG = P.geonet.minPrelimMag
MAX_PROCESS_DEPTH_KM = P.geonet.maxDepthKm
PLACEHOLDER_DEPTHS_KM = set(P.geonet.placeholderDepthsKm)
PROCESS_EVENT_TYPES = {"earthquake"}
NZ_BBOX = dict(zip(("lat_min", "lat_max", "lon_min", "lon_max"),
                   P.geonet.bbox))

# ---------------------------------------------------------------------------
# Task 5 — publication gate (auto_tdmt.cfg section 5)
# ---------------------------------------------------------------------------
PUBLISH_MIN_GRADE = P.publish.minGrade
PUBLISH_MIN_MW = P.publish.minMw
PUBLISH_MIN_PRED_DISP_M = P.publish.minDisplacementM
MAX_POSTS_PER_DAY = P.publish.maxPerDay
AFTERSHOCK_RADIUS_KM = P.publish.aftershockRadiusKm
AFTERSHOCK_WINDOW_DAYS = P.publish.aftershockWindowDays
AFTERSHOCK_MW_MARGIN = P.publish.aftershockMwMargin

# ---------------------------------------------------------------------------
# Task 2 — stations and waveforms (auto_tdmt.cfg section 2)
# ---------------------------------------------------------------------------
NETWORK = "NZ"
CHANNEL_PRIORITY = ("HH?", "BH?")  # broadband only; short-period useless at LP
station_min_dist_km = P.min_dist_km
station_max_dist_km = P.max_dist_km
# one-shot radius extension when the usable pool is thin (offshore events)
MIN_USABLE_BEFORE_EXTEND = P.invert.minStations + 1
RADIUS_EXTEND_KM = 100.0
band_candidates = P.bands_hz
window_tail_s = P.window_tail_s
FILTER_CORNERS = P.station.filterCorners


def band_tag(band_hz: tuple[float, float]) -> str:
    return f"band_{round(1/band_hz[1]):d}-{round(1/band_hz[0]):d}s"


def sector(azimuth: float) -> int:
    """Azimuth -> sector index (P.station.sectors equal sectors)."""
    width = 360.0 / P.station.sectors
    return int(azimuth // width) % P.station.sectors

# Sample spacing / windows — locked between data prep and the GF library,
# following the mttime example notebooks (Chiang) exactly: data trimmed
# origin-30 s .. origin+200 s at dt=1 s; GFs computed with npts=256 (FK
# needs a power of 2), vred=0 and t0=0 so they start at origin; station-
# table ts=30 samples. The inversion window per station ends at
# station.maxWindowS after origin at most (INV_NPTS).
DT = 1.0  # s
GF_NPTS = 256
TIME_BEFORE_S = 30
TIME_AFTER_S = 200
INV_NPTS = int(TIME_BEFORE_S + P.station.maxWindowS)
assert INV_NPTS <= GF_NPTS, (
    f"station.maxWindowS {P.station.maxWindowS:g} s exceeds the Green's "
    f"functions ({GF_NPTS - TIME_BEFORE_S} s after origin)")
# Response-removal pre-filter (Hz), from mttime example notebook 01.
RESPONSE_PRE_FILT = (0.004, 0.007, 10.0, 20.0)

# ---------------------------------------------------------------------------
# Green's function library grid (precomputed; see auto_tdmt.cfg footer)
# ---------------------------------------------------------------------------
GF_DIST_KM = list(range(10, 505, 5))
# Fine near the surface where depth discrimination happens (0.5 km to 5 km,
# 1 km to 10 km), 2 km through the crust, 4 km below 30 km (Fiordland slab).
GF_DEPTHS_KM = (
    [i * 0.5 for i in range(2, 11)]      # 1.0-5.0 km @ 0.5
    + list(range(6, 11))                 # 6-10 km @ 1
    + list(range(12, 31, 2))             # 12-30 km @ 2
    + list(range(34, 61, 4))             # 34-58 km @ 4
)
# Velocity models (models/<name>.d, citation inside each file):
# Ristau (2008) SRL 79(3) Table 1 — the models GeoNet's own regional CMT
# analysis was built on, so our solutions are directly comparable.
GF_MODELS = ("nz_south_ristau2008", "nz_north_ristau2008")
GF_VERSION = "v1"


def model_for_event(latitude: float, longitude: float) -> str:
    """Crude North/South Island split (Cook Strait)."""
    if latitude <= -40.5 and (longitude < 175.0 or longitude > 200.0):
        return "nz_south_ristau2008"
    return "nz_north_ristau2008"

# ---------------------------------------------------------------------------
# Task 3 — inversion (auto_tdmt.cfg section 3)
# ---------------------------------------------------------------------------
INVERSION_DEGREE = P.invert.degree
MIN_STATIONS_USED = P.invert.minStations
# Surface-wave group velocity used to convert per-station zcor into a
# velocity-model deviation percentage, dV% = zcor / (dist / v_group) * 100
# (EPS207 velocity-model analysis; ~Love-wave group velocity).
GROUP_VELOCITY_KMS = 3.0

# ---------------------------------------------------------------------------
# Task 4 — Okada forward model (auto_tdmt.cfg section 4)
# ---------------------------------------------------------------------------
SHEAR_MODULUS_PA = P.forward.shearModulusPa
POISSON_NU = P.forward.poissonRatio
NISAR_REPEAT_DAYS = 12

# ---------------------------------------------------------------------------
# Publishing transports
# ---------------------------------------------------------------------------
BLUESKY_HANDLE_ENV = "BLUESKY_HANDLE"
BLUESKY_APP_PASSWORD_ENV = "BLUESKY_APP_PASSWORD"
AUTO_PUBLISH = os.environ.get("AUTO_PUBLISH", "false").lower() == "true"
