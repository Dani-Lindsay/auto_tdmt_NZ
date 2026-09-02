"""Central configuration for the auto_tdmt_NZ pipeline.

Every threshold, filter band, grid definition and path lives here so nothing
is scattered or hardcoded in the processing modules. Values that are physics
(not preference) carry their derivation in a comment.
"""

import os
from pathlib import Path

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
CPS_BIN = Path(
    os.environ.get("AUTO_TDMT_CPS_BIN", OUTPUT_BASE / "cps" / "PROGRAMS.330" / "bin")
)
STATE_FILE = REPO_DIR / "events" / "index.json"

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

# NRT window: events older than this must use the archive FDSN services.
NRT_WINDOW_DAYS = 8

# ---------------------------------------------------------------------------
# Processing trigger (applied to GeoNet PRELIMINARY magnitude, which is a mix
# of M/MLv/mB summary magnitudes — deliberately loose; everything above this
# floor is processed and archived, publication is decided later on OUR Mw)
# ---------------------------------------------------------------------------
PROCESS_MIN_PRELIM_MAG = 4.0
# Restrict to events GeoNet can actually constrain (their network interest
# region); drop teleseisms tagged "outside of network interest".
PROCESS_EVENT_TYPES = {"earthquake"}
# Rough NZ bounding box guard (lat, lon); events outside are skipped.
NZ_BBOX = {"lat_min": -50.5, "lat_max": -33.0, "lon_min": 164.0, "lon_max": 182.5}

# ---------------------------------------------------------------------------
# Publication gate (applied to OUR inverted Mw after quality gates)
# ---------------------------------------------------------------------------
PUBLISH_MIN_MW = 5.0
PUBLISH_MIN_PRED_DISP_M = 0.01  # 1 cm predicted max surface displacement
MAX_POSTS_PER_DAY = 3
# Aftershock throttle: within this space/time window of an already-published
# event, a new event only publishes if within AFTERSHOCK_MW_MARGIN of the
# published mainshock Mw, or above PUBLISH_MIN_MW regardless.
AFTERSHOCK_RADIUS_KM = 75.0
AFTERSHOCK_WINDOW_DAYS = 14.0
AFTERSHOCK_MW_MARGIN = 0.5

# ---------------------------------------------------------------------------
# Waveform selection / pre-processing (EPS207 recipe; validated retrospectively)
# ---------------------------------------------------------------------------
NETWORK = "NZ"
CHANNEL_PRIORITY = ("HH?", "BH?")  # broadband only; short-period useless at LP
MAX_STATION_DIST_KM = 400.0
MIN_STATION_DIST_KM = 20.0
# Station distance must exceed ~3x source depth for the point-source /
# far-field assumptions used by TDMT (EPS207 §3.1).
MIN_DIST_DEPTH_RATIO = 3.0
MAX_STATIONS = 12  # nearest N passing QC; keeps inversion + plots tidy

# Zero-phase Butterworth passbands in Hz, ordered candidate lists per
# preliminary magnitude (BSL TDMT practice: a small menu of period bands,
# longer periods as magnitude grows). The pipeline tries each candidate and
# keeps the solution with the best variance reduction.
def band_candidates(prelim_mag: float) -> list[tuple[float, float]]:
    assert 0.0 < prelim_mag < 10.0, f"implausible magnitude {prelim_mag}"
    if prelim_mag < 4.5:
        return [(0.02, 0.10), (0.02, 0.05)]  # 10-50 s, 20-50 s
    if prelim_mag < 5.5:
        # 20-50 s, 10-50 s, 20-100 s
        return [(0.02, 0.05), (0.02, 0.10), (0.01, 0.05)]
    return [(0.01, 0.05), (0.01, 0.033)]  # 20-100 s, 30-100 s


def band_tag(band_hz: tuple[float, float]) -> str:
    return f"band_{round(1/band_hz[1]):d}-{round(1/band_hz[0]):d}s"


# Stations whose individual variance reduction falls below this after the
# first inversion pass are rejected and the inversion is rerun (EPS207 §3.1:
# drop persistently low-VR stations). With the SNR gate at Ristau's
# permissive 2.0, this is the filter that removes marginal stations that
# passed SNR but do not actually fit — keeping them dilutes %DC.
STATION_VR_FLOOR = 10.0

# Preferred-solution rule (EPS207 §3.3: VR alone is a weak depth
# discriminator; %DC is more diagnostic): among solutions whose VR is within
# this many percentage points of the maximum, prefer the highest %DC. Applied
# to the depth pick, the station-rejection reference and the band choice.
PREFER_DC_VR_TOLERANCE = 5.0

# Sample spacing / windows — locked between data prep and the GF library,
# following the mttime example notebooks (Chiang) exactly:
# data trimmed origin-30 s .. origin+200 s at dt=1 s; GFs computed with
# npts=256 (FK needs a power of 2), vred=0 and t0=0 so they start at origin;
# station-table ts=30 samples, inversion window npts=150.
DT = 1.0  # s
GF_NPTS = 256
INV_NPTS = 150
TIME_BEFORE_S = 30
TIME_AFTER_S = 200
FILTER_CORNERS = 3  # obspy corners, zerophase=True; notebook-02 values, applied
# identically to data and Green's functions (only consistency matters)
# Response-removal pre-filter (Hz), from mttime example notebook 01.
RESPONSE_PRE_FILT = (0.004, 0.007, 10.0, 20.0)

# SNR gate: min over components of RMS(signal, origin..+200 s) /
# RMS(noise, -120..-10 s), measured in the inversion passband. Ristau (2008):
# "a SNR higher than 2 is normally required to calculate a reliable moment
# tensor" — stations below this are dropped (loudly, recorded in provenance).
MIN_SNR = 2.0
# Tiered acceptance: when fewer than TIER_TARGET_STATIONS pass MIN_SNR,
# stations down to SNR_TIER_LOW are admitted (tagged snr_tier="low") so
# sparse/coda-contaminated events keep azimuthal coverage; the letter
# quality grade, not a hard gate, then tells the reader what it is worth.
SNR_TIER_LOW = 1.2
TIER_TARGET_STATIONS = 5

# ---------------------------------------------------------------------------
# Green's function library grid
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


def model_for_event(latitude: float, longitude: float) -> str:
    """Crude North/South Island split (Cook Strait); refine if offshore
    regions ever need their own model."""
    if latitude <= -40.5 and (longitude < 175.0 or longitude > 200.0):
        return "nz_south_ristau2008"
    return "nz_north_ristau2008"
GF_VERSION = "v1"

# ---------------------------------------------------------------------------
# Inversion / quality gates
# ---------------------------------------------------------------------------
INVERSION_DEGREE = 5  # deviatoric
MIN_STATIONS_USED = 3
MIN_VR_PUBLISH = 50.0  # % variance reduction floor for publication
MAX_AZ_GAP_DEG = 270.0  # flag (not fail) beyond this
# Depth search: full library depth range when GeoNet depth is a placeholder
# (their fixed depths are unreliable); otherwise +/- this margin.
DEPTH_SEARCH_MARGIN_KM = 20.0
PLACEHOLDER_DEPTHS_KM = {5.0, 12.0, 33.0}  # GeoNet fixed-depth values

# Surface-wave group velocity used to convert per-station zcor into a
# velocity-model deviation percentage, dV% = zcor / (dist / v_group) * 100
# (EPS207 velocity-model analysis; ~Love-wave group velocity).
GROUP_VELOCITY_KMS = 3.0

# ---------------------------------------------------------------------------
# Okada forward model
# ---------------------------------------------------------------------------
SHEAR_MODULUS_PA = 3.0e10
POISSON_NU = 0.25
FORWARD_GRID_HALFWIDTH_KM = 60.0
FORWARD_GRID_STEP_KM = 1.0
# NISAR L-band LOS geometry defaults (right-looking, ~34-48 deg incidence);
# per-track values refined in nisar_dates.py when a real granule is found.
NISAR_REPEAT_DAYS = 12

# ---------------------------------------------------------------------------
# Publishing
# ---------------------------------------------------------------------------
BLUESKY_HANDLE_ENV = "BLUESKY_HANDLE"
BLUESKY_APP_PASSWORD_ENV = "BLUESKY_APP_PASSWORD"
AUTO_PUBLISH = os.environ.get("AUTO_PUBLISH", "false").lower() == "true"
