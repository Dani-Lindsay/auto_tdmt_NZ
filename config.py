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
# Raw-waveform download cache (laptop iteration): miniSEED + StationXML per
# event, written on first download and reused by every later sweep, so
# repeated rule-testing does not re-download from GeoNet. DISPOSABLE —
# delete the whole directory at any time (rm -rf), nothing else references
# it; a README.txt inside says the same.
WF_CACHE_DIR = Path(os.environ.get("AUTO_TDMT_WFCACHE", OUTPUT_BASE / "wfcache"))
CPS_BIN = Path(
    os.environ.get("AUTO_TDMT_CPS_BIN", OUTPUT_BASE / "cps" / "PROGRAMS.330" / "bin")
)
STATE_FILE = REPO_DIR / "events" / "index.json"


def slugify(text: str) -> str:
    """Locality -> filesystem-safe slug for event directory names."""
    import re
    return re.sub(r"[^A-Za-z0-9-]+", "-", text).strip("-")


def event_dir_name(public_id: str, mw: float, depth_km: float,
                   locality: str) -> str:
    return (f"{public_id}_Mw{mw:.1f}_{depth_km:g}km_"
            f"{slugify(locality)[:40]}")


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

# NRT window: events older than this must use the archive FDSN services.
NRT_WINDOW_DAYS = 8

# ---------------------------------------------------------------------------
# Processing trigger (applied to GeoNet PRELIMINARY magnitude, which is a mix
# of M/MLv/mB summary magnitudes — deliberately loose; everything above this
# floor is processed and archived, publication is decided later on OUR Mw)
# ---------------------------------------------------------------------------
PROCESS_MIN_PRELIM_MAG = 4.0
# Events deeper than this cannot produce measurable surface displacement
# at these magnitudes (the purpose of this tool). GeoNet PLACEHOLDER depths
# (incl. 33 km) are exempt: their true depth is unknown and often shallow,
# so they are processed and the depth search decides.
MAX_PROCESS_DEPTH_KM = 30.0
# The distance > 3x depth far-field guard is a SHALLOW-source rule; beyond
# this depth every station is far-field and the rule is skipped.
DIST_DEPTH_RULE_MAX_DEPTH_KM = 40.0
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
MAX_STATION_DIST_KM = 300.0  # base ceiling; see station_max_dist_km
MIN_STATION_DIST_KM = 20.0  # near-field exclusion radius (larger events)


def station_min_dist_km(prelim_mag: float) -> float:
    """Near-field exclusion. Small shallow events put their information in
    the close stations (2026-09-03 review: 'we are missing out on info');
    an M4 rupture is ~1 km so the point-source assumption already holds at
    10 km. Larger events keep the 20 km exclusion."""
    return 10.0 if prelim_mag < 4.5 else MIN_STATION_DIST_KM


# When fewer than this many stations survive the peak/noise floor, the
# search radius is extended once by RADIUS_EXTEND_KM (offshore events:
# 2026p047833 review) and the annulus is fetched and processed too.
MIN_USABLE_BEFORE_EXTEND = 4
RADIUS_EXTEND_KM = 100.0

# Station clustering: dense sub-networks (e.g. the Ruapehu volcano ring)
# would let one site dominate the azimuth sector. Keep at most
# CLUSTER_MAX_STATIONS (the best by peak/noise) within CLUSTER_RADIUS_KM
# of each other; the rest are dropped with the cluster named in the reason.
CLUSTER_RADIUS_KM = 25.0
CLUSTER_MAX_STATIONS = 2


def station_max_dist_km(prelim_mag: float) -> float:
    """Magnitude-scaled search radius: small events attenuate below
    usefulness at far field, larger events still carry information there."""
    if prelim_mag < 4.0:
        return 120.0
    if prelim_mag < 4.5:
        return 180.0
    if prelim_mag < 5.0:
        return 250.0
    return MAX_STATION_DIST_KM


# Amplitude-consistency screen: a station whose distance-corrected peak
# amplitude (peak x distance) is more than this factor from the network
# median has broken response metadata or severe site pathology and would
# steer the least-squares moment — drop it before inversion.
AMPLITUDE_OUTLIER_FACTOR = 8.0
# Station distance must exceed ~3x source depth for the point-source /
# far-field assumptions used by TDMT (EPS207 §3.1).
MIN_DIST_DEPTH_RATIO = 3.0
MAX_STATIONS = 30  # pool safety cap; backward elimination prunes

# Zero-phase Butterworth passbands in Hz, ordered candidate lists per
# preliminary magnitude (BSL TDMT practice: a small menu of period bands,
# longer periods as magnitude grows). The pipeline tries each candidate and
# keeps the solution with the best variance reduction.
def band_candidates(prelim_mag: float) -> list[tuple[float, float]]:
    assert 0.0 < prelim_mag < 10.0, f"implausible magnitude {prelim_mag}"
    if prelim_mag < 4.5:
        # 10-50 s only: small events carry no coherent energy above ~20 s
        # period, and every reviewed case (2026-09-03) preferred 10-50 s;
        # the 20-50 s trials only ever fit noise (Mw inflation). Other
        # bands remain testable via run02 --band.
        return [(0.02, 0.10)]
    if prelim_mag < 5.5:
        # ORDERED PREFERENCE (first gate-passer wins below M5.5 — see
        # run02): 10-50 s first — the band where signal lives at these
        # magnitudes; 20-50 s is the fallback. 20-100 s pruned 2026-09-03
        # (won 1/19 in this bin and only ever fit noise on the losers).
        return [(0.02, 0.10), (0.02, 0.05)]
    return [(0.01, 0.05), (0.01, 0.033)]  # 20-100 s, 30-100 s


def band_tag(band_hz: tuple[float, float]) -> str:
    return f"band_{round(1/band_hz[1]):d}-{round(1/band_hz[0]):d}s"

# Greedy improvement elimination: a station above the floor is still
# dropped if removing it improves the joint fit by at least this many VR
# points (the manual "test around and drop what degrades" practice).
ELIMINATION_VR_GAIN = 2.0

# The %DC tie-break only engages when the fit is meaningful; below this
# VR maximum the event is junk-grade and DC differences are noise — take
# the plain VR maximum (2026p348732: VR max 7.7 made the whole grid a
# "plateau" and DC alone chose the depth).
DC_TIEBREAK_MIN_VR = 20.0

# Distance-adaptive inversion window (record length), ALL magnitudes: a
# close station's surface-wave train is over quickly regardless of event
# size, and fitting the empty tail only taxes VR. Per-station window =
# 30 s pre-origin + dist/group_vel + tail(M), clamped; the tail grows
# with magnitude because larger sources ring longer.
# 2.5 km/s (was 2.8): 2026p091845 review showed the tail of the
# surface-wave train clipped at the old velocity; the slower bound
# plus the longer small-event tail keeps the full packet in-window
WINDOW_GROUP_VEL_KMS = 2.5
WINDOW_MIN_S = 60.0


def window_tail_s(prelim_mag: float) -> float:
    if prelim_mag < 4.5:
        return 30.0
    if prelim_mag < 5.5:
        return 40.0
    return 60.0

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

# Peak-to-noise station quality (replaces the RMS SNR gate 2026-09-03;
# calibrated against manual keep/toss labelling of 2026p283255 + review of
# five further events). Metric: median over Z/R/T of
# peak|signal| / RMS(pre-event noise), signal measured ONLY inside the
# distance-adaptive window actually inverted — an impulsive surface-wave
# packet is a spike above background, which RMS-over-200s could not see.
# >= PEAK_NOISE_CORE: trusted core, inverted from the start.
# >= PEAK_NOISE_FLOOR: candidate ("yellow"): admitted only if the core
#    solution predicts its waveform (station VR >= CANDIDATE_STATION_VR_MIN
#    when added alone at the core's preferred depth).
# below the floor: dead channel, rejected outright (data kept for figures).
PEAK_NOISE_CORE = 5.0
PEAK_NOISE_FLOOR = 2.0
CANDIDATE_STATION_VR_MIN = 30.0
# sparse events (2026p189537 review: 3 stations is thin): when the
# core holds fewer than 5 stations, admission relaxes to this floor
# so azimuth constraint is not starved by one strict threshold
CANDIDATE_STATION_VR_MIN_SPARSE = 20.0
SPARSE_CORE_COUNT = 5
# azimuth-coverage cap (2026p033598 review: too many stations diluted
# DC): when more survive, keep the best-VR station per 45 deg sector
# plus the top-4 VR overall — coverage first, then fit
MAX_USED_STATIONS = 12

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
# (their fixed depths are unreliable). For real located depths the search
# is bounded by GeoNet's own depth uncertainty when QuakeML provides it —
# margin = clamp(2 x depthUncertainty, MIN, MAX) — else the default: a
# located GeoNet depth is typically good to ~+/-5 km, and an unconstrained
# search was occasionally preferring depths wildly inconsistent with it.
DEPTH_SEARCH_MARGIN_KM = 10.0
DEPTH_SEARCH_MARGIN_MIN_KM = 5.0
DEPTH_SEARCH_MARGIN_MAX_KM = 15.0
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
