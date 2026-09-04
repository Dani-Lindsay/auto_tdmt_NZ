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

# ---------------------------------------------------------------------------
# Code version stamped into every solution and catalogue row, so a catalogue
# with mixed-vintage rows is self-describing and any solution can be traced
# to the exact code that produced it (git tag selection-v3 = the pre-funnel
# state). Bump SELECTION_VERSION whenever selection or grading rules change.
# ---------------------------------------------------------------------------
SELECTION_VERSION = "v4"


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
# no coherent one (see invert.no_solution_record). Every consumer branches
# on is_solved() rather than assuming "preferred" exists.
STATUS_SOLVED = "solved"
STATUS_NO_SOLUTION = "no_coherent_solution"


def is_solved(solution: dict) -> bool:
    return solution.get("status", STATUS_SOLVED) == STATUS_SOLVED


def slugify(text: str) -> str:
    """Locality -> filesystem-safe slug for event directory names."""
    import re
    return re.sub(r"[^A-Za-z0-9-]+", "-", text).strip("-")


def event_dir_name(public_id: str, mw: float, depth_km: float,
                   locality: str) -> str:
    return (f"{public_id}_Mw{mw:.1f}_{depth_km:g}km_"
            f"{slugify(locality)[:40]}")


def no_solution_dir_name(public_id: str, locality: str) -> str:
    """Directory name for an event with no coherent solution — visibly
    different from a solved event's Mw/depth name."""
    return f"{public_id}_NOSOL_{slugify(locality)[:40]}"


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
# Lowered 4.0 -> 3.7 (2026-09-04) after the small-event review: GeoNet
# preliminary M runs ~0.2-0.3 above the final Mw for small events, so a
# 3.7 prelim floor reaches true Mw ~3.4-3.5 - the smallest this network
# + 10-50 s band demonstrably constrains (2026p216882 Mw 3.6 grade C,
# 2026p221690 Mw 3.77 grade B). A 3.5 prelim floor would mostly add
# unsolvable Mw~3.2 events (+30% compute) below Ristau's own coverage.
PROCESS_MIN_PRELIM_MAG = 3.7
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


# When fewer than this many stations enter the pool, the search radius is
# extended once by RADIUS_EXTEND_KM (offshore events: 2026p047833 review)
# and the annulus is fetched and processed too.
MIN_USABLE_BEFORE_EXTEND = 4
RADIUS_EXTEND_KM = 100.0

# Station clustering: dense sub-networks (e.g. the Ruapehu volcano ring)
# would let one site dominate an azimuth sector. Beyond CLUSTER_MAX_STATIONS
# within CLUSTER_RADIUS_KM, a station is TAGGED "cluster_surplus" — a
# DEMOTION, NEVER AN EXCLUSION (v4): it still gets its seat-earning test in
# the funnel, but must not reduce the joint fit to be admitted.
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
# The screen needs a network median to compare against; 3 is the minimum
# that gives one (was 4, which left the sparsest events unprotected —
# 2026p101368 used RDHZ at station VR -40 because the screen never ran).
AMP_SCREEN_MIN_STATIONS = 3
# Station distance below ~3x source depth violates the point-source /
# far-field heuristic (EPS207 §3.1) — but the CPS Green's functions are
# complete-wavefield solutions including the near-field terms, and the
# 2026-09-04 audit found this rule removing the CLOSEST station in every
# event it touched (94% of those ended C/D). It is now a DEMOTION TAG
# ("near_field"), NEVER AN EXCLUSION: the fit decides.
MIN_DIST_DEPTH_RATIO = 3.0
MAX_POOL_STATIONS = 30  # pool safety cap; the funnel prunes by fit

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

# (Removed 2026-09-04: the greedy "test-drop what improves joint VR" pass.
# It was evicting well-fitting stations to polish a number — 2026p091845
# lost THZ/MRZ/WRRZ at station VR 54-65 — which is VR vanity at the cost of
# azimuth coverage. Stations now leave only if they anti-fit.)

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

# ---------------------------------------------------------------------------
# Station selection v4 — the funnel (2026-09-04)
#
# Design (D. Lindsay): use ALL stations, drop the worst-VR ones keeping the
# majority (~10), re-run, choose the best 4-6 as the core, then incrementally
# add every other station back and decide which are useful, considering
# azimuth spread and individual fit. Nothing usable is deleted by a
# pre-filter; the fit decides. Motivation: the 693-event audit found 61% of
# all station exclusions were decided by the solution itself, and 69% of
# those vetoes were issued by cores whose own VR was below 20.
# ---------------------------------------------------------------------------

# Peak-to-noise station quality: median over Z/R/T of
# peak|signal| / RMS(pre-event noise), signal measured ONLY inside the
# distance-adaptive window actually inverted — an impulsive surface-wave
# packet is a spike above background, which RMS-over-200s could not see.
# Below PEAK_NOISE_DEAD the channel carries nothing and is hard-rejected
# (data still written for the figures); between DEAD and STRONG the station
# is TAGGED "weak_signal" and must earn its seat in the funnel.
PEAK_NOISE_DEAD = 1.2
PEAK_NOISE_STRONG = 5.0

# Time-shift (zcor) sanity. mttime's cross-correlation search is unbounded
# (a 60 s window can slide >100 s), which is exactly how a noise trace finds
# a chance alignment and earns undeserved VR. Archive calibration: stations
# used in grade A/B solutions have |zcor| <= 9 s at p95, while grade-D
# stations reach 0.88 of their whole travel time. Pass 1 therefore runs with
# NO shifts at all, and from pass 2 on a station whose solved shift exceeds
# this bound is rejected (its fit is not evidence).
ZCOR_MAX_S = 8.0

# Pass 1: survey inversion over the whole pool, coarse depth grid (every
# PASS1_DEPTH_STRIDE-th library depth plus both ends). Stations are ranked by
# the MEDIAN of their own VR over the contiguous VR plateau — chance
# alignment is depth-specific, real coherence is not.
PASS1_DEPTH_STRIDE = 2
PASS1_KEEP_N = 10          # "keep the majority"
PASS1_KEEP_MAX = 12        # after filling empty azimuth sectors
PASS1_KEEP_VR_MIN = 10.0   # a station must fit at least this well to survive

# Pass 2 -> core: the best few by own VR, big enough to constrain a
# mechanism, small enough to stay clean. The core must contain two stations
# at least AZ_PAIR_MIN_DEG apart or it cannot resolve a mechanism at all.
CORE_SIZE_MIN = 3
CORE_SIZE_MAX = 6
CORE_VR_MIN = 30.0

# Pass 3 — earn your seat: every non-core station (tagged ones included) is
# added at the core depth and kept if the core solution predicts its
# waveform. Sparse cores relax the floor: with 3-4 stations, extra azimuth
# coverage is worth a weaker individual fit (2026p189537 review).
ADMIT_VR_MIN = 30.0
ADMIT_VR_MIN_SPARSE = 20.0
SPARSE_CORE_COUNT = 5
# A station that fills an EMPTY azimuth sector is worth admitting on a
# weaker fit — geometry is the scarcer commodity in NZ.
ADMIT_SECTOR_VR_MIN = 10.0
# ...but no addition may cost more than this many joint VR points unless it
# brings a new sector.
ADMIT_MAX_JOINT_VR_DROP = 3.0
# Upper bound on the used set (2026p033598 review: an over-full network
# diluted %DC from 98 to 46).
MAX_USED_STATIONS = 12

# A station whose own VR is negative fits worse than silence: it only steers
# the tensor, so it is culled unconditionally (no joint-gain test, no
# sector protection — 2026p238013 NNZ sat at -41 behind a passing joint VR).
ANTIFIT_VR = 0.0

# Below this joint VR nothing coheres: the event is archived as
# "no coherent solution" rather than publishing a junk mechanism.
NO_SOLUTION_VR = 20.0

# Depth-pick guard: a VR maximum sitting on the first/last grid depth with a
# zero-width plateau is an artifact (the smoothest GFs absorb noise). Prefer
# the best INTERIOR local maximum within this tolerance — 2026p508890 rode a
# 58 km edge (VR 22.7) over the physical 8 km peak (VR 18.3, DC 88-96).
GRID_EDGE_VR_TOLERANCE = 5.0

# Mechanism geometry sanity: two used stations at least this far apart in
# azimuth. Every grade A/B solution in the archive already satisfies it;
# 330 of 559 grade-D solutions do not.
AZ_PAIR_MIN_DEG = 90.0

# Band escalation: if the inverted Mw overshoots the preliminary magnitude
# by this much, the event is bigger than the band menu assumed — run the
# next band up and prefer it if it fits at least as well (2026p336046:
# prelim 4.2 -> Mw 5.01 out of a 10-50 s-only run).
BAND_ESCALATE_DMW = 0.6

# ---------------------------------------------------------------------------
# Quality grades v2 (2026-09-04) — evidence, not station counts
#
# Station count and azimuthal gap are GONE as thresholds (D. Lindsay: "so
# long as you have stations 90 degrees from each other you can make a good
# solution"; 3 well-fitting stations can be grade A, as at BSL). What
# remains is measured evidence: how well the solution explains the data
# (VR), whether every used station actually fits (no passengers), whether
# the mechanism survives leaving a station out (jackknife rotation), and
# whether the depth is a real interior plateau rather than a grid-edge
# artifact. B additionally requires DC >= 60, which makes B exactly the
# BSL publishability rule (VR >= 60 and DC >= 60).
#
# Thresholds are the archive's own A/B statistics: min own-station VR p25
# 41 (A) / 35 (B); jackknife max rotation p90 12 deg (A), p75 14.5 (B).
# DC never enters SELECTION — a noise station can inflate it.
# ---------------------------------------------------------------------------
GRADE_RUBRIC = {
    "A": dict(vr=70.0, dc=60.0, min_own_vr=40.0, jk_rot_max=15.0,
              need_jackknife=True),
    "B": dict(vr=60.0, dc=60.0, min_own_vr=25.0, jk_rot_max=25.0,
              need_jackknife=False),
    "C": dict(vr=50.0, min_own_vr=10.0),
}


def sector(azimuth: float) -> int:
    """Azimuth -> one of 8 x 45 degree sectors."""
    return int(azimuth // 45) % 8

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
# The depth search always covers the full GF grid (INDEPENDENT BY CHOICE,
# 2026-09-03): bounding it around GeoNet's depth would force agreement and
# destroy the catalogue's value as an independent check.
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
