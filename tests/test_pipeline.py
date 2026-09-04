"""Anchor tests: fixed expectations that catch silent breakage of the
pipeline's load-bearing conventions."""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config  # noqa: E402
import greens  # noqa: E402
import invert  # noqa: E402
import okada_forward  # noqa: E402
import trigger  # noqa: E402
from geonet import Event  # noqa: E402


# --- config / selection -----------------------------------------------------

def test_band_candidates_ordering():
    small = config.band_candidates(4.0)
    assert small[0] == (0.02, 0.10), "small events lead with the 10-50 s band"
    large = config.band_candidates(6.0)
    assert all(hi <= 0.05 for _, hi in large), "large events exclude 10 s energy"


def test_band_tag():
    assert config.band_tag((0.01, 0.05)) == "band_20-100s"
    assert config.band_tag((0.02, 0.05)) == "band_20-50s"


def test_model_for_event():
    assert config.model_for_event(-44.45, 168.3) == "nz_south_ristau2008"  # Fiordland
    assert config.model_for_event(-36.8, 174.8) == "nz_north_ristau2008"  # Auckland
    assert config.model_for_event(-42.4, 173.7) == "nz_south_ristau2008"  # Kaikoura


def test_pick_preferred_prefers_dc_within_tolerance():
    # VRs within 5 points: pick max DC; far-below VR never wins on DC
    rows = [(80.0, 20.0), (78.0, 60.0), (60.0, 99.0)]
    assert invert.pick_preferred(rows) == 1
    assert invert.pick_preferred([(80.0, 20.0), (70.0, 90.0)]) == 0


def test_pick_preferred_junk_fit_skips_dc_tiebreak():
    # VR max below 20: everything is noise-grade, take plain VR max
    rows = [(6.8, 73.0), (3.5, 93.3), (7.7, 54.7)]
    assert invert.pick_preferred(rows) == 2


def test_pick_preferred_bimodal_vr_stays_contiguous():
    # deep lobe grazes the tolerance window but is disconnected from the
    # VR maximum: it must NOT steal the pick on DC (2026p660272 case)
    rows = [(91.0, 81.0), (90.0, 55.0), (88.0, 60.0), (30.0, 90.0),
            (87.0, 97.0)]
    assert invert.pick_preferred(rows) == 0
    # band selection (unordered) still considers all within tolerance
    assert invert.pick_preferred(rows, contiguous=False) == 4


# --- Green's functions ------------------------------------------------------

def test_nearest_grid_distance():
    assert greens.nearest_grid_distance(101.9) == 100
    assert greens.nearest_grid_distance(102.6) == 105
    with pytest.raises(AssertionError):
        greens.nearest_grid_distance(700.0)


@pytest.mark.parametrize("model", config.GF_MODELS)
def test_gf_library_integrity(model):
    root = greens.library_root(model)
    if not root.exists():
        pytest.skip(f"GF library for {model} not built on this machine")
    manifest = json.loads((root / "manifest.json").read_text())
    assert manifest["npts"] == config.GF_NPTS
    assert manifest["dt_s"] == config.DT
    assert manifest["model"] == model
    # spot-check one file on disk
    from obspy import read

    depth = greens.available_depths(model)[0]
    tr = read(str(greens.gf_file(model, depth, 100, "ZSS")), format="SAC")[0]
    assert tr.stats.npts == config.GF_NPTS
    assert abs(tr.stats.delta - config.DT) < 1e-6


# --- Okada forward ----------------------------------------------------------

def test_wells_coppersmith_anchor():
    length_m, width_m = okada_forward.wells_coppersmith_lw(5.6)
    assert 6.5e3 < length_m < 8.5e3
    assert 5.5e3 < width_m < 7.0e3


def test_forward_displacement_anchors():
    # Mw 5.5 thrust at 10 km: peak |u| of order 1-3 cm (calibration case
    # matching the user's ~2 cm InSAR-detectability intuition)
    r = okada_forward.predicted_displacement(
        mw=5.5, m0_dyne_cm=2.24e24, depth_km=10.0,
        strike=209.0, dip=65.0, rake=78.0)
    assert 0.005 < r["peak_abs_m"] < 0.05
    # a deep Mw 4 is far below detection
    r2 = okada_forward.predicted_displacement(
        mw=4.0, m0_dyne_cm=1.26e22, depth_km=20.0,
        strike=0.0, dip=45.0, rake=90.0)
    assert r2["peak_abs_m"] < 0.001


def test_okada_rake_side_convention():
    # N-S vertical LEFT-LATERAL (rake 0): east side must move north.
    # Guards the strike/dip/rake -> okada4py mapping end to end.
    import numpy as np
    r = okada_forward.predicted_displacement(
        mw=5.5, m0_dyne_cm=2.24e24, depth_km=5.0,
        strike=0.0, dip=89.9, rake=0.0)
    ix = int(np.argmin(np.abs(r["x_km"] - 5.0)))
    iy = int(np.argmin(np.abs(r["y_km"] - 0.0)))
    assert r["un_m"][iy, ix] > 0, "east side of left-lateral must move north"


def test_conjugate_planes_equivalent_far_field():
    # both nodal planes of a DC must give near-identical far-field statics
    import numpy as np
    a = okada_forward.predicted_displacement(
        mw=5.1, m0_dyne_cm=5.6e23, depth_km=8.0, strike=9, dip=50, rake=56)
    b = okada_forward.predicted_displacement(
        mw=5.1, m0_dyne_cm=5.6e23, depth_km=8.0,
        strike=235, dip=51, rake=123)
    L = a["fault"]["length_km"]
    X, Y = np.meshgrid(a["x_km"], a["y_km"])
    far = np.sqrt(X**2 + Y**2) > 3 * L
    corr = np.corrcoef(a["uz_m"][far], b["uz_m"][far])[0, 1]
    assert corr > 0.99, f"conjugate far fields diverge (corr {corr:.3f})"


# --- trigger / gates --------------------------------------------------------

def _event(mag=5.0, lat=-44.4, lon=168.3, depth=5.0, quality="best"):
    return Event(
        public_id="test", origin_time="2026-09-02T06:57:11Z",
        longitude=lon, latitude=lat, depth_km=depth,
        prelim_mag=mag, mag_type="M", locality="test", quality=quality,
    )


def test_processing_floor():
    assert trigger.passes_processing_floor(_event(mag=4.5))[0]
    assert not trigger.passes_processing_floor(_event(mag=3.5))[0]
    assert not trigger.passes_processing_floor(_event(quality="deleted"))[0]
    assert not trigger.passes_processing_floor(_event(lat=-20.0))[0]  # Tonga
    assert not trigger.passes_processing_floor(_event(depth=150.0))[0]  # slab
    assert not trigger.passes_processing_floor(_event(depth=45.0))[0]  # > cap
    assert trigger.passes_processing_floor(_event(depth=25.0))[0]
    # GeoNet placeholder depths are exempt: true depth unknown
    assert trigger.passes_processing_floor(_event(depth=33.0))[0]


def _solution(mw=5.2, vr=70.0, grade="A"):
    return {
        "event": {"public_id": "t", "latitude": -44.4, "longitude": 168.3},
        "preferred": {"mw": mw, "vr": vr},
        "quality": {"passed": grade in ("A", "B"), "grade": grade,
                    "checks": {}},
    }


def test_publish_gates():
    fwd_big = {"peak_abs_m": 0.05, "detectable": True}
    fwd_none = {"peak_abs_m": 0.0001, "detectable": False}

    d = trigger.publish_decision(_solution(mw=5.2), fwd_none, [])
    assert d["publish"], "Mw gate alone should publish"
    d = trigger.publish_decision(_solution(mw=4.6), fwd_big, [])
    assert d["publish"], "displacement gate alone should publish"
    d = trigger.publish_decision(_solution(mw=4.6), fwd_none, [])
    assert not d["publish"]
    d = trigger.publish_decision(_solution(mw=5.2, grade="C"), fwd_big, [])
    assert not d["publish"], "grade C/D blocks publication"
    d = trigger.publish_decision(_solution(mw=5.2, grade="B"), fwd_none, [])
    assert d["publish"], "grade B is emailable"


def test_quality_grades():
    """Grade v2: evidence only — station count and azimuthal gap are NOT
    thresholds (3 well-fitting stations spanning 90 deg can be an A)."""
    gates = invert.quality_gates

    def sol(vr, azimuths, own_vrs, dc=90.0, jk_rot=5.0, edge=False,
            geonet_depth=10.0):
        s = {
            "event": {"depth_km": geonet_depth},
            "preferred": {"vr": vr, "pdc": dc, "depth_km": 10.0},
            "stations_used": [
                {"azimuth": a, "final": {"own_vr": v}}
                for a, v in zip(azimuths, own_vrs)
            ],
            "depth_pick_flags": {"edge_artifact": edge},
        }
        if jk_rot is not None:
            s["jackknife"] = {"n_subsets": len(azimuths),
                              "max_tensor_rotation_deg": jk_rot}
        return s

    # A: strong fit, no passengers, stable, interior depth — and only
    # THREE stations, which the old count-based rubric could not grade A
    q = gates(sol(75, [0, 100, 200], [60, 55, 45]))
    assert q["grade"] == "A" and q["passed"]
    # B: a passenger drags min own VR below the A bar
    q = gates(sol(75, [0, 100, 200], [60, 55, 30]))
    assert q["grade"] == "B" and q["passed"]
    # B: jackknife skipped (too few stations to leave one out)
    q = gates(sol(75, [0, 100, 200], [60, 55, 45], jk_rot=None))
    assert q["grade"] == "B"
    # C: DC below the BSL publishability bar, however good the fit
    q = gates(sol(85, [0, 100, 200], [60, 55, 45], dc=40))
    assert q["grade"] == "C" and not q["passed"]
    # C: grid-edge depth artifact
    q = gates(sol(85, [0, 100, 200], [60, 55, 45], edge=True))
    assert q["grade"] == "C"
    # C: unstable mechanism
    q = gates(sol(85, [0, 100, 200], [60, 55, 45], jk_rot=40))
    assert q["grade"] == "C"
    # D: no 90-degree azimuth pair, however high the VR
    q = gates(sol(90, [0, 30, 60], [70, 70, 70]))
    assert q["grade"] == "D" and not q["passed"]
    # D: a station fitting worse than nothing
    q = gates(sol(90, [0, 100, 200], [70, 70, 5]))
    assert q["grade"] == "D"
    # C: our centroid depth cannot be 30 km from GeoNet's hypocentre and
    # still be published ("it is just not reasonable for the geonet
    # solution and ours to be completely different")
    q = gates(sol(85, [0, 100, 200], [60, 55, 45], geonet_depth=40.0))
    assert q["grade"] == "C" and not q["checks"]["depth_agrees_with_geonet"]
    # ...but a GeoNet PLACEHOLDER depth is not a measurement, so it is
    # never used to judge us
    q = gates(sol(85, [0, 100, 200], [60, 55, 45], geonet_depth=33.0))
    assert q["checks"]["depth_agrees_with_geonet"]


def test_aftershock_throttle():
    now = datetime.now(timezone.utc).isoformat()
    mainshock = {"public_id": "main", "mw": 5.8,
                 "latitude": -44.4, "longitude": 168.3, "published_utc": now}
    fwd_big = {"peak_abs_m": 0.05, "detectable": True}
    # a nearby much-smaller aftershock below the Mw gate is throttled
    d = trigger.publish_decision(_solution(mw=4.8), fwd_big, [mainshock])
    assert not d["publish"]
    # a comparable-size aftershock publishes
    d = trigger.publish_decision(_solution(mw=5.5), fwd_big, [mainshock])
    assert d["publish"]


# --- provenance -------------------------------------------------------------

def test_archived_solution_has_provenance():
    candidates = sorted(config.EVENTS_DIR.glob("*/solution.json"))
    if not candidates:
        pytest.skip("no archived solutions on this machine yet")
    sol = json.loads(candidates[0].read_text())
    prov = sol["provenance"]
    for key in ("velocity_model", "gf_version", "mttime_version",
                "obspy_version", "preferred_rule"):
        assert prov.get(key), f"provenance missing {key}"
    assert sol["quality"]["checks"], "quality gates missing"


# --- station selection rules (2026-09-03 review) ----------------------------

def test_small_events_single_band():
    # small events carry no coherent energy above ~20 s period: 10-50 s only
    assert config.band_candidates(4.0) == [(0.02, 0.10)]
    # mid events: ordered preference, 10-50 s FIRST (VR must not
    # arbitrate across bands below M5.5 - longer periods fit noise)
    assert config.band_candidates(5.0) == [(0.02, 0.10), (0.02, 0.05)]


def test_near_field_magnitude_dependent():
    # small shallow events keep their close stations (info lives there)
    assert config.station_min_dist_km(4.0) == 10.0
    assert config.station_min_dist_km(5.0) == config.MIN_STATION_DIST_KM


def test_selection_thresholds_ordered():
    # dead-channel floor below the strong-signal threshold; the admission
    # floors are real VR fractions; the funnel keeps a majority, not a
    # handful; grade rubric monotonic A -> C
    assert 0 < config.PEAK_NOISE_DEAD < config.PEAK_NOISE_STRONG
    assert 0 < config.ADMIT_VR_MIN_SPARSE <= config.ADMIT_VR_MIN < 100
    assert config.CORE_SIZE_MIN <= config.CORE_SIZE_MAX <= config.PASS1_KEEP_N
    assert config.PASS1_KEEP_N <= config.PASS1_KEEP_MAX
    assert config.MIN_STATIONS_USED <= config.CORE_SIZE_MIN
    r = config.GRADE_RUBRIC
    assert r["A"]["vr"] > r["B"]["vr"] > r["C"]["vr"]
    assert r["A"]["min_own_vr"] > r["B"]["min_own_vr"] > r["C"]["min_own_vr"]
    assert r["A"]["jk_rot_max"] < r["B"]["jk_rot_max"]
    # B is exactly the BSL publishability rule
    assert r["B"]["vr"] == 60.0 and r["B"]["dc"] == 60.0


# --- mechanism comparison metric (J. Townend recipe, 2026-08-20) ------------

def test_min_rotation_identity_and_conjugate():
    from obspy.imaging.beachball import aux_plane
    sdr = (40.0, 50.0, 60.0)
    assert invert.min_rotation_angle_deg(sdr, sdr) < 1e-3
    # the auxiliary plane parameterises the SAME double couple: angle 0
    conj = tuple(aux_plane(*sdr))
    assert invert.min_rotation_angle_deg(sdr, conj) < 0.5


def test_min_rotation_known_values():
    # two vertical strike-slips 30 deg apart in strike: rotation 30 deg
    a = invert.min_rotation_angle_deg((0, 90, 0), (30, 90, 0))
    assert abs(a - 30.0) < 0.5
    # symmetric in argument order
    assert abs(a - invert.min_rotation_angle_deg((30, 90, 0),
                                                 (0, 90, 0))) < 1e-6
