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
    gates = invert.quality_gates
    def sol(vr, stations, azimuths):
        return {
            "preferred": {"vr": vr, "depth_km": 10.0},
            "stations_used": [
                {"azimuth": a} for a in azimuths[:stations]
            ],
        }
    # A: strong VR, good coverage
    q = gates(sol(75, 6, [0, 60, 120, 180, 240, 300]))
    assert q["grade"] == "A" and q["passed"]
    # B: decent VR, 3 stations, wide gap
    q = gates(sol(65, 3, [0, 90, 180]))
    assert q["grade"] == "B" and q["passed"]
    # C: 2 stations
    q = gates(sol(80, 2, [0, 90]))
    assert q["grade"] == "C" and not q["passed"]
    # D: single station
    q = gates(sol(90, 1, [0]))
    assert q["grade"] == "D" and not q["passed"]


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
