"""The v5 rules, each pinned to a small example. If one of these fails,
a published rule has silently changed — check auto_tdmt.cfg first."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import config  # noqa: E402
import invert  # noqa: E402
import params  # noqa: E402
import waveforms  # noqa: E402
from config import P  # noqa: E402


# --- the parameter file ---------------------------------------------------

def test_cfg_and_defaults_agree():
    """Every key in auto_tdmt.cfg exists in params.DEFAULTS and vice versa,
    so a typo in the cfg cannot silently fall back to a default."""
    in_cfg = params.cfg_keys()
    in_code = {(t, k) for t, d in params.DEFAULTS.items() for k in d}
    assert in_cfg == in_code, (
        f"cfg only: {in_cfg - in_code}; code only: {in_code - in_cfg}")


def test_auto_resolves_and_types_are_kept(tmp_path):
    cfg = tmp_path / "t.cfg"
    cfg.write_text(
        "auto_tdmt.invert.depthWindowKm = 20   # x\n"
        "auto_tdmt.station.periodsMediumS = 10 50, 20 60\n"
        "auto_tdmt.invert.degree = auto\n"
        "auto_tdmt.publish.minGrade = A\n")
    Q = params.load(cfg)
    assert Q.invert.depthWindowKm == 20.0 and isinstance(Q.invert.degree, int)
    assert Q.station.periodsMediumS == [[10.0, 50.0], [20.0, 60.0]]
    assert Q.publish.minGrade == "A"
    assert Q.station.snrMin == params.DEFAULTS["station"]["snrMin"]


def test_unknown_parameter_is_an_error(tmp_path):
    cfg = tmp_path / "t.cfg"
    cfg.write_text("auto_tdmt.invert.depthWinowKm = 20\n")  # typo
    with pytest.raises(AssertionError):
        params.load(cfg)


def test_band_menu_is_ordered_and_magnitude_binned():
    assert P.bands_hz(4.0) == [(0.02, 0.10)]
    assert P.bands_hz(5.0)[0] == (0.02, 0.10)          # 10-50 s first
    assert all(hi <= 0.05 for _, hi in P.bands_hz(6.0))  # no 10 s energy
    assert config.band_tag((0.01, 0.05)) == "band_20-100s"


def test_window_cap_fits_the_greens_functions():
    assert config.INV_NPTS <= config.GF_NPTS
    assert config.INV_NPTS - config.TIME_BEFORE_S == P.station.maxWindowS == 200


# --- task 2: selection ------------------------------------------------------

def _row(name, az, snr, dist=100.0):
    return {"network": "NZ", "station": name, "location": "10",
            "azimuth": az, "snr_med": snr, "distance_km": dist}


def test_sector_round_robin_balances_when_it_can():
    # 6 stations in one sector, 1 each in two others, cap 4: the lone
    # stations get seats before the fourth-best of the crowded sector
    rows = ([_row(f"N{i}", 10.0, 90 - i) for i in range(6)]
            + [_row("E", 100.0, 20.0), _row("S", 200.0, 30.0)])
    sel, left = waveforms.select_by_sector(rows, max_stations=4)
    names = {r["station"] for r in sel}
    assert {"E", "S", "N0", "N1"} == names
    assert len(left) == 4


def test_sector_round_robin_never_starves_one_sided_geometry():
    # offshore: everything in one sector, cap 4 -> the best four, not two
    rows = [_row(f"N{i}", 10.0, 90 - i) for i in range(6)]
    sel, _ = waveforms.select_by_sector(rows, max_stations=4)
    assert [r["station"] for r in sel] == ["N0", "N1", "N2", "N3"]


# --- task 3: the loop, depth, grade ----------------------------------------

def test_station_verdict_fixed_floor_and_shift_cap():
    # a FIXED own-VR floor (25): the relative Clinton rule is off by
    # default because without sector refilling it ratchets to 3 stations
    assert invert.station_verdict(own_vr=26, shift_s=2, overall_vr=80) is None
    assert "poor fit" in invert.station_verdict(own_vr=20, shift_s=2,
                                                overall_vr=80)
    assert P.invert.stationVRDrop == 0
    # a chance alignment is rejected however well it "fits"
    assert "time shift" in invert.station_verdict(own_vr=90, shift_s=41,
                                                  overall_vr=60)


def test_depth_window_bounds_the_search_but_placeholders_get_everything():
    from geonet import Event

    def ev(depth):
        return Event(public_id="t", origin_time="2026-09-02T06:57:11Z",
                     longitude=168.3, latitude=-44.4, depth_km=depth,
                     prelim_mag=4.5, mag_type="M", locality="t", quality="best")
    try:
        d = invert.search_depths(ev(26.0), config.GF_MODELS[0])
    except AssertionError:
        pytest.skip("GF library not built on this machine")
    assert min(d) >= 26 - P.invert.depthWindowKm - 1e-6
    assert max(d) <= 26 + P.invert.depthWindowKm + 1e-6
    assert len(invert.search_depths(ev(33.0), config.GF_MODELS[0])) == \
        len(config.GF_DEPTHS_KM)


def test_depth_uncertainty_is_the_ten_percent_band():
    rows = [{"depth_km": d, "vr": v} for d, v in
            ((10, 50), (12, 60), (14, 66), (16, 63), (18, 55))]
    assert invert.depth_uncertainty(rows, pct=10) == (12, 16)


def test_ingv_grade_table():
    g = invert.grade_ingv
    # the VR bar falls as station count rises
    assert g(65, 3, 90, 5) == "C"    # A unreachable with 3 stations
    assert g(72, 3, 90, 5) == "B"
    assert g(50.8, 4, 97, 24.6) == "B"   # 2026p669681: "I'd be confident"
    assert g(62, 4, 90, 5) == "A"
    assert g(45, 6, 90, 5) == "B" and g(62, 6, 90, 5) == "A"
    assert g(35, 10, 90, 5) == "B" and g(52, 10, 90, 5) == "A"
    # dropping stations cannot buy a grade: 10 stations at 30 = B,
    # 2 stations at 65 = C
    assert g(30, 10, 90, 5) == "B" and g(65, 2, 90, 5) == "C"
    # below the C bar is D
    assert g(12, 6, 90, 5) == "D"
    # BSL rule: DC < 60 caps at C; unstable mechanism caps at C
    assert g(70, 6, 40, 5) == "C" and g(70, 6, 90, 40) == "C"
    # jackknife impossible (3 stations) is allowed for B
    assert g(72, 3, 90, None) == "B"


def test_quality_gates_shape():
    sol = {
        "event": {"depth_km": 8.0},
        "preferred": {"vr": 50.8, "pdc": 97.0, "depth_km": 6.0},
        "stations_used": [{"azimuth": a, "station_vr": v}
                          for a, v in ((0, 60), (90, 40), (180, 26), (300, 55))],
        "depth_search": [{"depth_km": d, "vr": 40} for d in (2, 4, 6, 8, 10)],
        "depth_pick_flags": {"depth_range_km": [4, 8]},
        "jackknife": {"n_subsets": 4, "max_tensor_rotation_deg": 24.6},
    }
    q = invert.quality_gates(sol)
    assert q["grade"] == "B" and q["passed"]
    assert q["n_stations_used"] == 4 and q["min_own_vr"] == 26
    assert set(q["checks"]) == {"vr_grade", "dc_floor", "jackknife_stable",
                                "min_stations"}


def test_reason_classes_cover_every_task():
    rc = invert.reason_class
    assert rc("no data: gaps") == "nodata"
    assert rc("low SNR: median component SNR 1.4 < 2") == "snr"
    assert rc("amplitude outlier: 12x") == "amp"
    assert rc("not selected: 16 seats filled") == "not_selected"
    assert rc("poor fit: own VR 10 < 25") == "fit"
    assert rc("time shift: +41 s exceeds 8 s") == "shift"
    assert rc("aborted: only 2 stations") == "abort"
