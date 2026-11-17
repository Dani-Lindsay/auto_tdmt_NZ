"""Selection v7: the depth rule, forced time shifts in mtinv.in, and the
cfg keys that drive the sweep."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import config  # noqa: E402
import invert  # noqa: E402
import sweep_select  # noqa: E402
from config import P  # noqa: E402
from geonet import Event  # noqa: E402


def _profile():
    # VR flat 84-90 from 1.5 to 14 km, %DC peaking at 4 km (2026p866048)
    return [dict(depth_km=z, vr=v, pdc=d, mw=4.8) for z, v, d in (
        (1.0, 84.4, 11), (2.0, 85.9, 69), (3.0, 87.4, 68), (4.0, 88.5, 99),
        (5.0, 89.3, 80), (8.0, 90.3, 49), (12.0, 89.0, 37), (20.0, 81.7, 31))]


def test_depth_rule_maxdc_picks_the_dc_peak_on_the_plateau(monkeypatch):
    monkeypatch.setattr(P.invert, "depthRule", "maxdc")
    monkeypatch.setattr(P.invert, "depthVRTolerance", 5.0)
    z, rule = sweep_select.pick_depth(_profile())
    assert z == 4.0 and "DC" in rule


def test_depth_rule_maxvr_is_mttimes_argmax(monkeypatch):
    monkeypatch.setattr(P.invert, "depthRule", "maxvr")
    z, _ = sweep_select.pick_depth(_profile())
    assert z == 8.0


def test_depth_rule_plateau_width_matters(monkeypatch):
    monkeypatch.setattr(P.invert, "depthRule", "maxdc")
    monkeypatch.setattr(P.invert, "depthVRTolerance", 1.0)   # only 89.3-90.3 qualify
    z, _ = sweep_select.pick_depth(_profile())
    assert z == 5.0


def test_write_mtinv_forces_shifts_and_turns_correlation_off(tmp_path):
    ev = Event("test", "2026-11-17T10:30:48Z", 171.776, -43.014, 5.0, 5.7, "M", "x", "best")
    rows = [dict(network="NZ", station="INZ", location="10", distance_km=42.0,
                 gf_distance_km=40.0, azimuth=320.0, longitude=171.4, latitude=-42.7,
                 window_end_s=76),
            dict(network="NZ", station="LTZ", location="10", distance_km=48.0,
                 gf_distance_km=50.0, azimuth=58.0, longitude=172.3, latitude=-42.8,
                 window_end_s=76)]
    path = invert.write_mtinv(ev, rows, [4.0], tmp_path, tmp_path / "greens",
                              shifts={"NZ.INZ.10": 2.0, "NZ.LTZ.10": -1.0})
    text = path.read_text().splitlines()
    assert any(line.startswith("correlate") and line.split()[-1] == "0" for line in text)
    table = [line.split() for line in text[12:]]
    cols = table[0]
    ts = {r[cols.index("station")]: int(r[cols.index("ts")]) for r in table[1:]}
    assert ts == {"NZ.INZ.10": config.TIME_BEFORE_S + 2, "NZ.LTZ.10": config.TIME_BEFORE_S - 1}


def test_sweep_parameters_exist():
    assert P.invert.selection in ("sweep", "clinton")
    assert P.invert.sweepMinStations >= P.invert.minStations
    assert 0 < P.invert.sweepVRTolerance <= P.invert.depthVRTolerance
    assert config.SELECTION_VERSION == "v7"
