"""Anchor tests for station selection v4 (the funnel).

Every helper here is pure — no mttime, no network — so the selection
LOGIC is pinned independently of the inversion. The cases are the real
events that motivated each rule (see docs/REVIEW_LEARNINGS.md).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config  # noqa: E402
import invert  # noqa: E402


def _row(station, az, dist=100.0, own1=50.0, own2=None, zok=True,
         tags=None):
    return {
        "network": "NZ", "station": station, "location": "10",
        "azimuth": az, "distance_km": dist, "pk_n": 6.0,
        "tags": tags or [], "tier": "demoted" if tags else "trusted",
        "pass1": {"own_vr": own1, "zcor_ok": zok},
        "pass2": {"own_vr": own1 if own2 is None else own2,
                  "zcor_ok": zok},
    }


# --- depth pick: the grid-edge guard ---------------------------------------

def test_grid_edge_guard_rejects_single_point_edge_maximum():
    # 2026p508890: the VR maximum sat alone on the 58 km grid edge while
    # the physical peak (with DC 88-96) was interior at 8 km
    vrs = [8.8, 7.2, 6.6, 9.7, 11.0, 18.3, 8.1, 14.4, 22.7]
    i, guarded = invert.grid_edge_reference(vrs)
    assert guarded and i == 5 and vrs[i] == 18.3


def test_grid_edge_guard_keeps_a_real_plateau_at_the_edge():
    # a shallow crustal event legitimately prefers the shallowest depths:
    # the edge maximum belongs to a plateau, so it stands
    vrs = [70.0, 68.0, 66.0, 40.0, 20.0]
    i, guarded = invert.grid_edge_reference(vrs)
    assert not guarded and i == 0


def test_grid_edge_guard_keeps_edge_when_nothing_interior_is_close():
    vrs = [5.0, 6.0, 4.0, 3.0, 40.0]
    i, guarded = invert.grid_edge_reference(vrs)
    assert not guarded and i == 4


def test_pick_preferred_still_honours_the_contiguous_plateau():
    # 2026p660272: a disconnected deep lobe grazing the tolerance must not
    # steal the pick on DC
    vr_pdc = [(60.0, 20.0), (58.0, 30.0), (20.0, 10.0), (57.0, 99.0)]
    assert invert.pick_preferred(vr_pdc) == 1


# --- geometry ---------------------------------------------------------------

def test_azimuth_pair_needs_ninety_degrees():
    assert invert.azimuth_pair_ok([0.0, 95.0])
    assert not invert.azimuth_pair_ok([0.0, 30.0, 60.0])
    assert invert.azimuth_pair_ok([350.0, 95.0])       # wraparound
    assert not invert.azimuth_pair_ok([10.0])


# --- the funnel: ranking, keeping the majority, the core --------------------

def test_rank_pass1_sinks_chance_alignments():
    """A station that only fits because its trace slid into a chance
    alignment must rank below every honest one, however good its VR."""
    rows = [_row("GOOD", 0, own1=40.0),
            _row("CHEAT", 90, own1=95.0, zok=False),
            _row("OK", 180, own1=55.0)]
    ranked = invert.rank_pass1(rows)
    assert [r["station"] for r in ranked] == ["OK", "GOOD", "CHEAT"]


def test_keep_after_pass1_keeps_the_majority_and_fills_sectors():
    rows = [_row(f"S{i}", i * 5.0, own1=90.0 - i) for i in range(14)]
    # one distant station alone in its sector, weaker than the top ten
    rows.append(_row("FAR", 200.0, own1=25.0))
    keep = invert.keep_after_pass1(invert.rank_pass1(rows))
    assert len(keep) <= config.PASS1_KEEP_MAX
    assert "FAR" in [r["station"] for r in keep]     # sector fill
    assert len([r for r in keep if r["station"] != "FAR"]) \
        == config.PASS1_KEEP_N


def test_keep_after_pass1_never_starves_below_the_minimum():
    rows = [_row(f"S{i}", i * 90.0, own1=2.0) for i in range(4)]
    keep = invert.keep_after_pass1(invert.rank_pass1(rows))
    assert len(keep) >= config.MIN_STATIONS_USED


def test_choose_core_enforces_ninety_degree_span():
    """Strong stations bunched on one side cannot resolve a mechanism: the
    core must reach past them for a weaker station that spans, even though
    the fit alone would never have chosen it."""
    rows = [_row(f"NEAR{i}", i * 10.0, own2=80.0 - i) for i in range(7)]
    rows.append(_row("WIDE", 150.0, own2=45.0))
    core, note = invert.choose_core(rows)
    assert len(core) == config.CORE_SIZE_MAX
    assert invert.azimuth_pair_ok([r["azimuth"] for r in core])
    assert "WIDE" in [r["station"] for r in core] and "swapped" in note


def test_choose_core_leaves_a_spanning_core_alone():
    rows = [_row("A", 0.0, own2=80.0), _row("B", 120.0, own2=70.0),
            _row("C", 240.0, own2=60.0)]
    core, note = invert.choose_core(rows)
    assert [r["station"] for r in core] == ["A", "B", "C"] and note == ""


def test_choose_core_relaxes_then_gives_up():
    rows = [_row("A", 0.0, own2=22.0), _row("B", 120.0, own2=21.0),
            _row("C", 240.0, own2=20.5)]
    core, note = invert.choose_core(rows)
    assert len(core) == 3 and "relaxed" in note
    core, note = invert.choose_core(
        [_row("A", 0.0, own2=5.0), _row("B", 120.0, own2=4.0)])
    assert core == [] and note == "no coherent core"


# --- the funnel: earning a seat ---------------------------------------------

def test_admission_accepts_a_fitting_station():
    ok, reason = invert.admission_verdict(
        52.0, 2.0, True, 55.0, 56.0, True, 4, [], sparse=False)
    assert ok and reason.startswith("admitted:") and "new sector" in reason


def test_admission_rejects_chance_alignment_however_good_the_fit():
    ok, reason = invert.admission_verdict(
        88.0, 41.0, False, 55.0, 60.0, True, 4, [], sparse=False)
    assert not ok and invert.reason_class(reason) == "not_admitted"
    assert "chance alignment" in reason


def test_admission_lets_a_weak_station_in_for_azimuth_coverage():
    """Geometry is the scarcer commodity: a station filling an empty
    sector is worth admitting on a weaker fit (2026p189537)."""
    ok, _ = invert.admission_verdict(
        15.0, 1.0, True, 50.0, 49.0, True, 3, [], sparse=False)
    assert ok
    # ...but not when it adds nothing geometrically
    ok, reason = invert.admission_verdict(
        15.0, 1.0, True, 50.0, 49.0, False, 3, [], sparse=False)
    assert not ok and "own VR" in reason


def test_admission_sparse_core_relaxes_the_floor():
    args = (22.0, 1.0, True, 50.0, 50.0, False, 3, [])
    assert not invert.admission_verdict(*args, sparse=False)[0]
    assert invert.admission_verdict(*args, sparse=True)[0]


def test_admission_cluster_surplus_must_earn_a_joint_gain():
    """A Ruapehu-ring station that merely duplicates its neighbours is
    admitted only if the joint fit actually improves."""
    ok, reason = invert.admission_verdict(
        45.0, 1.0, True, 55.0, 54.0, False, 4, ["cluster_surplus"],
        sparse=False)
    assert not ok and "cluster surplus" in reason
    ok, _ = invert.admission_verdict(
        45.0, 1.0, True, 55.0, 56.0, False, 4, ["cluster_surplus"],
        sparse=False)
    assert ok


def test_admission_respects_the_station_cap():
    ok, reason = invert.admission_verdict(
        80.0, 1.0, True, 60.0, 61.0, True, config.MAX_USED_STATIONS, [],
        sparse=False)
    assert not ok and "cap" in reason


def test_admission_refuses_a_station_that_spoils_the_joint_fit():
    ok, reason = invert.admission_verdict(
        45.0, 1.0, True, 55.0, 45.0, False, 4, [], sparse=False)
    assert not ok and "joint VR" in reason


# --- the shared reason vocabulary ------------------------------------------

def test_reason_class_covers_every_producer():
    cases = {
        "no data: download failed: 204": "nodata",
        "no data: gaps or <3 components": "nodata",
        "dead channel: peak/noise 0.9 < 1.2": "dead",
        "amplitude outlier: peak x dist 1e-3 vs network median 1e-5": "amp",
        "not admitted: own VR 12 < 30": "not_admitted",
        "not admitted: chance alignment (zcor 41 s > 8 s)": "not_admitted",
        "anti-fitting: own VR -38 at 12 km": "antifit",
        "aborted: final VR 8 < 20": "abort",
    }
    for reason, cls in cases.items():
        assert invert.reason_class(reason) == cls, reason


# --- the abort record -------------------------------------------------------

def test_no_solution_record_is_marked_and_publishes_nothing():
    class _Ev:
        public_id = "2026pTEST"
        locality = "15 km east of Nowhere"

        def to_dict(self):
            return {"public_id": self.public_id, "prelim_mag": 4.0,
                    "latitude": -41.0, "longitude": 174.0, "depth_km": 12.0,
                    "origin_time": "2026-01-01T00:00:00Z"}

    pool = [_row("A", 0.0), _row("B", 120.0)]
    rec = invert.no_solution_record(
        _Ev(), pool, [], "nz_south_ristau2008", (0.02, 0.10),
        "pass2", "no coherent solution: survey VR 8, majority VR 11",
        {"pass1": {"vr": 8.0}, "pass2": {"vr": 11.0}})
    assert not config.is_solved(rec)
    assert rec["quality"]["grade"] == "X" and not rec["quality"]["passed"]
    assert not rec["publish_decision"]["publish"]
    assert rec["stations_used"] == []
    assert rec["abort"]["best_vr"] == 11.0
    # every pool station is accounted for, with the abort reason
    assert len(rec["stations_dropped"]) == len(pool)
    assert all(invert.reason_class(d["reason"]) == "abort"
               for d in rec["stations_dropped"])
    assert config.no_solution_dir_name("2026pTEST", "15 km east of Nowhere") \
        == "2026pTEST_NOSOL_15-km-east-of-Nowhere"


def test_code_version_is_recorded():
    assert isinstance(config.code_version(), str)
    assert config.SELECTION_VERSION
