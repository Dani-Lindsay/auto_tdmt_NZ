"""3b. Station selection v7: exhaustive station sets + time-shift sweep + DC-on-plateau depth.

What the 2026-11-17 Arthur's Pass review (2026p866048) taught, now automated:

1. Every station SET is tried, not a drop loop. At each trial depth every
   subset of the pool (>= sweepMinStations, every member own VR >= stationVRFloor,
   >= sweepMinSectors azimuth sectors) is solved exactly from per-station
   normal-equation blocks (subsets.py, the v6 study kernel), so no drop
   order can decide the answer (Fukuyama et al. 1998).
2. Every station's time shift is SWEPT on a grid (-sweepShiftMaxS..+sweepShiftMaxS s)
   and refined by coordinate descent for the best sets. mttime's free
   cross-correlation picks the largest |correlation| peak; a station near a
   nodal plane then looks like bad data (own VR < 0) when a +2..+4 s
   velocity-model delay would fit it at own VR 80-95.
3. The set is chosen at the best-fitting depth: within sweepVRTolerance VR
   points of the maximum, most sectors first, then most stations, then VR
   (dropping stations always raises VR; Triantafyllis et al. 2016).
4. The chosen set with its shifts is inverted at EVERY depth (the smooth
   VR-depth curve the reviewer asked for), and the preferred depth is the
   maximum %DC among depths within depthVRTolerance VR points of the VR
   maximum (depthRule maxdc; the reviewer's EPS207 practice) - VR alone is
   flat over 10+ km for shallow events while %DC peaks sharply.
5. The leave-one-out jackknife keeps the shifts (invert.jackknife(shifts=)).

Cost: pool <= 16 stations -> <= 65k sets x 24 depths vectorised, then
sweepTopSets refinements per depth: ~1-3 min per band on a laptop.
"""
from __future__ import annotations

import itertools
import os
from pathlib import Path

import numpy as np

import config
import subsets
from config import P
from geonet import Event
from invert import NoCoherentSolution, run_inversion, station_id, station_reads, write_mtinv


def _shift_grid() -> list[int]:
    s = int(round(P.invert.sweepShiftMaxS / config.DT))
    return list(range(-s, s + 1))


def _tensor(inv, a, depth: float):
    from mttime.core.tensor import Tensor
    t = Tensor(inv._a2m(a) * 1e20, depth=depth, inversion_type=inv.config.inversion_type,
               components=inv.config.components, station_table=inv.config.station_table.iloc[:1],
               station_VR=np.array([0.0]), total_VR=0.0)
    t.decompose()
    fps = np.asarray(t.fps, dtype=float)
    return float(t.pdc), float(t.mw), fps


class _Blocks:
    """Per-station blocks for every forced shift at one depth: arrays (S, n, ...)."""

    def __init__(self, inv, depth: float, shifts: list[int]):
        n = inv.config.nsta
        S = len(shifts)
        self.depth, self.shifts = depth, shifts
        self.A = np.zeros((S, n, 5, 5)); self.b = np.zeros((S, n, 5)); self.c = np.zeros((S, n))
        self.U = np.zeros_like(self.A); self.u = np.zeros_like(self.b); self.e = np.zeros_like(self.c)
        inv._load_green(depth)
        inv._correlate()                                          # mttime's own shifts, the seed
        self.auto = (inv.config.station_table.ts.values - config.TIME_BEFORE_S).astype(int)
        for si, s in enumerate(shifts):
            inv.config.station_table["ts"] = config.TIME_BEFORE_S + s
            inv._data_to_d(); inv._weights()
            B = subsets.blocks_from_state(inv, depth)
            self.A[si], self.b[si], self.c[si] = B.A, B.b, B.c
            self.U[si], self.u[si], self.e[si] = B.U, B.u, B.e

    def solve(self, sta_idx, shift_idx: np.ndarray):
        """shift_idx (m, k) -> VR (m,), own VR (m, k), a (m, 5)."""
        m, k = shift_idx.shape
        As = np.zeros((m, 5, 5)); bs = np.zeros((m, 5)); cs = np.zeros(m)
        for j, st in enumerate(sta_idx):
            As += self.A[shift_idx[:, j], st]; bs += self.b[shift_idx[:, j], st]; cs += self.c[shift_idx[:, j], st]
        try:
            a = np.linalg.solve(As, bs[..., None])[..., 0]
        except np.linalg.LinAlgError:
            return np.full(m, -np.inf), np.full((m, k), -np.inf), np.zeros((m, 5))
        vr = (1 - (cs - 2 * np.einsum("mj,mj->m", a, bs) + np.einsum("mj,mjk,mk->m", a, As, a)) / cs) * 100
        own = np.zeros((m, k))
        for j, st in enumerate(sta_idx):
            Uj = self.U[shift_idx[:, j], st]; uj = self.u[shift_idx[:, j], st]; ej = self.e[shift_idx[:, j], st]
            own[:, j] = (1 - (ej - 2 * np.einsum("mj,mj->m", a, uj) + np.einsum("mj,mjk,mk->m", a, Uj, a)) / ej) * 100
        return vr, own, a

    def at_auto(self) -> subsets.Blocks:
        """A plain Blocks object with every station at its mttime shift (stage A)."""
        n = len(self.auto)
        si = np.array([self.shifts.index(int(np.clip(s, self.shifts[0], self.shifts[-1]))) for s in self.auto])
        ar = np.arange(n)
        return subsets.Blocks(self.depth, [str(i) for i in ar], self.A[si, ar], self.b[si, ar], self.c[si, ar],
                              self.U[si, ar], self.u[si, ar], self.e[si, ar], si)


def _refine(blk: _Blocks, sta_idx: list[int], sh: list[int], floor: float, passes: int = 4):
    """Coordinate descent on each station's shift index; objective VR, every own VR >= floor."""
    sh = list(sh)
    for _ in range(passes):
        changed = False
        for j in range(len(sta_idx)):
            grid = np.array([[*sh[:j], g, *sh[j + 1:]] for g in range(len(blk.shifts))])
            vr, own, _ = blk.solve(sta_idx, grid)
            vr = np.where(own.min(axis=1) >= floor, vr, -np.inf)
            g = int(np.argmax(vr))
            if np.isfinite(vr[g]) and g != sh[j]:
                sh[j] = g; changed = True
        if not changed:
            break
    vr, own, a = blk.solve(sta_idx, np.array([sh]))
    return sh, float(vr[0]), own[0], a[0]


def pick_depth(profile: list[dict]) -> tuple[float, str]:
    """depthRule maxdc: the highest %DC among depths within depthVRTolerance VR
    points of the VR maximum (plateau); maxvr: mttime's argmax VR."""
    best_vr = max(r["vr"] for r in profile)
    if P.invert.depthRule == "maxvr":
        z = max(profile, key=lambda r: r["vr"])["depth_km"]
        return z, "maximum VR"
    plateau = [r for r in profile if r["vr"] >= best_vr - P.invert.depthVRTolerance]
    z = max(plateau, key=lambda r: (r["pdc"], r["vr"]))["depth_km"]
    return z, (f"maximum %DC among depths within {P.invert.depthVRTolerance:g} VR points of the "
               f"VR maximum ({best_vr:.1f})")


def select(event: Event, pool: list[dict], depths: list[float], band_dir: Path, green_dir: Path):
    """Returns (inversion, used, rejected, info). Raises NoCoherentSolution."""
    from mttime import Configure, Inversion

    rows = [dict(r) for r in pool]
    n = len(rows)
    floor = P.invert.stationVRFloor
    min_n = max(P.invert.minStations, P.invert.sweepMinStations)
    if n < min_n:
        raise NoCoherentSolution("sweep", f"only {n} usable stations (minimum {min_n})", [])
    shifts = _shift_grid()
    sec = np.array([r["sector"] for r in rows])
    onehot = np.eye(P.station.sectors, dtype=int)[sec]

    cwd = os.getcwd(); os.chdir(band_dir)
    try:
        cfg = Configure(path_to_file=str(write_mtinv(event, rows, depths, band_dir, green_dir, correlate=True)))
        inv = Inversion(config=cfg)
        inv._load_data()
        inv.config.station_table.insert(loc=6 + inv.config.ncomp, column="weights", value=None)

        masks = np.array(list(itertools.product([False, True], repeat=n)), dtype=bool)
        masks = masks[masks.sum(1) >= min_n]
        nsec = (masks.astype(int) @ onehot > 0).sum(1)
        min_sec = P.invert.sweepMinSectors
        if (nsec >= min_sec).sum() == 0:                 # one-sided (offshore) pool
            min_sec = max(2, int(nsec.max()))
            print(f"sweep: pool spans {nsec.max()} sectors, minimum sectors relaxed to {min_sec}")

        # ---- stage A: every set at every depth with mttime's shifts -------------
        blocks: dict[float, _Blocks] = {}
        stage_a: dict[float, list[tuple[float, np.ndarray]]] = {}
        best_a = -np.inf
        for depth in depths:
            blk = _Blocks(inv, depth, shifts)
            blocks[depth] = blk
            B = blk.at_auto()
            vr = np.empty(len(masks)); minown = np.empty(len(masks))
            for i0 in range(0, len(masks), 200_000):
                Mc = masks[i0:i0 + 200_000]
                _, v, o = subsets.solve_many(B, Mc)
                vr[i0:i0 + len(Mc)] = v
                minown[i0:i0 + len(Mc)] = np.nanmin(np.where(Mc, o, np.inf), axis=1)
            ok = (minown >= floor) & (nsec >= min_sec) & np.isfinite(vr)
            top = np.flatnonzero(ok)[np.argsort(-vr[ok])[:P.invert.sweepTopSets]]
            stage_a[depth] = [(float(vr[j]), masks[j]) for j in top]
            if len(top):
                best_a = max(best_a, float(vr[top[0]]))
            print(f"  sweep z {depth:5.1f}: {int(ok.sum())} sets pass, best VR "
                  f"{vr[top[0]] if len(top) else float('nan'):5.1f}")
        if not np.isfinite(best_a):
            raise NoCoherentSolution("sweep", f"no station set passes own VR >= {floor:g} with "
                                     f">= {min_sec} sectors at any depth", [])

        # ---- stage B: refine shifts of the top sets near the best depths -------
        cands = []
        for depth in depths:
            if not stage_a[depth] or stage_a[depth][0][0] < best_a - 10.0:
                continue
            blk = blocks[depth]
            for vr0, mask in stage_a[depth]:
                idx = [int(i) for i in np.flatnonzero(mask)]
                sh0 = [shifts.index(int(np.clip(blk.auto[i], shifts[0], shifts[-1]))) for i in idx]
                sh, vr, own, a = _refine(blk, idx, sh0, floor)
                cands.append(dict(depth=depth, idx=idx, sh=sh, vr=vr, nsec=int((onehot[idx].sum(0) > 0).sum()),
                                  n=len(idx), own=own))
        vr_max = max(c["vr"] for c in cands)
        near = [c for c in cands if c["vr"] >= vr_max - P.invert.sweepVRTolerance]
        near.sort(key=lambda c: (c["nsec"], c["n"], c["vr"]), reverse=True)

        # ---- the chosen set over EVERY depth (fixed shifts) + depth rule ------
        chosen = None
        for c in near[:5]:
            profile = []
            for depth in depths:
                blk = blocks[depth]
                vr, own, a = blk.solve(c["idx"], np.array([c["sh"]]))
                pdc, mw, fps = _tensor(inv, a[0], depth)
                profile.append(dict(depth_km=depth, vr=float(vr[0]), pdc=pdc, mw=mw))
            z, rule = pick_depth(profile)
            dc_at_z = next(r["pdc"] for r in profile if r["depth_km"] == z)
            if chosen is None:
                chosen = (c, profile, z, rule)
            if dc_at_z >= P.invert.dcMinPublish:
                chosen = (c, profile, z, rule)
                break
        c, profile, z_pref, rule = chosen
        used = [rows[i] for i in c["idx"]]
        shift_s = {station_id(rows[i]): shifts[g] * config.DT for i, g in zip(c["idx"], c["sh"])}
        print(f"sweep: chosen {len(used)} stations / {c['nsec']} sectors at {c['depth']:g} km "
              f"(VR {c['vr']:.1f} of max {vr_max:.1f}); preferred depth {z_pref:g} km by {rule}")
        print("  shifts: " + " ".join(f"{k.split('.')[1]}{v:+.0f}" for k, v in shift_s.items()))

        # ---- the real mttime run of the chosen set, every depth, shifts fixed --
        inv_final = run_inversion(write_mtinv(event, used, depths, band_dir, green_dir,
                                              correlate=False, shifts=shift_s))
        zs = [float(m.depth) for m in inv_final.moment_tensors]
        inv_final.preferred_tensor_id = zs.index(z_pref)
    finally:
        os.chdir(cwd)

    pref = inv_final.moment_tensors[inv_final.preferred_tensor_id]
    reads = station_reads(pref, used)
    for r in used:
        own, shift = reads[station_id(r)]
        r["station_vr"], r["zcor_s"] = round(own, 1), round(shift, 1)
    chosen_ids = {station_id(r) for r in used}
    # own VR of each unused station under the chosen solution at the preferred depth
    blk = blocks[z_pref]
    rejected = []
    for i, r in enumerate(rows):
        if station_id(r) in chosen_ids:
            continue
        best_own = -np.inf
        for g in range(len(shifts)):
            grid = np.array([[*c["sh"], g]])
            _, own, _ = blk.solve(c["idx"] + [i], grid)
            best_own = max(best_own, float(own[0, -1]))
        rejected.append({**r, "station": station_id(r),
                         "reason": f"sweep: not in the best set (own VR {best_own:.0f} at best shift, "
                                   f"{z_pref:g} km)"})
    info = {
        "method": "sweep v7: every station set x shift grid; set chosen within "
                  f"{P.invert.sweepVRTolerance:g} VR of the maximum by sectors > stations > VR; "
                  "depth = " + rule,
        "chosen_at_depth_km": c["depth"], "set_vr_at_choice": round(c["vr"], 1),
        "vr_max_any_set": round(vr_max, 1), "n_candidate_sets": len(cands),
        "shifts_s": shift_s, "preferred_depth_km": z_pref, "depth_rule": rule,
        "profile": [{k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()} for r in profile],
    }
    return inv_final, used, rejected, info
