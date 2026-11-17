"""Station-set algebra on top of mttime (selection v7 kernel, from the v6 study).

Why any station set can be solved without re-running mttime
-----------------------------------------------------------
At one trial depth mttime (inversion.py, _single_inversion) solves the
weighted least-squares problem

    a = (G^T W G)^-1 G^T W d                                   (1)
    VR_total   = 1 - sum(w (d - G a)^2) / sum(w d^2)           (2)
    VR_station = 1 - sum_i(m (d - G a)^2) / sum_i(m d^2)       (3)

(m = 1 where w != 0), and the time shift of station i (_correlate)
cross-correlates station i's data with station i's Green's functions
only. So at a fixed depth station i contributes FIXED blocks

    A_i = G_i^T W_i G_i,  b_i = G_i^T W_i d_i,  c_i = d_i^T W_i d_i
    U_i = G_i^T M_i G_i,  u_i = G_i^T M_i d_i,  e_i = d_i^T M_i d_i

and for ANY set S

    a_S = (sum_S A_i)^-1 (sum_S b_i)                           (4)
    VR_total(S)  = 1 - (C - 2 a.B + a.A a) / C                 (5)
    VR_i(S)      = 1 - (e_i - 2 a.u_i + a.U_i a) / e_i         (6)

with A, B, C the sums over S. mttime's distance weights are distance /
min(distance in the set): a common factor, which cancels in (4) and (5),
so blocks built once on the whole pool are exact for every subset.

All mttime-specific work (reading Green's functions, time shifts, data
vector, weights) is done by mttime's own methods in blocks_from_mttime;
the rest of this module only sums blocks. tests/test_subsets.py checks
(4)-(6) against mttime itself.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SINGULAR_COND = 1e12   # cond(A) above this: the set cannot resolve the tensor


@dataclass
class Blocks:
    """Per-station normal-equation blocks at one trial depth."""
    depth: float
    stations: list[str]
    A: np.ndarray        # (n, k, k) weighted G^T W G
    b: np.ndarray        # (n, k)    weighted G^T W d
    c: np.ndarray        # (n,)      weighted d^T W d
    U: np.ndarray        # (n, k, k) unweighted (fitted components only)
    u: np.ndarray        # (n, k)
    e: np.ndarray        # (n,)
    ts: np.ndarray       # (n,) mttime ts (samples) at this depth


def blocks_from_mttime(inv, depth: float) -> Blocks:
    """Run mttime's own per-depth preparation (Green's functions, time
    shift, data vector, weights), then slice it per station."""
    inv._load_green(depth)
    if inv.config.correlate:
        inv._correlate()
    inv._data_to_d()
    inv._weights()
    return blocks_from_state(inv, depth)


def blocks_from_state(inv, depth: float) -> Blocks:
    """Slice the G, d, w that mttime has just used at ``depth`` (called
    straight after its _single_inversion, so nothing is recomputed)."""
    G, d, w = inv.G, inv.d, inv.w
    m = (w != 0).astype(float)
    n, k = inv.config.nsta, G.shape[1]
    A = np.zeros((n, k, k)); b = np.zeros((n, k)); c = np.zeros(n)
    U = np.zeros((n, k, k)); u = np.zeros((n, k)); e = np.zeros(n)
    for i in range(n):
        s = slice(inv.config.index1[i, 0], inv.config.index2[i, -1])
        Gi, di, wi, mi = G[s], d[s], w[s], m[s]
        A[i] = Gi.T @ (wi[:, None] * Gi); b[i] = Gi.T @ (wi * di); c[i] = wi @ (di * di)
        U[i] = Gi.T @ (mi[:, None] * Gi); u[i] = Gi.T @ (mi * di); e[i] = mi @ (di * di)
    assert np.all(e > 0), "a station has no fitted data"
    return Blocks(float(depth), list(inv.config.station_table.station), A, b, c,
                  U, u, e, inv.config.station_table.ts.values.copy())


def solve_subset(B: Blocks, idx) -> tuple[np.ndarray, float, np.ndarray]:
    """Equations (4)-(6) for the stations idx: (a, VR_total %, VR_i % per idx)."""
    idx = np.asarray(idx, dtype=int)
    A = B.A[idx].sum(0); b = B.b[idx].sum(0); c = B.c[idx].sum()
    a = np.linalg.solve(A, b)
    vr = (1 - (c - 2 * a @ b + a @ A @ a) / c) * 100
    return a, float(vr), station_vr(B, a, idx)


def set_vr(B: Blocks, idx) -> float:
    """VR_total of one set, -inf when its normal matrix is singular at
    this depth (the set cannot resolve the tensor there)."""
    mask = np.zeros(len(B.c), dtype=bool)
    mask[np.asarray(idx, dtype=int)] = True
    return float(solve_many(B, mask[None, :])[1][0])


def station_vr(B: Blocks, a: np.ndarray, idx) -> np.ndarray:
    """Equation (6) for stations idx under solution a (in or out of the set)."""
    idx = np.asarray(idx, dtype=int)
    quad = np.einsum("j,ijk,k->i", a, B.U[idx], a)
    return (1 - (B.e[idx] - 2 * B.u[idx] @ a + quad) / B.e[idx]) * 100


def solve_many(B: Blocks, masks: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Equations (4)-(6) for many sets at once. masks: (m, n) bool.
    Returns a (m, k), VR_total (m,), VR_i (m, n) with nan outside each set.
    Same algebra as solve_subset, written as matrix products."""
    n, k = B.b.shape
    M = masks.astype(float)
    A = (M @ B.A.reshape(n, k * k)).reshape(-1, k, k)
    b = M @ B.b
    c = M @ B.c
    # a set whose normal matrix is singular cannot resolve the tensor
    # (mttime's gaussj would fail on it too): it is solved as nothing and
    # marked VR = -inf so it can never pass a gate
    try:
        a = np.linalg.solve(A, b[..., None])[..., 0]
        bad = np.zeros(len(A), dtype=bool)
    except np.linalg.LinAlgError:
        bad = ~(np.linalg.cond(A) < SINGULAR_COND)
        a = np.zeros((len(A), k))
        if (~bad).any():
            a[~bad] = np.linalg.solve(A[~bad], b[~bad][..., None])[..., 0]
    aa = (a[:, :, None] * a[:, None, :]).reshape(-1, k * k)       # a_j a_k
    vr = (1 - (c - 2 * np.einsum("mj,mj->m", a, b)
               + np.einsum("mj,mj->m", aa, A.reshape(-1, k * k))) / c) * 100
    quad = aa @ B.U.reshape(n, k * k).T                          # a.U_i a
    own = (1 - (B.e[None] - 2 * a @ B.u.T + quad) / B.e[None]) * 100
    vr[bad] = -np.inf
    own[bad] = -np.inf
    own[~masks] = np.nan
    return a, vr, own


def condition_number(B: Blocks, idx) -> float:
    """Condition number of the weighted design matrix W^1/2 G of the set,
    sqrt(cond(A_S)) (Bayesian ISOLA's resolvability measure; the value
    depends on how the basis is scaled, so a threshold must be calibrated
    in this basis)."""
    A = B.A[np.asarray(idx, dtype=int)].sum(0)
    return float(np.sqrt(np.linalg.cond(A)))
