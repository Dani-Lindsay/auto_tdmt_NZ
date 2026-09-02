"""Forward-model predicted surface displacement from a moment tensor
solution using a rectangular Okada (1992) dislocation via okada4py.

For each nodal plane (the MT cannot distinguish them):
  - fault dimensions from Wells & Coppersmith (1994) "all" regressions
      log10(RLD) = -2.44 + 0.59 Mw   (subsurface rupture length, km)
      log10(RW)  = -1.01 + 0.32 Mw   (downdip width, km)
  - uniform slip from M0 = mu * L * W * s
  - fault plane centred on the centroid (epicentre lon/lat, inverted depth)
  - displacement on a local east/north grid around the epicentre

Outputs peak |u|, peak vertical and the full grids for plotting.
"""

from __future__ import annotations

import numpy as np

import config


def wells_coppersmith_lw(mw: float) -> tuple[float, float]:
    """(length_m, width_m) from Wells & Coppersmith (1994), all slip types."""
    assert 3.0 < mw < 9.0, f"Mw {mw} outside regression range"
    length_km = 10 ** (-2.44 + 0.59 * mw)
    width_km = 10 ** (-1.01 + 0.32 * mw)
    return length_km * 1e3, width_km * 1e3


def m0_dyne_cm_to_nm(m0: float) -> float:
    return m0 * 1e-7


def predicted_displacement(
    mw: float, m0_dyne_cm: float, depth_km: float,
    strike: float, dip: float, rake: float,
    halfwidth_km: float | None = None, step_km: float | None = None,
) -> dict:
    """Surface displacement grids (east, north, up in m) for one nodal plane.

    Grid is local: x east, y north, in km relative to the epicentre.
    """
    import okada4py

    halfwidth_km = halfwidth_km or config.FORWARD_GRID_HALFWIDTH_KM
    step_km = step_km or config.FORWARD_GRID_STEP_KM

    length_m, width_m = wells_coppersmith_lw(mw)
    m0_nm = m0_dyne_cm_to_nm(m0_dyne_cm)
    slip_m = m0_nm / (config.SHEAR_MODULUS_PA * length_m * width_m)

    # keep the top of the fault at or below the surface
    depth_m = depth_km * 1e3
    dip_rad = np.radians(dip)
    top_m = depth_m - 0.5 * width_m * np.sin(dip_rad)
    if top_m < 0:
        depth_m -= top_m  # push centroid down so the plane stays buried
    assert depth_m > 0, "fault centre above surface"

    ss = slip_m * np.cos(np.radians(rake))
    ds = slip_m * np.sin(np.radians(rake))

    n = int(round(2 * halfwidth_km / step_km)) + 1
    axis_m = np.linspace(-halfwidth_km, halfwidth_km, n) * 1e3
    gx, gy = np.meshgrid(axis_m, axis_m)
    xs, ys = gx.ravel(), gy.ravel()
    zs = np.zeros_like(xs)

    u, d, s, flag, flag2 = okada4py.okada92(
        xs, ys, zs,
        np.array([0.0]), np.array([0.0]), np.array([depth_m]),
        np.array([length_m]), np.array([width_m]),
        np.array([float(dip)]), np.array([float(strike)]),
        np.array([ss]), np.array([ds]), np.array([0.0]),
        config.SHEAR_MODULUS_PA, config.POISSON_NU,
    )
    assert (np.asarray(flag) == 0).all(), "okada4py flagged bad geometry"
    u = np.asarray(u).reshape(len(xs), 3)
    ue = u[:, 0].reshape(gx.shape)
    un = u[:, 1].reshape(gx.shape)
    uz = u[:, 2].reshape(gx.shape)

    mag = np.sqrt(ue**2 + un**2 + uz**2)
    return {
        "x_km": axis_m / 1e3,
        "y_km": axis_m / 1e3,
        "ue_m": ue,
        "un_m": un,
        "uz_m": uz,
        "peak_abs_m": float(mag.max()),
        "peak_uz_m": float(np.abs(uz).max()),
        "fault": {
            "strike": strike, "dip": dip, "rake": rake,
            "length_km": length_m / 1e3, "width_km": width_m / 1e3,
            "slip_m": slip_m, "centroid_depth_km": depth_m / 1e3,
            "scaling": "Wells & Coppersmith 1994 (all types)",
        },
    }


def forward_both_planes(solution: dict) -> dict:
    """Run the forward model for both nodal planes of a solution.json dict."""
    pref = solution["preferred"]
    out = {}
    for name in ("plane1", "plane2"):
        p = pref[name]
        out[name] = predicted_displacement(
            mw=pref["mw"], m0_dyne_cm=pref["m0_dyne_cm"],
            depth_km=pref["depth_km"],
            strike=p["strike"], dip=p["dip"], rake=p["rake"],
        )
    out["peak_abs_m"] = max(out["plane1"]["peak_abs_m"],
                            out["plane2"]["peak_abs_m"])
    out["detectable"] = out["peak_abs_m"] >= config.PUBLISH_MIN_PRED_DISP_M
    return out
