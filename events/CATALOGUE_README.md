# catalogue.csv — column reference

One row per automated moment tensor solution. The CSV is regenerated from
the per-event `solution.json` archives after every processed event — it is
derived data, never hand-edited; `events/<dir>/solution.json` is the full
record (stations used/dropped with reasons, depth-search table, band
search, per-station zcor/VR, provenance) if you need more than a row.

All solutions are PRELIMINARY, deviatoric-only (no isotropic component)
and produced without human review; do not interpret mechanisms in
volcanic/geothermal settings. See docs/METHOD.md for the method and
README.md for citations.

## Origin (from GeoNet)

| Column | Meaning |
|---|---|
| `PublicID` | GeoNet event identifier (links to geonet.org.nz) |
| `Date` | Origin time, UTC (GeoNet) |
| `Latitude`, `Longitude` | Epicentre, degrees (GeoNet) |
| `GeoNet_M` | GeoNet preliminary magnitude (mixed types: M/MLv/mB) |
| `GeoNet_depth` | GeoNet initial depth, km (5/12/33 are fixed placeholders) |
| `GeoNet_depth_unc` | GeoNet depth uncertainty from QuakeML, km (blank if unavailable). For located depths our depth search is bounded to clamp(2x this, 5-15 km) around GeoNet_depth; placeholder depths search the full grid |

## Our solution

| Column | Meaning |
|---|---|
| `strike1,dip1,rake1` / `strike2,dip2,rake2` | Nodal planes of the DC part, degrees, Aki & Richards convention; the data cannot distinguish which plane is the fault |
| `Mw` | Moment magnitude from the inversion |
| `Depth` | Our centroid depth, km — the FINAL pick (VR-first, %DC tie-break on the contiguous plateau). The full depth grid is always searched; nothing is constrained toward GeoNet |
| `Depth_lo`, `Depth_hi` | Range of depths whose VR is within 10% of the maximum — the depth uncertainty (Vallée et al. 2011) |
| `Depth_VRmax` | Depth of the maximum variance reduction, km |
| `dZ_GeoNet` | Our centroid depth minus the GeoNet hypocentre depth, km (the search is bounded to ± 30 km around it, never forced) |
| `Depth_DCmax` | Depth of the maximum %DC, km |
| `Plateau_km` | Width of the near-VR-max plateau, km; large values mean depth is weakly constrained by the waveforms |
| `Mo` | Scalar moment, dyne-cm |
| `NS` | Number of stations used in the final solution (sets the row of the grade table) |
| `AzGap` | Largest azimuthal gap between used stations, degrees (informational — NOT a grade threshold; the geometry requirement is simply that two stations are >= 90 deg apart) |
| `Grade` | INGV table on (VR, N stations): the VR bar falls as N rises (3 stations: B ≥ 70; 4: B ≥ 40, A ≥ 60; 5-8: B ≥ 40, A ≥ 60; >8: B ≥ 30, A ≥ 50); A/B also need DC ≥ 60 and jackknife rotation ≤ 25°. X = no coherent solution. Only A/B are emailed |
| `DC`, `CLVD` | Percent double-couple / compensated linear vector dipole of the deviatoric solution |
| `MinStaVR` | Own VR of the worst-fitting used station |
| `VR` | Total variance reduction, percent (distance-weighted) |

## Stability (leave-one-station-out jackknife at the preferred depth)

| Column | Meaning |
|---|---|
| `Jk_n` | Number of jackknife subsets (blank if < 4 stations) |
| `Jk_Mw_std` | Standard deviation of Mw across subsets |
| `Jk_DC_std` | Standard deviation of %DC across subsets |
| `Jk_rot_deg` | Maximum DC-tensor rotation of any subset vs the full solution, degrees (plane-flip independent) |

## Displacement forecast (the purpose of this tool)

| Column | Meaning |
|---|---|
| `PredDisp_cm` | Peak predicted surface displacement, cm — Okada (1992) uniform-slip rectangle, Wells & Coppersmith (1994) dimensions, maximum over both nodal planes |
| `Detectable` | True when `PredDisp_cm` >= 1 cm (InSAR-interesting) |

## Moment tensor elements

`Mxx, Mxy, Mxz, Myy, Myz, Mzz` — deviatoric tensor in the mttime XYZ
basis, units of **1e20 dyne-cm** (matching the published NZ regional CMT
solutions CSV so the catalogues are directly comparable).

## Processing / bookkeeping

| Column | Meaning |
|---|---|
| `Band` | Chosen filter passband (e.g. `20-50s`); candidates tried per magnitude, best kept (VR first, %DC tie-break) |
| `Model` | 1-D velocity model used (Ristau 2008 North/South Island) |
| `Status` | `solved`, or `no_coherent_solution` for an event the network could not constrain (mechanism columns are then empty and the grade is `X`) |
| `Selection` | Station-selection/grading version that produced the row (e.g. `v4`). A catalogue mixing vintages is self-describing |
| `Code` | Short git commit of the code that produced the row, so any solution can be reproduced exactly |
| `quality_flag` | `True`, or the failed checks (`VR_below_grade_B`, `low_DC`, `unstable_mechanism`, `few_stations`, or `no_coherent_solution:<stage>`) |
| `publish_flag` | `True`, or why the email gate declined (`grade_C`, `too_small_no_disp`, `aftershock`, `daily_cap`) |
| `published` | Whether this solution was emailed to the list |

## Companion tables

- `not_published.csv` — every processed event that did NOT email, with
  `Why_not` (the publish decision's reasons) and `Tags`.
- `station_ledger.csv` — one row per station per event: distance,
  azimuth, sector, SNR per component, amplitude ratio, window, whether it
  was used, the reason class and reason if not, its own VR and time
  shift, and the event's band/depth/Mw/VR/grade. The year-one learning
  table; `station_performance.csv` is its per-station aggregate.
