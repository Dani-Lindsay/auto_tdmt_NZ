# USGS NEIC — operational MT practice and QC

Findings gathered 2026-09-05. Full-text obtained for Hayes, Rivera &
Kanamori (2009), Hayes (2011), Patton et al. (2016) and the operational
SynDepth source; Herrmann, Benz & Ammon (2011) is closed access and
**only its abstract was read** — its "codification of data quality
tests" is therefore NOT verified here.

## The numbers that matter to us

**Time shifts (c).** NEIC's SynDepth cross-correlates observed against
synthetic and then **hard-clips the shift to +/-6, 8 or 10 s** by
magnitude bin (`ShiftAmount = 6.0, 8.0, 10., 10., 10.` for
Mw < 5.5 / 5.5-6 / 6-6.5 / 6.5+). Source:
https://code.usgs.gov/ghsc/neic/algorithms/syndepth (software DOI
10.5066/P924LDLT), `asp.json`.
=> **Our 8 s absolute bound sits inside published operational practice.**

**Signal quality (b).** SynDepth's SNR is a **peak-amplitude ratio, per
channel, measured after demean/taper/bandpass**: max|signal| in a window
from (P - 10 s) forward, over max|noise| in the pre-signal trace, with
**magnitude-binned thresholds 2.0 / 2.25 / 2.5 / 3.0**
(`SNRCutoff` in `asp.json`; `checkSNR_PeakValue()` in
`syndepth/util/processData.py`). If nothing passes, the run returns
nothing.
=> **Our peak-to-noise metric is the same kind of measurement** (peak
amplitude over noise, per channel, post-filter). Our dead-channel floor
of 1.2 is more permissive than their 2.0-3.0; the difference is that we
then make every station earn a seat by fit, which they do not.

**Depth (e).** Depth is **never tied to the hypocentre** anywhere in
NEIC's MT stack: Mww centroid depths come off a 1 km grid (empirically
offset by 0.5 km, with 11.5 km a common shallow floor), the SLU/USGS
regional MT grid-searches ~1-29 km at 1 km spacing, and SynDepth
searches 4-750 km. The hypocentre supplies only the STARTING point.
**But** SynDepth does let the starting depth constrain the search
window: `CMTSecondCheck = 15`, `SecondMinimumDepth = 9` — "if start
depth is larger than CMTSecondCheck then the minimum depth is set to
SecondMinimumDepth **to prevent very shallow depths which can cause
local minimums**".
=> A precedent for using the hypocentre to bound the SEARCH, not the
answer — the opposite direction to our problem (we get depths too deep,
they were guarding against too shallow), but the same principle.

**Grid-edge depths (e).** SynDepth prints quality checks for every
solution: **number of local minima** in the fit surface, the prominence
of the chosen minimum on each side, and **"On Max Edge" / "On Min Edge"**
flags (`qualityChecks()` in `syndepth/util/plotData.py`). Hayes et al.
(2009) apply the same idea to the W-phase centroid grid: "This initial
grid size is increased (in both grid searches) **if the solution is
within one cell of the grid edge**."
=> **Our grid-edge guard is independently standard practice.** Their
response (widen the grid) differs from ours (prefer the best interior
local maximum); worth considering.

**Mislocation from time shifts (c) — WORTH ADOPTING.** The SLU/USGS
regional MT solution pages do not clip the shift; they **fit its
azimuthal pattern** to detect a location error:
"Time_shift = A + B cos(Azimuth) + C sin(Azimuth)", and report the
implied shift in origin time and epicentral coordinates. Per-trace time
shift AND per-trace variance reduction are published for every station
of every event. Example page (read in full):
https://www.eas.slu.edu/eqc/eqc_mt/MECH.NA/MECH.NA.2025/20250109035217/HTML.REG/index.html
=> We currently only clip. Fitting the pattern would tell us whether a
systematic offset is a MISLOCATION rather than bad stations — directly
relevant to our depth disagreements.

**Band and magnitude floor (b/d).** Herrmann, Benz & Ammon (2011)
abstract: "**Using the 0.02-0.10 Hz passband**, we can usually
determine ... moment tensor solutions for **earthquakes with Mw as
small as 3.7**", with the threshold "significantly influenced by the
density of stations, the location of the earthquake relative to the
seismic stations and, of course, the signal-to-noise ratio."
=> **Our 10-50 s band IS 0.02-0.10 Hz, and our processing floor is
prelim M 3.7.** Both match the published North American regional
practice exactly. This is the closest analogue to our pipeline.

**Record window (d).** SLU/USGS regional MT cuts
`o DIST/3.3 -40  o DIST/3.3 +40` — a **3.3 km/s group-velocity window,
-40 to +40 s**, with a 0.05-0.15 Hz filter for an M4.
=> Ours is 2.5 km/s plus a magnitude-dependent tail, envelope-extended.
Theirs is tighter and symmetric about the group arrival.

**Publication gating (f) — the important negative result.**
There is **no published numeric threshold anywhere in NEIC's MT stack**:
no minimum variance reduction, no minimum station count, no azimuthal
gap limit. What exists is (i) named but unquantified gates in Hayes et
al. 2009 (an SNR screen that can abort the inversion; station count;
the condition number of the solution), (ii) diagnostics printed for
analyst judgement, and (iii) a state machine — solutions publish
automatically as `preliminary` and are promoted to `reviewed` /
`confirmed` by a duty seismologist, every manual inclusion/exclusion
logged in Hydra's "Passport" (Patton et al. 2016).
=> **We cannot cite a threshold; we have to derive and defend our own.**
Our grade rubric is doing something NEIC does not publish.

**Empirical, from the products themselves** (analysis of 39 Mww
products Jan-Mar 2025, not a published claim): every Mww channel had
weight 1 and offset 0 — **W-phase applies no per-station time shift at
all** — while Mwb and Mwr products carry per-station offsets and
fractional weights. Minimum retained per-channel waveform fit across
those Mww solutions was 0.62, suggesting a post-hoc cutoff near 0.6.
Variance reduction 68-92%, channel counts 12-119.

## Key references

- Hayes, G. P., Rivera, L., & Kanamori, H. (2009). Source inversion of
  the W-phase: real-time implementation and extension to low magnitudes.
  SRL 80(5), 817-822. doi:10.1785/gssrl.80.5.817
- Herrmann, R. B., Benz, H., & Ammon, C. J. (2011). Monitoring the
  earthquake source process in North America. BSSA 101(6), 2609-2625.
  doi:10.1785/0120110095  [ABSTRACT ONLY - closed access]
- Yeck, W. L., Herrmann, R. B., Patton, J., Barnhart, W. D., & Benz,
  H. M. (2025). Estimating earthquake source depth using teleseismic
  broadband waveform modeling at the USGS NEIC. SRL 96(6), 3643-3655.
  doi:10.1785/0220240372  [abstract only; operational code public]
- Patton, J. M., et al. (2016). Hydra - the NEIC's 24/7 seismic
  monitoring ... tool suite. USGS OFR 2016-1128.
  doi:10.3133/ofr20161128
- Hayes, G. P. (2011). Rapid source characterization of the 2011 Mw 9.0
  Tohoku earthquake. EPS 63(7), 529-534. doi:10.5047/eps.2011.05.012
- SynDepth source + parameters:
  https://code.usgs.gov/ghsc/neic/algorithms/syndepth
  (software DOI 10.5066/P924LDLT, Green's functions DOI 10.5066/P996KCPK)
