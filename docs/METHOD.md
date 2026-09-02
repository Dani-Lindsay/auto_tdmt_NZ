# auto_tdmt_NZ — Method and Implementation

Documentation of the automated regional moment tensor processing flow, for
review. Danielle Lindsay (danielle.lindsay@earthsciences.nz), September 2026.

Companion to the [README](../README.md) (setup/operations). This document
covers the science choices, module by module, and worked examples from the
2026-09-02 Fiordland sequence.

## 1. Motivation and scope

GeoNet computes regional moment tensors manually and publishes them in a
CSV updated roughly monthly ([GeoNet/data moment-tensor](https://github.com/GeoNet/data/tree/main/moment-tensor));
there is no automated or near-real-time MT product for New Zealand. This
pipeline produces PRELIMINARY automated solutions minutes-to-hours after an
event, using the same class of method (Dreger-style time-domain MT with 1-D
Green's functions) and the same velocity models (Ristau 2008) as GeoNet's
own regional CMT analysis, so solutions are directly comparable to the
manual catalogue when it appears.

The inversion engine is **mttime** (Chiang, LLNL;
github.com/LLNL/mttime) — the Python implementation of Dreger's TDMT used
operationally in California — driven end-to-end by this repository. We do
not reimplement any inversion mathematics.

## 2. Processing flow

```
GeoNet quake API poll (cron, 10 min) ........................ run01_watch.py
  processing floor: prelim M >= 4.0, NZ bbox, not "deleted"     trigger.py
     |
event processing ............................................ run02_process.py
  1. event metadata (quake API)                                 geonet.py
  2. station inventory + waveforms (FDSN NRT/archive)           waveforms.py
  3. pre-processing -> SAC displacement (cm), ZRT, 1 sps        waveforms.py
  4. Green's functions from precomputed CPS library             greens.py
  5. mttime deviatoric inversion + depth search                 invert.py
     for each filter band in a magnitude-dependent menu;
     low-VR station rejection pass; VR+%DC preferred pick
  6. quality gates                                              invert.py
  7. Okada forward model of predicted surface displacement      okada_forward.py
  8. NISAR pass timing at the epicentre (NASA CMR)              nisar_dates.py
  9. figures (mttime waveform fits + cartopy maps)              figure.py
 10. publication gate -> email                                  trigger.py, publish.py
```

Every solution is archived (`events/<publicID>/`) whether or not it is
published, with full provenance in `solution.json`.

## 3. Data selection and pre-processing (`waveforms.py`)

- **Stations**: network NZ broadbands (HH? preferred over BH?; location
  codes taken from the station inventory, never assumed), 20–400 km
  epicentral distance, nearest 12 that pass QC. Station distance must
  exceed 3x source depth (point-source/far-field assumption; waived for
  GeoNet placeholder depths, which are unreliable).
- **Waveform windows**: origin−150 s to origin+230 s downloaded; final cut
  origin−30 s to origin+200 s.
- **Response removal**: to displacement with pre-filter
  (0.004, 0.007, 10, 20) Hz, then rotation to ZNE and NE→RT along the
  great-circle back-azimuth.
- **SNR gate**: per station, min over the three components of
  RMS(signal: origin→+200 s) / RMS(noise: −120→−10 s), measured in the
  inversion passband. Threshold 2.0, after Ristau (2008): "a SNR higher
  than 2 is normally required to calculate a reliable moment tensor."
  Dropped stations are recorded with their SNR in `solution.json`.
- **Filtering**: zero-phase 3-corner Butterworth, band from the menu in
  §5; 5% taper; decimate/resample to 1 sps; convert m→cm (TDMT
  convention). Green's functions receive the identical filter.
- Data are never padded, interpolated across gaps, or substituted: any
  deficiency drops the station with a recorded reason (fail loud).

The chain follows the mttime example notebooks (Chiang) and matches the
preparation described in Ristau (2008) §"Preparation of the observed
waveforms".

## 4. Green's functions (`greens.py`)

- **Velocity models** (`models/*.d`, model96 format): Ristau (2008), SRL
  79(3), Table 1 — North Island (8 layers) and South Island (5 layers,
  Moho 39 km) models with Qp=400/Qs=200, transcribed verbatim from the
  paper (doi:10.1785/gssrl.79.3.400). Model chosen per event by a
  North/South Island split at Cook Strait (`config.model_for_event`);
- **Computation**: CPS 3.30 (Herrmann) wavenumber integration, run once
  locally per model:
  `hprep96 -EQEX | hspec96 | hpulse96 -D -i | f96tosac -B`,
  producing the ten Herrmann fundamental-source time series
  (ZDD RDD ZDS RDS TDS ZSS RSS TSS ZEX REX) that mttime's
  `green="herrmann"` mode expects (displacement, cm, 1e20 dyne-cm source).
- **Library grid**: distances 10–500 km at 5 km spacing; source depths
  2–30 km at 2 km and 34–58 km at 4 km (Fiordland subduction events
  exceed crustal depths); dt = 1 s, npts = 256, vred = 0 (traces start at
  origin time). Runtime lookup takes the nearest grid distance
  (≤2.5 km error, absorbed by the per-station time-shift search).
- Libraries are stored unfiltered and band-passed per event to match the
  data exactly. A manifest records the model file's SHA-256, grid, and
  build time. CI never runs CPS — it downloads the library tarball.

## 5. Inversion (`invert.py`)

- **Engine**: mttime deviatoric inversion (degree 5), ZRT components,
  inverse-distance weighting, per-station cross-correlation time shifts
  (zcor; the pandas<3 pin is required for this — mttime issue #15).
- **Depth search**: all library depths when GeoNet's depth is a fixed
  placeholder (5/12/33 km); otherwise library depths within ±20 km of the
  hypocentre. Inversion window npts=150 at 1 sps, ts=30.
- **Filter-band menu** (BSL practice: a small menu of period bands, longer
  periods for larger events; the pipeline tries each and picks by the rule
  below):
  | preliminary M | candidate bands |
  |---|---|
  | < 4.5 | 10–50 s, 20–50 s |
  | 4.5–5.5 | 20–50 s, 10–50 s, 20–100 s |
  | ≥ 5.5 | 20–100 s, 30–100 s |
- **Station rejection**: after a first inversion, stations with individual
  VR < 10% at the preferred depth are removed and the inversion is rerun
  once (cf. dropping persistently low-VR stations, Ristau 2008; Dreger &
  Helmberger 1993). This is the counterweight to the permissive SNR gate:
  stations that passed SNR 2.0 but do not actually fit are removed here
  rather than being allowed to dilute the %DC. Rejections are recorded
  with reasons.
- **Preferred solution rule**: among solutions (depths; and then bands)
  whose total VR is within 5 percentage points of the maximum, take the
  highest %DC. Rationale: VR is a weak discriminator with depth (it often
  climbs monotonically) while spurious CLVD grows where the depth/model is
  wrong; pure VR-max picks produced solutions with implausible ~70% CLVD
  that the %DC-aware rule resolves (see worked examples). Both the VR-max
  and %DC-max depths are recorded and a disagreement flag is set.

## 6. Quality gates and publication (`invert.py`, `trigger.py`)

Solution quality gates (must all pass for publication):
- ≥ 3 stations used;
- total VR ≥ 50%;
- azimuthal gap ≤ 270°.

Recorded as warnings (do not block): preferred depth at the search-grid
edge (a shallow crustal event legitimately prefers the shallowest depth).

Publication gate, applied to **our inverted Mw** (never GeoNet's mixed
ML-type preliminary magnitudes): publish if Mw ≥ 5.0 OR the Okada-predicted
peak surface displacement ≥ 1 cm (the InSAR-interesting case). Anti-spam:
max 3 emails/day; within 75 km and 14 days of an already-published event,
a smaller event must be within 0.5 Mw of it (or exceed the Mw gate) to
publish. Everything processed is archived regardless.

## 7. Deformation forward model (`okada_forward.py`)

For each nodal plane (the MT cannot distinguish them): rectangular uniform-
slip Okada (1992) dislocation via okada4py; length/width from Wells &
Coppersmith (1994, all-type regressions); slip from M0 = μLW·s with
μ = 30 GPa; centroid at the inverted depth (plane pushed down if its top
would breach the surface); displacement evaluated on a 121×121 km grid.
Calibration anchors: Mw 5.5 at 10 km ≈ 1.6 cm peak |u|; Mw 4.0 at 5 km ≈
0.4 mm (hence the displacement gate passes only shallow M≳4.7 events).

## 8. Outputs and figures (`figure.py`, `diagnostics.py`)

Per event: `solution.json` (all numbers + provenance), `draft_email.txt`,
mttime waveform-fit and depth-search figures, and three outward figures for
the email: (1) mttime's waveform-fit page, untouched (full Deviatoric =
DC + CLVD decomposition, per-station fits, VR/%DC); (2) the modelling
figure — Mercator map panels of the station geometry with the moment
tensor beachball drawn from the actual tensor elements (obspy/mopad
full-MT rendering, so any CLVD is represented honestly rather than
collapsed to the closest DC) and the Okada-predicted E/N/U displacement
over the modelled area (cmcrameri 'vik', modelled nodal plane outlined),
with the NISAR pass table and provenance line; (3) a depth-sensitivity
summary (VR, %DC/%CLVD, Mw vs depth) that also displays the
selection-rule window. All plotting is matplotlib (one stack; maps
via cartopy with ocean/land/coastline features and gridlines).

`--debug` additionally saves stage-by-stage QC figures per band: raw
counts, response-removed displacement, filtered ZRT record sections, and a
station-geometry map.

Every `solution.json` records: velocity model + GF version, mttime/obspy
versions, filter band, all depth-search rows, stations used (with SNR,
distance, azimuth) and dropped (with reasons), the band-search table,
quality gates/warnings, forward-model parameters, NISAR passes, and the
publication decision with reasons.

## 9. Worked examples: Fiordland sequence, 2026-09-02 (Milford Sound)

Four events at the southern termination of the Alpine Fault / Fiordland
subduction corner, all with GeoNet fixed 5 km placeholder depths. Processed
with the South Island model. See `docs/examples/` for the figures; numbers
below are from the final run (SNR ≥ 2.0).

*(Table completed from the archived solution.json files — see
docs/examples/README.md.)*

Key observations:
1. **Velocity model dominates solution quality.** With GIL7 (California)
   the M5.6 mainshock inverted to 25% DC; with the Ristau South Island
   model, 81% DC at the same VR — the spurious CLVD was model error.
2. **Depth recovery from placeholders.** Both M5.5+ events prefer ~8 km
   over the 5 km fixed depth, with clean single-peaked VR(depth) curves.
3. **ML–Mw offset.** GeoNet preliminary magnitudes exceed our Mw by
   0.4–0.6 units across the sequence — the recurring NZ ML–Mw discrepancy
   discussed by Ristau (2008) — which is why all gating uses our Mw.
4. **Gates catch the marginal cases.** The M4.8 aftershock lost all but
   2 stations to the SNR gate in its first-choice band and was blocked by
   the 3-station minimum rather than emailing a poorly-constrained
   solution.

## 10. Known limitations / review questions

- 1-D models; no offshore-specific model yet (Ristau 2008 anticipates
  offshore-region issues north of North Island and at Puysegur).
- North Island model transcribed but not yet exercised on a real event.
- The North/South model split is a crude Cook Strait rule.
- Depth grid starts at 2 km; very shallow events sit at the grid edge
  (flagged, not blocked).
- No isotropic term (degree 5); volcanic/geothermal events with real ISO
  components will be forced deviatoric.
- Single-pass station rejection; no jackknife uncertainty yet (the EPS207
  jackknife scheme is the natural next addition).
- Validation against GeoNet's manual CMT catalogue is implemented as a
  comparison hook but awaits their solutions for overlapping events
  (their CSV updates ~monthly; note their MT units are 1e20 dyne-cm).

## References

- Chiang, A. MTtime: Time Domain Moment Tensor Inversion in Python.
  LLNL-CODE-814839. github.com/LLNL/mttime
- Dreger, D., & Helmberger, D. (1993). Determination of source parameters
  at regional distances with three-component sparse network data. JGR 98.
- Herrmann, R. B. (2013). Computer Programs in Seismology. SRL 84.
- Okada, Y. (1992). Internal deformation due to shear and tensile faults
  in a half-space. BSSA 82.
- Ristau, J. (2008). Implementation of routine regional moment tensor
  analysis in New Zealand. SRL 79(3), 400–415. doi:10.1785/gssrl.79.3.400
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships
  among magnitude, rupture length, rupture width, rupture area, and
  surface displacement. BSSA 84.
