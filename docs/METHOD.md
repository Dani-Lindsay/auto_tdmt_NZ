# auto_tdmt_NZ — Method and Implementation

Documentation of the automated regional moment tensor processing flow, for
review. Danielle Lindsay (danielle.lindsay@earthsciences.nz), September 2026.

Companion to the [README](../README.md) (setup/operations). This document
covers the science choices, module by module, and worked examples from the
2026-09-02 Fiordland sequence.

## 1. Motivation and scope

This is a personal, external project. Its purpose is rapid InSAR response:
within minutes of a New Zealand earthquake, decide whether it has likely
produced measurable surface displacement and whether action is needed
(acquisition tasking, processing preparation, checking the next NISAR
pass). A moment tensor is the prerequisite for that displacement forecast,
so the pipeline produces PRELIMINARY automated solutions
minutes-to-hours after an event using the Dreger-style time-domain method
with 1-D Green's functions and the published Ristau (2008) NZ velocity
models; solutions are therefore directly comparable to the published NZ
regional CMT solutions
([GeoNet/data moment-tensor](https://github.com/GeoNet/data/tree/main/moment-tensor)),
which serve as the validation reference. The MT solutions are a useful
by-product for other scientists; the displacement field is the point.
The project makes no representation about any organisation's internal or
operational systems.

The inversion engine is **mttime** (Chiang, LLNL;
github.com/LLNL/mttime) — the Python implementation of Dreger's TDMT used
operationally in California — driven end-to-end by this repository. We do
not reimplement any inversion mathematics.

## 2. Processing flow

```
GeoNet quake API poll (cron, 10 min) ........................ run01_watch.py
  processing floor: prelim M >= 3.7, NZ bbox, not "deleted"     trigger.py
     |
event processing ............................................ run02_process.py
  1. event metadata (quake API)                                 geonet.py
  2. station inventory + waveforms (FDSN NRT/archive)           waveforms.py
  3. pre-processing -> SAC displacement (cm), ZRT, 1 sps        waveforms.py
  4. Green's functions from precomputed CPS library             greens.py
  5. mttime deviatoric inversion + depth search                 invert.py
     station selection by FIT (the funnel, 3.1), grid-edge
     depth guard, ordered filter-band preference
  6. quality gates -> letter grade, or "no coherent solution"   invert.py
  7. Okada forward model of predicted surface displacement      okada_forward.py
  8. NISAR pass timing at the epicentre (NASA CMR)              nisar_dates.py
  9. figures (mttime waveform fits + cartopy maps)              figure.py
 10. publication gate -> email                                  trigger.py, publish.py
```

Every solution is archived (`events/<publicID>/`) whether or not it is
published, with full provenance in `solution.json`.

## 3. Data selection and pre-processing (`waveforms.py`)

### 3.1 Station selection (task 2, `waveforms.py`)

Selection v5 (2026-09) replaced two earlier designs — a pre-filter
cascade (v3) that starved itself, and a four-pass funnel (v4) that was
correct but complex — with rules taken from published operational
systems (the review is in `docs/lit_review/`). Nothing here is an
invention; every rule is one line in `auto_tdmt.cfg` §2 with its source.

1. **Distance window** scaled by magnitude: 120 / 180 / 250 / 300 km at
   M < 4 / < 4.5 / < 5 / ≥ 5, near limit 10 km below M4.5 and 20 km
   above. Gisola uses 10–250 km for M4–5 and 40–300 km for M5–5.5
   (Triantafyllis et al. 2022). One radius extension of 100 km when
   fewer than four stations are usable (offshore events).
2. **Unusable data** rejected: no waveform, no response, gaps, wrong
   sample count.
3. **Signal quality**: per-component SNR = RMS of the 200 s after the P
   arrival over RMS of the 200 s before it, measured on the filtered
   traces; a station is usable when the median of its three components
   is ≥ 2 (BMKG: Halauwet et al. 2024, GJI 239; AutoBATS uses 2.0, INGV
   5). Every published system that defines an SNR screens per
   component, not per station.
4. **Broken responses**: reject a station whose peak × distance exceeds
   three times the network median. W-phase screens 0.1×–3× on both
   sides (Duputel et al. 2012); the low side is deliberately not
   screened because a station near a nodal plane genuinely has small
   amplitude, and Duputel et al. name exactly that failure mode.
5. **Azimuth balance**: rank the usable stations by SNR and take them
   round-robin across eight 45° sectors until 16 (Gisola's sector
   design, SCARDEC's best-SNR-per-bin; relaxed from Gisola's hard cap of
   two per sector so a one-sided offshore geometry fills up from its
   best stations rather than being starved).
6. **Filter bands by magnitude**, an ordered preference: 10–50 s below
   M4.5; 10–50 s then 20–50 s to M5.5; 20–100 s then 30–100 s above.
   The first band whose solution passes its gate wins, because VR is
   not comparable across bands (a longer period is smoother and scores
   higher even when it fits noise). Near-identical to SCSN's menu
   (Clinton, Hauksson & Solanki 2006) and INGV's practice of pushing
   small events to higher frequencies.
7. **Record window** per station: 30 s before origin to distance /
   2.5 km/s + a magnitude-dependent tail after it, extended to where the
   smoothed envelope decays back toward the pre-event level (slow
   Hikurangi paths), capped at 200 s after origin. The earlier 120 s cap
   cut 300 km stations mid-train.

Everything that enters and leaves is recorded with a reason string in a
shared vocabulary (`invert.reason_class`), so the all-station waveform
figure and `events/station_ledger.csv` show why.

### 3.2 Pre-processing

- **Far-field guard**: a station inside 3x the source depth is TAGGED
  `near_field`, not excluded (§3.1). The Green's functions include the
  near-field terms, so the fit is better evidence than the heuristic.
- **Waveform windows**: origin−150 s to origin+230 s downloaded; final
  cut origin−30 s to origin+200 s.
- **Inversion record length**: distance-adaptive at every magnitude —
  each station'''s fitted window is 30 s pre-origin + distance/2.5 km/s +
  a magnitude-dependent tail (30 s below M4.5, 40 s to M5.5, 60 s
  above), clamped to 60–150 s. A close station'''s train is over quickly
  regardless of event size; fitting the empty tail only taxes VR. This
  is a deterministic kinematic cut, not an amplitude-based one, for
  reproducibility.
- **Response removal**: to displacement with pre-filter
  (0.004, 0.007, 10, 20) Hz, then rotation to ZNE and NE→RT along the
  great-circle back-azimuth.
- **Peak-to-noise measurement**: per station, median over the three
  components of peak|signal| / RMS(pre-event noise), signal window = the
  distance-adaptive inversion window above, noise window −115 to −2 s.
  This replaced an RMS-ratio SNR gate (threshold after Ristau 2008)
  that could not order stations the way visual inspection does: an
  emergent band-limited packet is a localised peak, and RMS over a
  200 s record dilutes it with empty tail. Dropped stations are
  recorded with their peak/noise in `solution.json`.
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
  1.0–5.0 km at 0.5 km (fine where shallow depth discrimination
  happens), 6–10 km at 1 km, 12–30 km at 2 km, and 34–58 km at 4 km
  (Fiordland subduction events exceed crustal depths); dt = 1 s,
  npts = 256, vred = 0 (traces start at origin time). Runtime lookup
  takes the nearest grid distance (≤2.5 km error, absorbed by the
  per-station time-shift search).
- Libraries are stored unfiltered and band-passed per event to match the
  data exactly. A manifest records the model file's SHA-256, grid, and
  build time. CI never runs CPS — it downloads the library tarball.

## 5. Inversion (task 3, `invert.py`)

- **Engine**: mttime deviatoric inversion (degree 5), ZRT components,
  distance weighting (F-net's published VR definition weights stations
  proportional to distance), per-station cross-correlation time shifts.
  mttime's own control file is written and read; nothing of the
  inversion mathematics is reimplemented.
- **Depth grid**: the GeoNet hypocentre ± 30 km, clipped to the library;
  the full library when GeoNet's depth is a placeholder. The search is
  bounded, the answer is never forced — the same design as F-net
  (JMA ± 30 km), Gisola (± 31 km), W-phase and SCARDEC (± 50 km). An
  unbounded grid let VR climb monotonically into the smoothest Green's
  functions (2026p666955: 50 km against a 26 km hypocentre).
- **Preferred depth**: the maximum VR, mttime's own pick and every
  operational system's. Recomputing the pick over 343 archived events
  with a Ristau reference depth showed a %DC tie-break adds nothing
  (median |ΔZ| 8.0 km either way). **Depth uncertainty** is the range of
  depths whose VR is within 10% of the maximum (SCARDEC: Vallée et al.
  2011; Bernardi et al. 2004), reported as `Depth_lo`/`Depth_hi`.
- **Station selection by fit — the Clinton loop** (Clinton, Hauksson &
  Solanki 2006, BSSA 96): invert every selected station; drop any whose
  own VR is below max(overall VR − 10, 25) **or whose solved time shift
  exceeds 8 s**; re-invert; repeat until nothing drops or three stations
  remain (at most five rounds). mttime's shift search is unbounded, so
  the cap is what stops a noise trace sliding into a chance alignment
  (NEIC SynDepth clips at 6/8/10 s by magnitude). No other selection
  logic exists.
- **Jackknife**: leave-one-station-out at the preferred depth; the
  largest minimum-rotation angle (Townend et al. 2012) of any subset
  from the full solution is the stability evidence. NIED state in print
  that the jackknife is the only reliable detector of a station gone
  bad, because the misfit does not localise on the offending station
  (Fukuyama et al. 1998).

## 6. Quality grades and publication (`invert.py`, `trigger.py`)

The grade is INGV's table (terremoti.ingv.it/en/help#TDMT), in which
the VR needed for a grade FALLS as the station count rises:

| N stations | D | C | B | A |
|---|---|---|---|---|
| 3 | VR < 20 | 20–70 | ≥ 70 | — |
| 4 | < 20 | 20–40 | 40–60 | ≥ 60 |
| 5–8 | < 15 | 15–40 | 40–60 | ≥ 60 |
| > 8 | < 15 | 15–30 | 30–50 | ≥ 50 |

A and B additionally require %DC ≥ 60 (the BSL publishability rule;
INGV's `a` suffix) and a jackknife rotation ≤ 25° when the jackknife is
possible. The shape of the table is the point: because VR is
anti-correlated with station count, a flat "VR > X" gate is gameable by
the selection itself — Triantafyllis et al. (2016) show ten stations at
VR 0.7 with a condition number of 5 becoming two stations at VR 0.9
with a condition number above 10. Under this table grade A is
unreachable with three stations by construction. F-net's flat rule
(VR > 50%, M > 3.5) is the simpler alternative; the table agrees with
it for typical station counts and is stricter for small ones.

**No coherent solution.** When fewer than three stations survive the
loop the event is archived with `"status": "no_coherent_solution"` and
grade `X`: the full station ledger and every round's numbers, but no
mechanism, magnitude or depth. F-net simply does not publish below its
floor; GeoNet fall back to USGS W-phase for events they cannot
constrain. Warnings (recorded, never blocking): the preferred depth
touching the ± 30 km window edge; the jackknife being impossible.

Publication gate, applied to **our inverted Mw**: grade A or B and
(Mw ≥ 5.0 or the Okada-predicted peak displacement ≥ 1 cm). Anti-spam:
at most three emails a day; within 75 km and 14 days of a published
event, a smaller event must be within 0.5 Mw of it to publish.
Everything processed is archived regardless, and every event that did
not publish is listed with its reason in `events/not_published.csv`.

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
- Validation against the published NZ CMT solutions is implemented as a
  comparison hook but awaits their solutions for overlapping events
  (their CSV updates ~monthly; note their MT units are 1e20 dyne-cm).

### Verification aids

Nodal-plane and fault-geometry conventions were cross-checked numerically
(conjugate-plane far-field equivalence, rake/side acid tests — now anchor
tests in the suite) and visually against the USGS finite-fault event pages
(black-outlined plane, red up-dip edge convention, adopted here) and the
interactive Focal Mechanism Explorer at https://eq.comoglu.com/bb/.

## Validation metric

Mechanism agreement with reference catalogues (`run05_validate.py`) and
jackknife stability are measured as the **minimum rotation angle**
between two double couples (J. Townend, pers. comm. 2026-08-20):
each mechanism is expressed as a rotation matrix of its principal axes
with respect to geographic coordinates (Walsh et al. 2009, eqs 1-3),
and the angle is arccos((tr(R1^T R2) - 1)/2), minimised over the
double-couple symmetry group so the result is independent of which
nodal plane parameterises either mechanism (cf. Kagan 1991). The
approach follows Townend et al. (2012), supplement eq. 1
(`docs/Townend_etal_2012_supplement.pdf`).

## References

- Walsh, D., Arnold, R., Townend, J. (2009). A Bayesian approach to
  determining and parameterising earthquake focal mechanisms. GJI 176,
  235-255.
- Townend, J., et al. (2012). Three-dimensional variations in present-day
  tectonic stress along the Australia-Pacific plate boundary in New
  Zealand. EPSL, doi:10.1016/j.epsl.2012.08.003 (supplement in docs/).
- Kagan, Y. Y. (1991). 3-D rotation of double-couple earthquake sources.
  GJI 106, 709-716.
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
