# Station selection and quality control in operational moment-tensor systems

A review of how twelve operational regional (and, for context, two
teleseismic) moment-tensor systems select stations, measure signal
quality, bound time shifts, choose depth and gate publication — and
what this pipeline adopted from each. Reviewed September 2026.

Status markers used below: **ADOPTED** (implemented in the current
pipeline, selection v5), **ADAPTED** (implemented with a stated
modification), **DEFERRED** (published practice deliberately left for
a later version; listed in `docs/METHOD.md` §10), **NOT APPLICABLE**
(teleseismic or otherwise outside this pipeline's scope).

Evidence markers: **[verified]** read in the primary source, code or
configuration; **[secondary]** via a citing paper; **[not verified]**
could not be confirmed.

---

## 1. Berkeley / SCSN — the TDMT lineage this pipeline descends from

The method is Dreger & Helmberger (1993), *JGR* 98, 8107–8125,
doi:10.1029/93JB00023, and Dreger (2003), *TDMT_INV*, IASPEI Handbook
81B, p. 1627; the automated Berkeley system is Pasyanos, Dreger &
Romanowicz (1996), *BSSA* 86(5), 1255–1269 (full text not obtained).
Berkeley's own operational thresholds are not documented publicly
beyond "automatically produced solutions of high quality" being
published to the web [verified, BSL annual report].

The fully documented automation of the same code is **Clinton,
Hauksson & Solanki (2006)**, *An evaluation of the SCSN moment tensor
solutions*, *BSSA* 96(5), 1689–1705, doi:10.1785/0120050241
[verified]. Its rules:

- Pre-inversion screens: sensor corner period ≥ 100 s; epicentral
  distance 45–700 km; peak amplitude below 80% of the clip level. No
  SNR threshold.
- Geometry: six azimuthal sectors, subdivided until six are populated;
  the initial pick per sector is the station nearest 60 km.
- Filter bands by ML: < 4.2 → 10–50 s; 4.2–5.5 → 20–50 s; > 5.5 →
  20–100 s.
- Depth: fixed trial set 5, 8, 11, 15, 18, 21 km; maximum overall VR.
- The selection loop: gather data, select by azimuth at optimal
  distance, invert, reject the stations with the poorest individual
  fits, reselect to maximise azimuthal distribution, repeat; relax the
  required quality if the station list is exhausted.
- Grades on overall VR (OVR) and individual-station VR (IVR): **A**
  (published without review) six stations, OVR > 60%, reject IVR <
  max(OVR − 10, 50); **B** four stations, OVR > 40%, reject IVR <
  max(OVR − 15, 25), Mw distributed but the mechanism withheld; **C**
  not distributed.
- Time shifts vetted against a distance regression: ZCOR_expected =
  0.13 r + 1.42 s, station removed when |observed − expected| > 9 s.

**ADOPTED**: the selection loop; the filter-band menu (near-identical);
the quality-B own-VR floor of 25 as a fixed bar. **ADAPTED**: the
relative rule (IVR < OVR − 10) is off by default — without SCSN's
reselection step it ratchets, because every removal raises the overall
VR and therefore the bar (observed: a 16-station M4.9 reduced to 3
stations in three rounds); it remains available as `stationVRDrop`.
**DEFERRED**: the distance-regressed time-shift residual; the clipping
screen.

## 2. GeoNet / Ristau — the reference catalogue

Ristau (2008), *SRL* 79(3), 400–415, doi:10.1785/gssrl.79.3.400, and
Ristau (2013), *BSSA* 103(4), 2520–2533, doi:10.1785/0120120339. The
full texts could not be obtained for this review [not verified], so no
selection thresholds are attributed to them here. What is verifiable
from the GeoNet/data README [verified]: moment tensors have been
computed since August 2003 for M > ~4, currently by an analyst;
solutions from 2003-08-21 to 2020-06-18 used the Dreger/Berkeley code
(Method 1), and all solutions since 2020-06-18 use Herrmann's Computer
Programs in Seismology (Method 2); the published quality fields are
station count, %DC and VR, with no stated publication threshold; for
Dusky Sound 2009 and Kaikōura 2016 the catalogue carries USGS W-phase
solutions because a regional solution could not be obtained. Note that
the README's DOI for Ristau (2013) is incorrect (it gives Dreger &
Helmberger 1993); the correct DOI is 10.1785/0120120339.

**ADOPTED**: the velocity models (Ristau 2008, Table 1); "no solution
below the floor" as a legitimate outcome. Obtaining the Ristau papers'
own selection rules remains the most valuable open item.

## 3. Herrmann / CPS — the engine under both GeoNet and this pipeline

Herrmann (2013), *SRL* 84, 1081–1088, doi:10.1785/0220110096, with the
CPS moment-tensor course notes and the routine SLU solution pages
[verified]. Regional traces to 700 km with interactive per-trace
quality control (P polarity consistent on Z and R, no P on T, Rayleigh
particle motion); redundant stations at the same distance and azimuth
removed; depth grid-searched at 1 km with the maximum fit taken; per-
station time shifts and VR published for every solution, the shifts
fitted to A + B cos(az) + C sin(az) to diagnose a mislocation; a
factor-of-two amplitude misfit flagged for inspection.

Herrmann, Benz & Ammon (2011), *Monitoring the earthquake source
process in North America*, *BSSA* 101(6), 2609–2625,
doi:10.1785/0120110095 (abstract only [not verified]) reports routine
regional moment tensors in the **0.02–0.10 Hz** passband down to
**Mw 3.7**, limited by station density and signal-to-noise. This is the
closest published analogue to the present pipeline.

**ADOPTED**: the 0.02–0.10 Hz band for small events; the Mw 3.7
processing floor; the group-velocity record window (SLU cut
dist/3.3 km/s ± 40 s; here 2.5 km/s plus a magnitude-dependent tail
for slower NZ paths); per-station time shift and VR as published
diagnostics. **DEFERRED**: the azimuthal fit of the shifts as a
mislocation test.

## 4. USGS NEIC — regional practice, and two ideas from the teleseismic stack

NEIC's regional moment tensors (Mwr, 0–10°) follow Herrmann, Benz &
Ammon (2011) above. The W-phase (Kanamori & Rivera 2008; Hayes, Rivera
& Kanamori 2009, *SRL* 80(5), 817–822; **Duputel, Rivera, Kanamori &
Hayes 2012**, *GJI* 189(2), 1125–1147,
doi:10.1111/j.1365-246X.2012.05419.x [verified]) and the SynDepth
depth-modelling tool (Yeck et al. 2025, *SRL* 96(6),
doi:10.1785/0220240372; code doi:10.5066/P924LDLT [verified]) are
**teleseismic** and **NOT APPLICABLE** as methods here. Two design
elements transfer:

- A strictly **two-stage** screen: rejections that need no forward
  model first (instrument response fit, pre-event noise against a noise
  model, completeness, a median amplitude screen rejecting traces below
  0.1× or above 3× the event median), then iterative rejection by
  misfit after a first inversion. Duputel et al. note the amplitude
  screen's known failure mode: it can reject a nodal station.
- Magnitude-scaled **time-shift caps** in SynDepth: 6, 8, 10 s for
  Mw < 5.5, 5.5–6.0, ≥ 6.0, applied by clipping the cross-correlation
  search; and printed depth diagnostics including "on grid edge" flags.

NEIC publishes no numeric quality threshold for any moment tensor;
solutions are released automatically as preliminary and promoted after
analyst review. See `neic.md` for the details.

**ADOPTED**: the two-stage architecture; the 8 s time-shift cap.
**ADAPTED**: the amplitude screen is applied on the high side only (a
station far above the network median indicates a broken response; a
station far below may be nodal and is kept). **DEFERRED**: a
noise-model (PSD) screen.

## 5. INGV Italy

Scognamiglio, Tinti & Michelini (2009), *BSSA* 99(4), 2223–2242,
doi:10.1785/0120080104 (full text not obtained); Scognamiglio et al.
(2012), *Ann. Geophys.* 55(4), doi:10.4401/ag-6159 [verified];
Scognamiglio et al. (2016), *GJI* 206(2), 792–806,
doi:10.1093/gji/ggw173 [verified]; the INGV TDMT quality legend
[verified]. Automatic TDMT at ML ≥ 3.5 within 6–10 minutes, analyst-
reviewed for the catalogue. Stations chosen from eight 45° sectors with
a distance–magnitude weight (closer stations below M4.5, more distant
above, to avoid tilt and saturation), refined by a geometric-
uniformity optimisation; SNR > 5 on a 500 s window; a common time shift
across the three components of a station; no per-station VR cut-off
(a reviewed solution retained a station at VR −1.2%); small events
inverted at 0.02–0.10 Hz, larger at 0.02–0.05 Hz. The quality code's
first character depends jointly on station count and VR:

| N stations | D | C | B | A |
|---|---|---|---|---|
| 1 | VR < 50% | ≥ 50% | — | — |
| 2 | < 30% | 30–70% | ≥ 70% | — |
| 3 | < 20% | 20–70% | ≥ 70% | — |
| 4 | < 20% | 20–40% | 40–60% | ≥ 60% |
| 5–8 | < 15% | 15–40% | 40–60% | ≥ 60% |
| > 8 | < 15% | 15–30% | 30–50% | ≥ 50% |

with a second character `a`/`b` for %DC ≥ 60 / < 60. The VR bar falls
as the station count rises, so a solution cannot improve its grade by
shedding stations.

**ADOPTED**: the grade table (rows 3, 4, 5–8, > 8) with the `a` rule as
a requirement for A/B; eight azimuth sectors; higher frequencies for
small events. **ADAPTED**: sectors are used to order candidates
round-robin rather than to cap them (see §6).

## 6. The ISOLA family — ISOLA, Bayesian ISOLA, Gisola, scisola, BMKG

Zahradník & Sokos (2018), *ISOLA code for multiple-point source
modeling — review*, doi:10.1007/978-3-319-77359-9_1 [verified];
Vackář, Burjánek, Gallovič, Zahradník & Clinton (2017), *Bayesian
ISOLA*, *GJI* 210(2), 693–705, doi:10.1093/gji/ggx158 [verified];
Triantafyllis, Sokos, Ilias & Zahradník (2016), *scisola*, *SRL* 87(1),
157–163; Triantafyllis et al. (2022), *Gisola*, *SRL* 93(2A), 957–966
(configuration and source [verified]); Halauwet et al. (2024), *GJI*
239(2), 1000–1020, doi:10.1093/gji/ggae309 (BMKG) [verified].

Key practice: magnitude-scaled distance windows (Gisola: 10–250 km for
M4–5, 40–300 km for M5–5.5, up to 700 km for M6+); eight 45° sectors
with a minimum number populated and a maximum per sector (Gisola: at
least 3 sectors, at most 2 per sector, 24 stations overall); a depth
grid of ± 31 km about the hypocentral depth (Gisola) or a free
space–time grid (Bayesian ISOLA); no SNR threshold in Bayesian ISOLA,
which instead weights by a noise covariance matrix; a defined SNR in
BMKG (RMS of 200 s after P over 200 s before, per component, threshold
2, components below it ignored); the MouseTrap detector for long-period
instrumental disturbances (Vackář, Burjánek & Zahradník 2015, *SRL*
86(2A), 442–450); a trust criterion in Bayesian ISOLA of VR > 0.5,
condition number < 8, DC > 50% and bounded posterior spreads; letter
grades from VR × station count in Gisola and BMKG, with the condition
number as a second digit.

Two cautions from this literature shaped the pipeline. Triantafyllis
et al. (2016) demonstrate that repeatedly re-inverting while keeping
only high-VR stations degrades the solution: ten stations at VR ≈ 0.7
with condition number ≈ 5 become two stations at VR 0.8–0.9 with
condition number > 10. And the "ten minimum shear wavelengths" rule
(Zahradník & Sokos 2018) limits the usable high-frequency corner with
distance.

**ADOPTED**: the magnitude-scaled distance window; the 16-station cap;
the per-component RMS SNR with threshold 2; the ± 30 km depth window;
the principle that a VR gate must depend on station count (via the
INGV table). **ADAPTED**: the SNR window is the record window actually
inverted rather than a fixed 200 s (a fixed 200 s diluted an M3.7
event's 60 s train to noise), and a station is usable when its best
component clears the bar, so a nodal component does not veto the
station; sectors order candidates round-robin without a per-sector
cap, so a one-sided offshore geometry is not starved. **DEFERRED**:
MouseTrap; the condition number; distance-dependent bandwidth;
noise-covariance weighting in place of thresholds.

## 7. NIED F-net (Japan)

Fukuyama, Ishida, Dreger & Kawai (1998), *Zisin* 51(1), 149–156,
doi:10.4294/zisin1948.51.1_149 [verified, Japanese]; Fukuyama & Dreger
(2000), *Earth Planets Space* 52, 383–392 [verified]; Kubo, Fukuyama,
Kawai & Nonomura (2002), *Tectonophysics* 356, 23–48 (full text not
obtained); the F-net method page [verified]. The most minimal
long-running design: hypocentral distance 50–400 km; at most three
stations, taken in order of increasing distance among those with good
data (maximum amplitude and completeness checks; no SNR); deliberately
no azimuthal requirement; time shifts by cross-correlation against the
Green's functions with no published bound; depth searched every 3 km
within ± 30 km of the JMA hypocentre; VR weighted proportional to
epicentral distance; publication when M > 3.5 and VR > 50%, calibrated
in 1998 against first-motion mechanisms. Fukuyama et al. (1998) state
that automatically detecting a station gone bad is difficult because
the misfit does not necessarily localise on that station, and that the
jackknife is the reliable way to remove anomalous data.

**ADOPTED**: the ± 30 km depth window about the hypocentre; distance
weighting (`weight = distance` in mttime); the jackknife as the
stability evidence in the grade; withholding a solution below the
floor rather than publishing it.

## 8. GFZ / gempa AUTOMT (SeisComP)

Documented defaults [verified, docs.gempa.de]: distance ≤ 70°; minimum
six stations; minimum overall fit 0.3; a per-station fit floor of 0.3
below which stations are removed **unless the minimum station count
has been reached**; SNR minima of 2–3 by wave type (definition not
documented); time-shift caps of 10 s (body), 30 s (surface), 45 s
(mantle), 60 s (W-phase); depth searched coarse-to-fine, unconstrained
by default; no azimuthal-gap gate.

**ADOPTED**: a fixed per-station fit floor that pruning cannot breach
below the minimum station count.

## 9. SED (ETH), SCARDEC (IPGP), AutoBATS (Taiwan), Ibero-Maghreb

- Bernardi, Braunmiller, Kradolfer & Giardini (2004), *GJI* 157(2),
  703–716: a first inversion at a fixed depth removes traces with
  normalised variance ≥ 0.8 or large re-alignment; depth uncertainty as
  the range within a 10% variance increase.
- Vallée et al. (2011), *SCARDEC*, *GJI* 184(1), 338–358 (teleseismic):
  one station per 10° azimuth bin chosen on SNR; depth free within
  ± 50 km of the reference; uncertainty as misfit within 10% of the
  optimum.
- Jian, Tseng, Liang & Huang (2018), *AutoBATS*, *BSSA* 108,
  doi:10.1785/0120170231: a defined spectral SNR (three-component
  average, 150 s windows either side of P, threshold 2.0); stations
  beyond 30 km; time shifts ± 2 s; depth ± 12 km about the hypocentre;
  three parallel station selections (azimuth, SNR, distance) whose
  agreement is the stability test; publication gates on ISO ≤ 20%,
  CLVD ≤ 30%, non-DC ≤ 40%.
- Stich, Ammon & Morales (2003), *JGR* 108(B3): excluded stations'
  predictions are still computed as a compatibility check; station
  permutations tested to detect a controlling station.

**ADOPTED**: depth uncertainty as the range of depths within 10% of the
maximum VR (Bernardi et al.; Vallée et al.); best-SNR ordering within
azimuth bins. **DEFERRED**: predictions for excluded stations; the
three-selection agreement test; a stricter shift cap.

## 10. Network geometry and resolvability

Johnson et al. (2016), *GJI* 206(1), 525–556, doi:10.1093/gji/ggw141
(synthetic and real tests with CPS: focal-plane constraint depends on
geometry and faulting style more than on station count; stations
within ~30° of each other add little; Dufumier & Cara 1995 recommend at
least three stations 60° apart); Zahradník & Custódio (2012), *BSSA*
102(3), 1235–1254 (resolvability maps from geometry, model and band
alone); Ford, Dreger & Walter (2010), *BSSA* 100(5A), 1962–1970
(network sensitivity solutions); Sokos & Zahradník (2013), *SRL* 84(4),
656–665 (FMVAR/STVAR); Cesca & Heimann (2018),
doi:10.1007/978-3-319-77359-9_7 (ISO/CLVD trade-off at shallow depth).

## Comparison across the six questions

| | Pre-inversion vs fit-based rejection | Signal metric | Time shift | Distance and azimuth | Depth | Publication gate |
|---|---|---|---|---|---|---|
| SCSN | both; sensor corner, distance, clipping; then IVR loop with reselection | none (structural) | residual vs 0.13r + 1.42 s, drop if > 9 s | 45–700 km, six sectors | fixed set, max VR | A: 6 stations, OVR > 60 |
| GeoNet | [not verified] | [not verified] | [not verified] | [not verified] | [not verified] | NS/DC/VR reported; analyst |
| Herrmann/CPS | interactive per-trace QC | physical checks; amplitude factor 2 flag | fitted to A + B cos + C sin | < 700 km | 1–29 km at 1 km, max fit | analyst |
| NEIC W-phase (teleseismic) | two-stage: response, PSD, completeness, median 0.1–3×; then misfit ladder | PSD vs noise model | single centroid delay | 5–90°; N ≥ 30, gap ≤ 270° for the full search | inverted, ± 50 km, 12 km floor | analyst review |
| INGV | sectors + uniformity; no per-station cut | SNR > 5 | common per station, unbounded | 8 sectors | grid search, max VR | VR × N table |
| Bayesian ISOLA | pre only: MouseTrap, gaps, < 2 km | none; noise covariance | centroid time | r < 2^(2 ML) km | free grid | VR > 0.5, CN < 8, DC > 50 |
| Gisola | pre only: clip, PPSD, MouseTrap | optional STA/LTA + SNR 5 | centroid time grid | 10–700 km by M; ≥ 3 of 8 sectors, ≤ 2 each | ± 31 km about hypocentre | VR × N letter, CLVD digit |
| BMKG | SNR < 2 per component, MouseTrap | RMS 200 s after/before P | ± 3 s | 8 sectors, 3–8 stations | 5–40 km, 15 km radius | VR × N letter, CN digit |
| F-net | amplitude, completeness | none | cross-correlation, unbounded | 50–400 km, ≤ 3 stations, no azimuth rule | JMA ± 30 km at 3 km | M > 3.5, VR > 50 |
| gempa | post: station fit < 0.3, floored at 6 stations | SNR 2–3 | 10–60 s by wave type | ≤ 70°, no gap gate | free, coarse-to-fine | fit ≥ 0.3, ≥ 6 stations |
| AutoBATS | SNR ≥ 2 | spectral, 150 s windows | ± 2 s | > 30 km; 3–7 stations, three selections | ± 12 km at 1 km | ISO/CLVD/non-DC limits, misfit classes |
| **this pipeline (v5)** | unusable data, SNR, amplitude (high side); then the fit loop | RMS over the inverted window, best component ≥ 2 | ± 8 s cap | 120–300 km by M; 8 sectors round-robin to 16 | GeoNet ± 30 km, max VR, 10% band | INGV table + DC ≥ 60 + jackknife ≤ 25° |

## What was not verified

Ristau (2008, 2013) full texts; the "VR > 65%" filter used in some NZ
studies (not traceable to Ristau's text); INGV's depth grid; Kubo et
al. (2002) full text; Bernardi et al. (2004) Table 2; Pasyanos et al.
(1996); gempa's SNR definition; Herrmann, Benz & Ammon (2011) beyond
the abstract. Readers with access to these are invited to correct this
document.
