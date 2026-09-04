# Station selection and QC across operational MT systems

Reviewed 2026-09-05. Twelve systems: Berkeley/SCSN (the lineage this
pipeline descends from), GeoNet/Ristau, Herrmann CPS, USGS W-phase,
INGV, ISOLA / Bayesian ISOLA / Gisola / scisola / BMKG, NIED F-net,
GFZ-gempa, SED ZUR_RMT, SCARDEC, AutoBATS Taiwan.

[verified] = read in the primary source, code or config; [secondary] =
via a citing paper; [not verified] = could not confirm, flagged rather
than guessed.

---

## THE PAPER WE SHOULD HAVE READ FIRST

**Clinton, Hauksson & Solanki (2006), BSSA 96(5), 1689-1705,
doi:10.1785/0120050241** — the Dreger TDMT code, automated for southern
California, with every rule written down. Its architecture is almost
exactly our funnel, published twenty years ago. [verified]

Their loop, verbatim: *"(1) gather all available data and perform
initial quality control; (2) select stations for the initial inversion
by choosing waveforms from different azimuths at optimal distance for
best signal-to-noise ratio; (3) perform inversion; (4) if the inversion
results do not satisfy the desired quality, reject stations with the
poorest waveform fits, and select new stations subject to the
constraint of maximizing the available azimuthal distribution; repeat
(2) and (3); (5) if the station list is exhausted... the quality
required is relaxed."*

- Pre-inversion hard screens: sensor corner period >= 100 s; distance
  **45-700 km**; peak amplitude **< 80% of the clip level**.
  **No SNR threshold at all** - handled structurally.
- Geometry: **six azimuthal sectors**, subdivided until six are
  populated; first pick per sector is the station **nearest 60 km**.
- Bands by ML: <4.2 -> 10-50 s; 4.2-5.5 -> 20-50 s; >5.5 -> 20-100 s.
- Depth: fixed trial set 5, 8, 11, 15, 18, 21 km, max overall VR.
- Grades (OVR = overall VR, IVR = individual station VR):
  **A+** 6 stations, OVR > 85%, reject IVR < max(OVR-10, 75);
  **A** 6 stations, OVR > 60%, reject IVR < max(OVR-10, 50)
  -> PUBLISHED WITHOUT REVIEW;
  **B** 4 stations (max IVR per quadrant), OVR > 40% -> Mw only,
  mechanism *"not considered stable enough for distribution"*;
  **C** 4 highest-IVR stations, azimuth ignored -> not distributed.

**Their time-shift rule is better than ours.** They regress the
cross-correlation shift against distance over all quality-A solutions:
**ZCOR_expected = 0.13 r + 1.42 s** (r in km), and *"if the difference
between the observed and expected individual station ZCOR is greater
than 9 sec, that station is automatically removed"*. Motive: teleseismic
energy from a distant large event can otherwise be modelled as a local
source (their case: a local ML 3.0 given Mw 4.45 at 39% VR because of
an Mw 5.7 some 2000 km away).

---

## THE DEPTH QUESTION, ANSWERED

Most systems do NOT tie centroid depth to the hypocentre - but nearly
every one BOUNDS THE SEARCH around it, and one ties it outright.

| system | depth search |
|---|---|
| NIED F-net | **JMA hypocentral depth +/- 30 km**, 3 km steps (ties it) |
| Gisola | **+/-31 km about the hypocentral depth**, 2 km, floored 1 km |
| USGS W-phase | 3-D grid about the PDE, **+/-50 km**, min centroid **12 km** |
| SCARDEC | z_n **+/-50 km** about NEIC, floor 12 km (4-6 km after 2016) |
| AutoBATS | **+/-12 km about the CWB hypocentre**, 1 km |
| BMKG | absolute **5-40 km**, 2.5 km |
| SCSN | fixed set 5/8/11/15/18/21 km |
| SLU / Herrmann | 1-29 km at 1 km |
| INGV | grid search, range/step not published |
| **ours** | **full 1-58 km, unbounded** |

We are more permissive than every operational system reviewed. The
resolution of "independent vs forced to GeoNet" is therefore: **bound
the search, never force the answer**. A bounded window is the norm, not
a compromise of independence. Our free grid plus a >8 km downgrade flag
remains defensible, but it is a deliberate DEPARTURE from practice and
should be presented as one.

Caveats worth knowing: W-phase's 12 km floor is a known artefact source
(shallow events pile up on it; an Mww depth of 11.5 km means "at or
below the floor", not "resolved"). And at very shallow depths ISO and
CLVD become nearly indistinguishable (Cesca & Heimann 2018), so a %DC
gate does not mean the same thing at 2 km as at 20 km.

---

## THE SHARPEST CRITICISM OF OUR DESIGN

**Triantafyllis, Sokos, Ilias & Zahradnik (2016), SRL 87(1), 157-163**
(scisola), verbatim:

> *"Stations with large amplitudes may have an instrumental
> disturbance; these signals could be well fitted, but actually just
> these should be removed from the inversion (Zahradnik and Plesinger,
> 2005, 2010). This is why we do not support the idea of repeating the
> automatic inversion by simply keeping only stations with large VR."*

Their demonstration, a shallow Mw ~4 with 10 stations inside 100 km:

| stations | VR | condition number | FMVAR | verdict |
|---|---|---|---|---|
| 10 | ~0.7 | ~5 | ~10 deg | good |
| 2 | **0.8-0.9** | **>10** | 20-40 deg | bad |
| 1 | ~0.98 | - | - | *"always dangerous"* |

**VR is anti-correlated with station count, so any gate phrased purely
as "VR > X" is gameable by our own funnel.**

**The cleanest published defence is INGV\'s quality table**, read from
their legend (terremoti.ingv.it/en/help#TDMT) - the VR bar FALLS as
station count rises:

| N stations | D | C | B | A |
|---|---|---|---|---|
| 1 | VR < 50% | >= 50% | - | - |
| 2 | < 30% | 30-70% | >= 70% | - |
| 3 | < 20% | 20-70% | >= 70% | - |
| 4 | < 20% | 20-40% | 40-60% | >= 60% |
| 5-8 | < 15% | 15-40% | 40-60% | >= 60% |
| > 8 | < 15% | 15-30% | 30-50% | >= 50% |

Second character: `a` = %DC >= 60, `b` = %DC < 60. Structurally,
**grade A is unreachable with <= 3 stations and B with 1 station**.
Gisola encodes the same idea in code (a 2-station solution needs
VR >= 0.7 for B; a 5-station solution needs 0.6 for A), and gempa puts
a floor under pruning: `minimumFinalStationFit = 0.3`, *"Stations are
removed below this threshold unless the minimum station count has been
reached"*.

Since our funnel deliberately shrinks the station set, our flat VR gate
is the wrong shape.

---

## WHERE WE ALREADY MATCH PUBLISHED PRACTICE

1. **Two-stage architecture** (hard rejects needing no forward model,
   then fit-based pruning) = Duputel et al. 2012 sec 3.2, Clinton 2006.
2. **The funnel** = Clinton\'s five-step loop; described almost word for
   word for SED\'s Dreger workflow in Vackar et al. 2017: *"select the
   optimal set of stations and event depth that produces the best MT,
   taking into account the Variance Reduction and the % Double Couple,
   whilst retaining as many stations as possible."*
3. **Worst-station VR as a gate** = Clinton\'s IVR floor and gempa\'s
   `minimumFinalStationFit = 0.3`.
4. **Depth searched, not fixed** - universal.
5. **Core forced to span >= 90 deg** = a weaker form of SCSN\'s six
   sectors, Gisola\'s 3-of-8, INGV/scisola/BMKG\'s 8 sectors.
6. **Grading on VR + %DC + jackknife** = Gisola\'s letter+digit and
   ISOLA\'s VR/CN/FMVAR/STVAR.
7. **10-300 km magnitude-scaled** brackets Gisola\'s 10-250 km (M4-5)
   and 40-300 km (M5.1-5.5).
8. **Cross-correlation off in the survey pass** is mttime\'s own default
   (`correlate 0` in the shipped BSL example), and necessary because
   `Inversion._correlate()` searches the entire lag range unbounded.
9. **Distance weighting** (`weight = distance`): F-net\'s published VR
   definition weights stations *"proportional to the epicentral
   distance"* - our choice has an operational precedent.
10. **Small events pushed to HIGHER frequencies**: INGV inverts
    ML >= 3.8 at 0.02-0.05 Hz but *"Lower magnitude earthquakes were
    inverted in the frequency band of 0.02-0.1 Hz"* - exactly our
    10-50 s for the smallest events.
11. **The jackknife as the instrument for detecting a bad station**:
    NIED state in print that automatic detection is unsolved *because
    the misfit does not necessarily localise on the offending station*,
    and conclude the jackknife is the only way (Fukuyama et al. 1998
    sec 5). Strong support for grading on jackknife stability rather
    than per-station VR alone.
12. **"No coherent solution" as an outcome**: GeoNet themselves fall
    back to USGS W-phase for Dusky Sound 2009 and Kaikoura 2016
    *"as a reliable regional moment tensor solution could not be
    calculated"*.

---

## RANKED CHANGES WORTH MAKING

1. **Make the VR threshold a function of station count** (INGV table
   above; Gisola in code). One-line change, citable precedent, and the
   direct antidote to our funnel inflating VR by shrinking the set.
2. **Replace the absolute 8 s shift bound with a residual-against-
   expectation bound** (SCSN): fit our own ZCOR-vs-distance regression
   on accepted solutions, drop stations whose residual exceeds ~9 s.
   Then fit survivors to **A + B cos(az) + C sin(az)** (Herrmann\'s
   routine practice = GeoNet\'s current code family) and report the
   implied origin-time and epicentre correction - turning time shifts
   from a nuisance into a mislocation diagnostic. Reference points:
   gempa 10 s body / 30 s surface; AutoBATS **+/-2 s "following Dreger
   (2003)"**; ISOLA treats a shift *"> 10 seconds"* as an indication the
   velocity model is inapplicable.
3. **Adopt a defined per-COMPONENT SNR and retire peak/noise 1.2.**
   Every system that specifies a unit screens per component, not per
   station - NEIC routinely keeps a vertical and drops both horizontals
   *"the horizontal components usually being noisier"*. Best specified:
   AutoBATS (3-component average point-by-point spectral SNR, 150 s
   windows either side of P, 5-point smoothing, **threshold 2.0**);
   simplest: BMKG (**RMS 200 s after P / RMS 200 s before P, per
   component, threshold 2**); INGV uses **SNR > 5** on a 500 s window.
   Ours at 1.2 is "signal 20% above noise", which no published system
   would accept.
4. **Make the usable bandwidth DISTANCE-DEPENDENT** - the most repeated
   advice in the ISOLA literature and absent from our pipeline.
   Zahradnik & Sokos: waveforms cannot be modelled beyond ~**10 minimum
   shear wavelengths**, so at beta = 3 km/s, 0.1 Hz is usable only to
   ~300 km. Bayesian ISOLA: for distance > 100 km, cap the high corner
   so the minimum wavelength is no shorter than distance/5. Pooling
   10-300 km in one fixed band under-uses our near stations and
   over-trusts our far ones.
5. **Put a floor under the funnel that pruning cannot breach**
   (gempa\'s rule), because of the scisola criticism.
6. **Add the CONDITION NUMBER** - the one diagnostic our grading lacks,
   and precisely the one that detects "high VR from too few or too
   clustered stations". CN > 5-10 => ill-posed (Zahradnik & Sokos);
   Bayesian ISOLA requires **CN < 8**, trusted solutions cluster at
   2-4; BMKG\'s top grade is CN < 5. One SVD of G.
7. **Long-period disturbance detection (MouseTrap).** Step-like
   instrumental disturbances survive bandpass filtering, carry huge
   false amplitudes and FIT WELL - so neither our amplitude screen nor
   our VR-based funnel catches them. Every ISOLA-family system runs it
   by default. Open source Python. Vackar, Burjanek & Zahradnik 2015
   SRL 86(2A) 442-450, doi:10.1785/0220140168.
8. **Clipping screen** - absent from our hard rejects. SCSN: reject if
   peak amplitude exceeds **80% of the clip level**, because
   *"broadband sensors can have a nonlinear response significantly
   short of their expected clip level"*. F-net swaps to a strong-motion
   sensor instead. Matters for our 10 km stations.
9. **Two-side the amplitude screen?** Duputel et al. 2012 reject if
   peak-to-peak **< 0.1x or > 3x the event median**. NOTE THE TENSION:
   the authors state its known failure mode is *"can accidentally
   reject some good data (e.g. a nodal station)"* - which is exactly
   the objection that made us go one-sided on 2026-09-04. Both
   positions are defensible and the nodal case is acknowledged in the
   literature, not denied. Needs a decision, not a default.
10. **Near-field stations: invert both WITH and WITHOUT** (ISOLA\'s
    prescription) rather than only tagging. Our 10 km inner limit is
    far closer than SCSN\'s 45 km, AutoBATS\'s 30 km *"to reduce the
    mislocation effect"*, F-net\'s 50 km or Bayesian ISOLA\'s 2 km.
11. **Decluster by azimuth, not distance.** SCARDEC: *"When several
    stations are present in a 10 deg azimuthal range, we only select
    the one with the best signal-to-noise ratio."* Gisola: max 2
    stations per 45 deg sector. INGV goes further with an explicit
    geometric-uniformity optimisation (graph theory, minimising the
    spread of distances and adjacent azimuth gaps; Scognamiglio et al.
    2012, Ann. Geophys. 55(4), doi:10.4401/ag-6159). This would also
    strengthen our 90 deg span requirement - 90 deg spanned by two
    clusters is not 90 deg spanned by two stations.
12. **Keep excluded stations as PREDICTIONS.** Stich et al. 2003:
    excluded stations *"predictions are still calculated to confirm a
    basic compatibility with the obtained moment tensor solution"*.
    Nearly free (we have the Green\'s functions) and turns every
    exclusion into an independent check - ideal for a system that has
    to justify itself in an email.
13. **Bound the depth grid** (table above).
14. **Two cheap validation ideas**: AutoBATS runs **three station
    selections in parallel - best azimuth, best SNR, shortest
    distance - and treats their agreement as the stability test**;
    SCARDEC and Bernardi et al. both define parameter uncertainty as
    *"misfit not exceeding the optimum by more than 10%"*.

---

## SYSTEM NOTES

**NIED F-net** - the most minimal design, and instructive for that.
Hypocentral distance **50-400 km hard**, then *"we used at most three
station whose hypocentral distance is between 50 km and 400 km and
whose data quality is good"*, chosen by increasing distance. **No
azimuthal requirement, deliberately** - the 1998 conclusion is that the
solution is recovered accurately even when stations do not surround the
source, which is why it suits offshore trench events. Screening is
maximum amplitude plus completeness; no SNR. Time shift by 8
cross-correlations against the Green\'s function components, taking the
largest absolute lag; no published maximum. Publication gate, explicit:
*"We only show reliable solutions which satisfy the criteria that
magnitude is greater than 3.5 and quality (variance reduction) is
greater than 50 %"* - calibrated in 1998 against first-motion
mechanisms via P-wave radiation-pattern cross-correlation (for M >= 4
and VR >= 50%, correlation >= 0.7). NIED\'s own caveat: *"When you use
automatically determined moment tensor solutions, due to unexpected
noise the result might happen to be completely wrong."*

**INGV** - AUTO-TDMT at ML >= 3.5, solution in 6-10 min, *"automatically
published on the World Wide Web for solution qualities exceeding a
predefined threshold"*, plus an analyst-reviewed REV-TDMT catalogue.
8 x 45 deg sectors with a distance-magnitude weight that *"avoids very
close stations that may be affected by equipment tilting and
significant centroid location errors, and privileges stations at
greater distances for larger magnitudes to avoid using saturated
waveforms"*. **No per-station VR cut-off**: a reviewed 7-station
solution retained a station at **VR -1.2%** - our unconditional
anti-fitting cull is stricter than INGV. SNR > 5 on a 500 s window.
Time shift *"forced to be the same for the three station components"*,
no published maximum.

**GeoNet** [verified from the GeoNet/data README] - MT computed since
2003-08 for M > ~4, and *"At the moment these are computed manually"*:
we are automating something GeoNet does by hand. **Method 2, all
solutions since 2020-06-18, is Herrmann\'s CPS** - the same family as
our Green\'s functions, so Herrmann\'s documented practice is the
closest thing to a house style we have. Published quality fields are
only NS, DC and VR, with no stated publication thresholds.
WARNING: the README prints the WRONG DOI for Ristau (2013) - it gives
10.1029/93JB00023, which is Dreger & Helmberger 1993. The correct DOI
is **10.1785/0120120339**.

**USGS W-phase** - the most rigorously specified screening pipeline,
strictly two-stage. Pre-inversion: instrument response fit within 3%;
pre-event PSD rejected against the New High Noise Model over 1-10 mHz;
completeness; **median screening, reject if peak-to-peak < 0.1x or >
3x the event median**. Post-inversion: iterative misfit rejection at
**three thresholds 3.0 -> 2.0 -> 1.0**, discarding on average **50% of
channels**. The rejection unit is the CHANNEL, never the station.
Single global centroid delay; no per-station shifts.

---

## BIGGEST GAP - ACTION FOR DANI

**Ristau (2008) SRL 79(3) 400-415 and Ristau (2013) BSSA 103(4)
2520-2533 could not be obtained** (SRL/BSSA 403; the VUW open-access
thesis and its DSpace mirror behind Cloudflare; the Springer 2018
chapter behind auth). We therefore cannot say what our own reference
catalogue\'s published station-selection and quality rules actually
are. With institutional access this is the single highest-value half
hour available. Also unverified: the "VR > 65%" filter used by
downstream NZ studies could not be traced to Ristau\'s own text.

Other gaps: Kubo et al. 2002 full text; Bernardi et al. 2004 Table 2;
Pasyanos et al. 1996 (the Berkeley PDF is 404); gempa\'s SNR definition
(thresholds documented, metric not); INGV\'s depth range and step.

---

## KEY REFERENCES

- Clinton, Hauksson & Solanki (2006). BSSA 96(5), 1689-1705.
  doi:10.1785/0120050241
- Dreger & Helmberger (1993). JGR 98, 8107-8125. doi:10.1029/93JB00023
- Pasyanos, Dreger & Romanowicz (1996). BSSA 86(5), 1255-1269.
- Dreger (2003). TDMT_INV. IASPEI Handbook 81B, 1627.
- Ristau (2008). SRL 79(3), 400-415. doi:10.1785/gssrl.79.3.400
- Ristau (2013). BSSA 103(4), 2520-2533. doi:10.1785/0120120339
- Herrmann (2013). CPS. SRL 84, 1081-1088. doi:10.1785/0220110096
- Herrmann, Malagnini & Munafo (2011). BSSA 101(3), 975-993.
  doi:10.1785/0120100184
- Duputel, Rivera, Kanamori & Hayes (2012). GJI 189(2), 1125-1147.
  doi:10.1111/j.1365-246X.2012.05419.x
- Scognamiglio, Tinti & Michelini (2009). BSSA 99(4), 2223-2242.
  doi:10.1785/0120080104
- Scognamiglio et al. (2012). Ann. Geophys. 55(4). doi:10.4401/ag-6159
- Scognamiglio, Magnoni, Tinti & Casarotti (2016). GJI 206(2), 792-806.
  doi:10.1093/gji/ggw173
- Fukuyama, Ishida, Dreger & Kawai (1998). Zisin 51(1), 149-156.
  doi:10.4294/zisin1948.51.1_149
- Fukuyama & Dreger (2000). Earth Planets Space 52, 383-392.
- Kubo, Fukuyama, Kawai & Nonomura (2002). Tectonophysics 356(1-3),
  23-48. doi:10.1016/S0040-1951(02)00375-X
- Zahradnik & Sokos (2018). Springer. doi:10.1007/978-3-319-77359-9_1
- Vackar, Burjanek, Gallovic, Zahradnik & Clinton (2017). GJI 210(2),
  693-705. doi:10.1093/gji/ggx158
- Triantafyllis et al. (2016). scisola. SRL 87(1), 157-163.
- Triantafyllis et al. (2022). Gisola. SRL 93(2A), 957-966.
- Halauwet et al. (2024). GJI 239(2), 1000-1020. doi:10.1093/gji/ggae309
- Bernardi, Braunmiller, Kradolfer & Giardini (2004). GJI 157(2),
  703-716. doi:10.1111/j.1365-246X.2004.02215.x
- Vallee, Charlety, Ferreira, Delouis & Vergoz (2011). GJI 184(1),
  338-358. doi:10.1111/j.1365-246X.2010.04836.x
- Jian, Tseng, Liang & Huang (2018). AutoBATS. BSSA 108.
  doi:10.1785/0120170231
- Stich, Ammon & Morales (2003). JGR 108(B3).
- Johnson, Hayes, Herrmann, Benz, McNamara & Bergman (2016). GJI
  206(1), 525-556. doi:10.1093/gji/ggw141
- Zahradnik & Custodio (2012). BSSA 102(3), 1235-1254.
  doi:10.1785/0120110216
- Ford, Dreger & Walter (2010). BSSA 100(5A), 1962-1970.
  doi:10.1785/0120090140
- Sokos & Zahradnik (2013). SRL 84(4), 656-665. doi:10.1785/0220130002
- Vackar, Burjanek & Zahradnik (2015). MouseTrap. SRL 86(2A), 442-450.
  doi:10.1785/0220140168
- Zahradnik & Plesinger (2005). BSSA 95(5). doi:10.1785/0120040210
- Cesca & Heimann (2018). Springer. doi:10.1007/978-3-319-77359-9_7
