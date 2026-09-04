# Review learnings — 2026-09-03 figure audit

Synthesis of the manual review of the 2026 catalogue (Danielle Lindsay,
with Claude assistance). Three parts: what we now know is true, what to
watch for when reviewing any solution, and the fix backlog with the events
that motivate each item.

## 1. Learnings (established, with evidence)

1. **Metrics must measure what the eye measures.** RMS SNR could not
   reproduce expert keep/toss calls (2026p283255: keeps and drops
   interleaved 0.7–1.2); peak-to-noise inside the inverted window could
   (16/18 agreement). Same lesson at band level: VR is not comparable
   across bands — smoother long-period data posts higher VR while fitting
   noise (2026p033598: "terrible" 20–50 s at VR 61 beat signal-fitting
   10–50 s at VR 40, with Mw +0.22 inflation).
2. **Every automated judgment is only as good as its reference
   solution.** Admission, culling and per-station VR are computed against
   the current model; when the core solution is junk (depth at grid edge,
   VR < ~20), verdicts invert — the best station gets culled and a
   dead channel kept (2026p300334: PYZ culled at VR −149, flat-Z WHZ
   used; 2026p189537: four stations with visibly coherent pulses all
   rejected with negative VR against a 46 km core).
3. **Robustness hierarchy of the outputs**: strike/dip/rake is the robust
   deliverable; depth is soft; %DC is the softest — it is a derivative of
   depth (the DC-vs-depth see-saw converts ±3 km of depth scatter into
   ±40 DC points; Te Kaha twins: DC 52 at 3 km vs DC 98 at 7 km, 12°
   apart in mechanism) and inflatable by noise. Mw inflates when the band
   holds no signal (noise fitting), so Mw ≫ prelim is a red flag, not a
   discovery.
4. **Isolated depth spikes at velocity-model layer boundaries are not
   solutions.** The Green's functions come from a 1-D layered model, and
   a source moved across a layer interface changes its excitation
   abruptly — mode shapes, amplitude ratios and depth phases all jump.
   The depth search samples straight across these jumps, and the
   per-station time-shift search (zcor) can slide traces into chance
   alignment on one side of an interface, producing a single depth with
   anomalously good VR (or a DC that snaps to ~100) surrounded by much
   worse neighbours (2026p101368 and others). These spikes are artifacts
   of the model's roughness, not information about the source. The
   interface depths where they congregate: South Island model 3, 12, 27,
   39 km; North Island model 1, 3, 5, 15, 33, 64 km (from
   `models/*.d`). Rule of thumb: a depth is believable when it sits on a
   CONTIGUOUS plateau whose mechanisms morph smoothly through it (the
   depth-sensitivity strip shows this at a glance); a one-point spike —
   especially adjacent to an interface depth or at the grid edge — is to
   be distrusted no matter how good its VR or DC looks. The
   plateau-contiguity pick rule automates part of this, but reviewers
   should still check spikes against the interface list when a depth
   looks surprising. This is also why %DC see-saws ±40 points over ±3 km
   near interfaces (the depth-quantization CLVD effect).
5. **The letter grades quarantine correctly.** Every pathology found in
   review was already C/D; nothing broken could have reached the email
   gate. Archive honesty, not publication safety, is what review
   improves.
6. **Aftershock sequences are consistent when the station set is.** The
   five good-grade Te Kaha solutions (11–19° from the mainshock) share a
   spine (KUZ, HAZ, KRBZ, GRZ, OPRZ); the deviants (75–100°) used almost
   disjoint sets, selected during the mainshock's coda. Minutes-to-an-hour
   after a big event, long-period coda poisons both selection and moment
   (2026p336046 at +5 min: prelim 4.2 → Mw 5.01, VR 24).
7. **Sparse events currently disable their own safety rails.** The
   ≥3-station floor blocks the greedy loop (and with it the anti-fitting
   cull), and the amplitude screen needs ≥4 stations — so the weakest
   events get the least protection (2026p101368: RDHZ used at VR −40 with
   zcor at the +100 s limit).
8. **Physics limits are not bugs.** The closest stations can carry the
   sharpest-looking records and still be unfittable: at 10–50 s a 36 km
   station is ~one wavelength out, where waveforms are most sensitive to
   model error (2026p238013 BSWZ/WMVZ at VR −73/−20 while 55 km stations
   fit). And an offshore M4 with no coherent energy at any station cannot
   be inverted at all (2026p111636) — the honest product is "no
   solution", not a grade-D artefact.
9. **Selection rules compound.** 3×depth exclusion + magnitude-scaled
   radius + dead-channel floor + cluster thinning can jointly starve an
   easy event (2026p355073: WEL/BSWZ/WLCZ — the best stations — all
   removed by 3×depth for an 18 km-deep event; 2026p348732 similar).
   The 3×depth rule is a body-wave far-field heuristic; the CPS Green's
   functions are complete (near-field terms included), so 3× is likely
   far too conservative.

## 2. Reviewer watch-list (per-solution red flags)

- Depth at the grid edge, especially with Mw ≫ prelim → band/noise
  mismatch; distrust everything downstream (DC, depth, Mw).
- All-station figure: candidates rejected *en masse* with negative VRs
  while visibly sharing coherent arrivals → the core solution is wrong,
  not the stations.
- Any used station with own VR ≤ 0 (a "passenger") — check how it
  survived (sparse floor, sector-sole protection).
- zcor at the search limit (±100 s) → chance alignment, not fit.
- DC differing wildly between neighbouring events or depths → the
  depth/CLVD see-saw, not source physics.
- A single-depth VR or DC spike adjacent to a model interface (SI: 3,
  12, 27, 39 km; NI: 1, 3, 5, 15, 33, 64 km) or at the grid edge →
  layer-boundary artifact; trust plateaus with smoothly-morphing
  mechanisms, never one-point winners.
- Event within ~1 h and ~100 km of an M5.5+ → coda contamination of both
  selection and moment.
- Beautiful close-station records rejected at long period → near-field
  physics, not selection error.
- Sequence beachballs "looking different" → check tensor angle first;
  plane-1/plane-2 labelling flips make conjugates look unlike.
- One-component-only signal (per-component peaks in the corner of the
  all-station figure) — the shared per-station scale makes the other
  components look flat; that is real relative amplitude, not a bug.

## 3. Fix backlog (motivating events in brackets)

**Status 2026-09-05: items 1-6, 8, 9, 13 and 14 are IMPLEMENTED in
station selection v4** (see METHOD §3.1). The funnel replaced the
pre-filter architecture entirely: nothing usable is deleted, the
3x-depth rule and cluster thinning became demotion tags, the
dead-channel floor dropped to 1.2, the amplitude screen became
one-sided (nodal-plane insight), the greedy VR-polishing pass was
removed, bad cores can no longer veto data (pass 1 has no time shifts
and the core must span 90 deg), the grid-edge depth guard is in,
signal-aware windows are in, band menus escalate, and events that
cannot be constrained now abort with "no coherent solution" instead of
producing a junk mechanism. Remaining: 7 (sequence mode), 10 (degree-6
/ source-type screen), 11 (coda guard), 12 (component marking on the
figures).

**New learning (2026-09-05): a station that sees no signal is not
necessarily a bad station.** If it sits near a nodal plane its
amplitude is genuinely small, and that small amplitude is information
about the mechanism. Two consequences, both now implemented: the
amplitude screen only rejects stations far ABOVE the network median
(broken responses), never below; and low peak/noise only demotes,
never excludes.

**New learning (2026-09-05): mttime's time-shift search is unbounded**
(a 60 s window can slide past 100 s), which is the mechanism by which a
noise trace earns undeserved VR. Two consequences: the survey pass runs
with shifts disabled, so its ranking cannot be gamed; and from then on
a station whose solved shift exceeds 8 s is rejected. Archive
calibration: stations in grade A/B solutions sit at |zcor| <= 9 s
(p95), grade-D stations reach 0.88 of their whole travel time.

**New learning (2026-09-05): what a careful analyst catalogue looks
like.** Ristau's NZ CMT catalogue uses a median of 7 stations
(quartiles 4-11, sometimes 1-3), runs at a median VR of 72 with a 10th
percentile of 57, and achieves a median VR of 71 even on Mw<4.5 events
with 6 stations. His centroid depths sit a median 4 km from the GeoNet
hypocentre (57% within 5 km, 82% within 10 km). So: few stations are
normal, small events CAN be fit well (ours at VR 20-30 are
underperforming rather than at a physical limit), and a centroid depth
far from the hypocentre is claiming something a careful catalogue
rarely does.


1. **Mw-triggered band escalation**: if the inverted Mw exceeds the
   menu's bracket (e.g. ≥4.7 from a 10–50 s-only run), run the next
   band up and prefer it [336046, 348732].
2. **Relax the 3×depth near-field rule** (complete GFs justify ~1.5–2×,
   or cap the exclusion radius) [355073, 348732, 238013 context].
3. **Safety rails for sparse events**: let the anti-fitting cull override
   the 3-station floor (better 2 honest stations, or an abort, than 3
   with poison); run the amplitude screen at ≥3 [101368].
4. **Bad-core recovery**: when most candidates return negative VR,
   re-form the core (or run one joint search with everything) instead of
   trusting the initial core; declare failure if nothing coheres
   [189537, 300334].
5. **No-signal abort**: when no station clears a modest peak/noise bar,
   archive "no coherent long-period signal" instead of inverting noise
   [111636].
6. **Greedy own-VR shield**: never evict a station fitting above ~40–50
   just to polish joint VR; revisit sector-sole protection (HAZ cost
   2026p338701 DC 56→72, VR 60→67) [091845, 338701].
7. **Sequence mode**: within the aftershock window of a well-graded
   solution, inherit its station spine and use its depth as a prior
   [Te Kaha family].
8. **Grade rubric v2 — "no passengers" replaces station count**: A needs
   VR ≥ 70, gap ≤ 180, every used station fitting (own VR ≥ ~30), and
   jackknife rotation ≤ 15° when available; 3-station solutions can be
   A (BSL practice). Retroactive regrade is cheap (all inputs stored).
9. **Soften the dead-channel floor** toward ~1.5 now that admission
   exists to catch what slips through [348732: MRZ at 1.9].
10. **Degree-6 + source-type screen** for low-DC / TVZ-path events
    [336118].
11. **Coda guard**: within X h / Y km of M≥5.5, require stronger
    peak/noise or defer processing [336046, 336118].
12. **Figure**: mark the signal-carrying component per station (per-comp
    peak/noise annotation or bolding); keep the shared per-station scale
    but say so on the figure.

13. **Signal-aware window end**: East Coast offshore paths (Hikurangi
    accretionary prism) deliver 10–50 s trains at effective group
    velocities of ~1.2–1.7 km/s — the kinematic cut (dist/2.5 + tail)
    bisects the train at every station [385521: PUZ/TKGZ strongest
    arrivals inside the grey]. Fix: keep the kinematic value as the
    window-end MINIMUM, extend to where the smoothed envelope decays
    back toward pre-event level, cap at 150 s — deterministic,
    data-derived, path-adaptive.

14. **Grid-edge depth guard**: a VR maximum on the first/last grid depth
    with a zero-width plateau is an artifact (the smoothest GFs absorb
    noise); prefer the best INTERIOR local maximum within the VR
    tolerance, applying the normal plateau+DC rule there; keep the edge
    only when no interior candidate comes close [508890: edge 58 km at
    VR 22.7 beat the physical 8 km peak (VR 18.3, DC max 88–96 at
    9–10 km, GeoNet ~8 km); the contiguity guard, built against
    disconnected lobes, protected the artifact instead].

Items 1–3 are the highest-leverage: they address every "should be easy
but looks bad" case found in review; item 13 joins them for East Coast
events, and item 14 catches the grid-edge runaways (508890, 348732,
300334 all rode a 58 km edge max).
