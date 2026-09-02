# Worked examples — Fiordland sequence, 2026-09-02

Four events near Milford Sound at the southern termination of the
Alpine Fault, processed automatically with the Ristau (2008) South
Island model. GeoNet initial depths were 5 km placeholders for all
four. Full archives (solution.json with complete provenance, all
figures, per-band results) are under `events/<publicID>/`.

| publicID | GeoNet prelim M / depth | our Mw / depth | VR | %DC/%CLVD | plane1 s/d/r | stations | band | pred. peak disp | email gate |
|---|---|---|---|---|---|---|---|---|---|
| 2026p660242 | M5.6 / 5 km | **5.03** / 7 km | 85% | 97/3 | 359/55/31 | 6 | 30-100s | 0.47 cm | YES |
| 2026p660160 | M5.5 / 5 km | **5.07** / 8 km | 72% | 66/34 | 7/45/46 | 10 | 10-50s | 0.50 cm | YES |
| 2026p660272 | M4.8 / 5 km | **5.03** / 58 km | 87% | 97/3 | 272/72/156 | 1 | 20-100s | 0.01 cm | no — quality gates failed: {'min_stations': False, 'vr_floor': Tr |
| 2026p660321 | M4.6 / 5 km | **4.09** / 4 km | 80% | 94/6 | 279/68/-174 | 8 | 20-50s | 0.03 cm | no — Mw 4.09 < 5.0 and predicted displacement 0.03 cm < 1 cm |

Per event, the email carries three figures (copied here):
`*_waveform_fits.jpg` (mttime fits with the Deviatoric = DC + CLVD
decomposition), `*_stations_displacement_field.jpg` (station map + Okada E/N/U
predicted displacement + NISAR passes), and
`*_depth_sensitivity.jpg` (VR/%DC/Mw vs depth).
