# Worked examples — Fiordland sequence, 2026-09-02

Four events near Milford Sound at the southern termination of the
Alpine Fault, processed automatically with the Ristau (2008) South
Island model. GeoNet initial depths were 5 km placeholders for all
four. Full archives (solution.json with complete provenance, all
figures, per-band results) are under `events/<publicID>/`.

| publicID | GeoNet prelim M / depth | our Mw / depth | VR | %DC/%CLVD | plane1 s/d/r | stations | band | pred. peak disp | email gate |
|---|---|---|---|---|---|---|---|---|---|
| 2026p660242 | M5.6 / 5 km | **5.03** / 9 km | 78% | 92/8 | 3/53/38 | 11 | 30-100s | 0.32 cm | no — daily cap 3 reached |
| 2026p660160 | M5.5 / 5 km | **5.05** / 5 km | 72% | 53/47 | 10/36/57 | 16 | 20-50s | 1.23 cm | no — daily cap 3 reached |
| 2026p660272 | M4.8 / 5 km | **4.22** / 5 km | 47% | 73/27 | 178/63/-24 | 7 | 10-50s | 0.05 cm | no — quality grade D (email requires A or B): {'min_stations': Tr |
| 2026p660321 | M4.6 / 5 km | **4.21** / 1 km | 65% | 88/12 | 181/86/-61 | 11 | 20-50s | 1.73 cm | no — daily cap 3 reached |

Per event, the email carries three figures (copied here):
`*_waveform_fits.jpg` (mttime fits with the Deviatoric = DC + CLVD
decomposition), `*_stations_displacement_field.jpg` (station map + Okada E/N/U
predicted displacement + NISAR passes), and
`*_depth_sensitivity.jpg` (VR/%DC/Mw vs depth).
