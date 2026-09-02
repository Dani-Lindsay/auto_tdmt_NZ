"""Publishing: email the solution (figure attached) to the distribution list.

SMTP credentials come from the environment (GitHub Actions secrets in CI):
    SMTP_HOST, SMTP_PORT (587), SMTP_USER, SMTP_PASS, MAIL_FROM, MAIL_TO
MAIL_TO is ideally ONE address — a Google Group / list address — so no
subscriber addresses live in this public repo.

Bluesky/LinkedIn can be added later behind the same publish() interface.
"""

from __future__ import annotations

import os
import smtplib
from email.message import EmailMessage
from pathlib import Path

import config


def draft_text(solution: dict, forward: dict, passes: list[dict]) -> tuple[str, str]:
    """(subject, body) for one solution."""
    ev = solution["event"]
    pref = solution["preferred"]
    p1, p2 = pref["plane1"], pref["plane2"]

    subject = (
        f"[auto-MT NZ] Mw {pref['mw']:.1f} {ev['locality']} "
        f"({ev['public_id']}), depth {pref['depth_km']:g} km"
    )

    lines = [
        f"Automated regional moment tensor solution — PRELIMINARY",
        "",
        f"Event: {ev['public_id']}  {ev['origin_time']}",
        f"Location: {ev['locality']} ({ev['latitude']:.3f}, {ev['longitude']:.3f})",
        f"GeoNet preliminary magnitude: M {ev['prelim_mag']:.1f}",
        "",
        f"Mw: {pref['mw']:.2f}",
        f"Centroid depth: {pref['depth_km']:g} km (GeoNet initial: {ev['depth_km']:g} km)",
        f"Variance reduction: {pref['vr']:.0f}%  |  DC {pref['pdc']:.0f}% / CLVD {pref['pclvd']:.0f}%",
        f"Filter band: {1/solution['filter_band_hz'][1]:.0f}-"
        f"{1/solution['filter_band_hz'][0]:.0f} s  |  "
        f"Velocity model: {solution['provenance']['velocity_model']}",
        f"Nodal plane 1 (strike/dip/rake): {p1['strike']:.0f}/{p1['dip']:.0f}/{p1['rake']:.0f}",
        f"Nodal plane 2 (strike/dip/rake): {p2['strike']:.0f}/{p2['dip']:.0f}/{p2['rake']:.0f}",
        f"Stations used: {solution['quality']['n_stations_used']}",
        "",
    ]
    if forward["detectable"]:
        lines.append(
            f"Predicted peak surface displacement (Okada, both nodal planes): "
            f"{forward['peak_abs_m'] * 100:.1f} cm — potentially detectable "
            f"with InSAR."
        )
    else:
        lines.append(
            f"Predicted peak surface displacement {forward['peak_abs_m'] * 100:.2f} cm "
            f"— no detectable surface deformation expected."
        )
    if passes:
        lines.append("")
        lines.append("NISAR passes over the epicentre:")
        for p in passes:
            lines.append(
                f"  track {p['track']:03d} {p['direction']}: last "
                f"{p['last_utc']}, next ~{p['next_utc']}"
            )
    else:
        lines.append("")
        lines.append("No NISAR coverage found at the epicentre yet.")

    prov = solution["provenance"]
    lines += [
        "",
        "--",
        f"Automated solution (mttime {prov['mttime_version']}, "
        f"model {prov['velocity_model']}, GF {prov['gf_version']}). "
        f"Not reviewed by a human. Solutions archived at the auto_tdmt_NZ repo.",
    ]
    return subject, "\n".join(lines)


def send_email(subject: str, body: str, attachments: list[Path]) -> None:
    host = os.environ.get("SMTP_HOST")
    assert host, "SMTP_HOST not set — cannot send email"
    port = int(os.environ.get("SMTP_PORT", "587"))
    user = os.environ["SMTP_USER"]
    password = os.environ["SMTP_PASS"]
    mail_from = os.environ.get("MAIL_FROM", user)
    mail_to = os.environ["MAIL_TO"]

    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = mail_from
    msg["To"] = mail_to
    msg.set_content(body)
    for path in attachments:
        assert path.exists(), f"attachment missing: {path}"
        msg.add_attachment(
            path.read_bytes(), maintype="image", subtype="jpeg",
            filename=path.name,
        )

    with smtplib.SMTP(host, port, timeout=60) as s:
        s.starttls()
        s.login(user, password)
        s.send_message(msg)
    print(f"emailed '{subject}' to {mail_to} "
          f"({len(attachments)} attachments)")
