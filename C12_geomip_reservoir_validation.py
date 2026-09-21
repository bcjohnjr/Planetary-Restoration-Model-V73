#!/usr/bin/env python3
"""V73 GeoMIP-constrained aerosol-reservoir validation.

Purpose
-------
Stress-test the V72 feedback controller against aerosol persistence absent from
its simple one-year-lag actuator.  The controller is NOT retuned.  A first-order
aerosol/forcing reservoir is inserted between the delayed SAI command and the
two-layer temperature-response plant.  Reservoir e-folding times span the
approximate 8-17 month sulfur-lifetime range reported across four Earth system
models in G6-1.5K-SAI (Lee et al., Atmos. Chem. Phys. 26, 7463-7483, 2026;
doi:10.5194/acp-26-7463-2026).

This is a reduced-complexity robustness benchmark, not a coupled-Earth-system
simulation and not an operational SAI design.  The SO2-injection and
precipitation diagnostics are literature-response envelopes only.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
BASE = pd.read_csv(ROOT / "D20_controller.csv")
UNC = pd.read_csv(ROOT / "D13_ctrl_uncert.csv")
BENCH = pd.read_csv(ROOT / "D31_geomip_benchmark.csv")

YEARS = BASE.year.to_numpy(int)
CO2 = BASE.co2_ppm.to_numpy(float)
TBASE = BASE.baseline_temperature_c.to_numpy(float)
CAP = BASE.temperature_cap_c.to_numpy(float)
DEPLOY_START = 2035
KP = 5.0
KI = 0.12
MAX_SAI_CMD = 3.5
MAX_DCMD = 0.20
PHASEOUT_START_PPM = 300.0
PHASEOUT_END_PPM = 280.0
TAUS = [0.67, 1.0, 1.5]
FORCING_ZERO_THRESHOLD = 0.01  # W m-2 magnitude
RNG = np.random.default_rng(7301)


def phaseout_limit(co2_ppm: float) -> float:
    if co2_ppm >= PHASEOUT_START_PPM:
        return MAX_SAI_CMD
    return MAX_SAI_CMD * np.clip(
        (co2_ppm - PHASEOUT_END_PPM) / (PHASEOUT_START_PPM - PHASEOUT_END_PPM),
        0.0, 1.0,
    )


def sai_equilibrium_forcing(command_mag: float, efficacy: float, alpha: float) -> float:
    """Equilibrium cooling-forcing magnitude for a sustained command."""
    c = max(0.0, float(command_mag))
    return efficacy * c / (1.0 + alpha * c)


def simulate_reservoir(
    tau_years: float,
    C1: float = 7.0,
    C2: float = 100.0,
    lamb: float = 1.2,
    gamma: float = 0.7,
    efficacy: float = 1.0,
    alpha: float = 0.12,
    logistics_lag_years: int = 1,
    noise: np.ndarray | None = None,
) -> pd.DataFrame:
    """Run unchanged V72 PI controller with an added aerosol reservoir.

    A command chosen after observing year t becomes available after the V72
    logistics lag.  Instead of becoming instantaneous effective forcing, it is
    treated as the equilibrium forcing of a first-order aerosol reservoir:

        A[t+1] = exp(-dt/tau) A[t] + (1-exp(-dt/tau)) F_eq(command_in)

    The physical radiative forcing applied to the plant is -A.
    """
    if tau_years <= 0:
        raise ValueError("tau_years must be positive")
    n = len(YEARS)
    temp = np.empty(n, float)
    cmd = np.zeros(n, float)
    reservoir = np.zeros(n, float)
    realized = np.zeros(n, float)
    d1 = d2 = integ = prev = 0.0
    queue = [0.0] * max(1, int(logistics_lag_years))
    decay = math.exp(-1.0 / float(tau_years))
    A = 0.0

    for i, year in enumerate(YEARS):
        incoming = queue.pop(0)
        eq_mag = sai_equilibrium_forcing(incoming, efficacy, alpha)
        A = decay * A + (1.0 - decay) * eq_mag
        reservoir[i] = A
        realized[i] = -A

        rf = -A
        nd1 = d1 + (rf - lamb * d1 - gamma * (d1 - d2)) / C1
        nd2 = d2 + gamma * (d1 - d2) / C2
        d1, d2 = nd1, nd2
        temp[i] = TBASE[i] + d1 + (0.0 if noise is None else float(noise[i]))

        if year < DEPLOY_START:
            new = 0.0
        else:
            cap = CAP[i]
            err = 0.0 if np.isnan(cap) else temp[i] - cap
            if err > 0:
                integ = min(5.0, integ + err)
            else:
                integ = max(0.0, integ + 0.25 * err)
            desired = max(0.0, KP * max(0.0, err) + KI * integ)
            lim = phaseout_limit(CO2[i])
            desired = min(desired, lim, MAX_SAI_CMD)
            new = float(np.clip(
                desired,
                max(0.0, prev - MAX_DCMD),
                min(MAX_SAI_CMD, prev + MAX_DCMD),
            ))
            new = min(new, lim)
        cmd[i] = new
        queue.append(new)
        prev = new

    return pd.DataFrame({
        "year": YEARS,
        "co2_ppm": CO2,
        "baseline_temperature_c": TBASE,
        "temperature_cap_c": CAP,
        "controlled_temperature_c": temp,
        "sai_command_magnitude_wm2": cmd,
        "aerosol_reservoir_magnitude_wm2": reservoir,
        "sai_realized_forcing_wm2": realized,
    })


def first_year_at_or_below(x: np.ndarray, threshold: float, start: int = DEPLOY_START):
    m = (YEARS >= start) & (x <= threshold)
    idx = np.where(m)[0]
    return None if len(idx) == 0 else int(YEARS[idx[0]])


def first_year_after_last_active(x: np.ndarray, threshold: float):
    idx = np.where(np.abs(x) > threshold)[0]
    if len(idx) == 0:
        return int(YEARS[0])
    j = int(idx[-1]) + 1
    return None if j >= len(YEARS) else int(YEARS[j])


def window_slope(y: np.ndarray, start_year, years: int = 20):
    if start_year is None or start_year + years > int(YEARS[-1]):
        return np.nan
    m = (YEARS >= start_year) & (YEARS <= start_year + years)
    if m.sum() < 2:
        return np.nan
    xx = YEARS[m].astype(float)
    yy = y[m].astype(float)
    return float(np.polyfit(xx - xx[0], yy, 1)[0])


def metrics(df: pd.DataFrame) -> dict:
    t = df.controlled_temperature_c.to_numpy(float)
    cmd = df.sai_command_magnitude_wm2.to_numpy(float)
    rf = df.sai_realized_forcing_wm2.to_numpy(float)
    mask = YEARS >= DEPLOY_START
    cap = CAP[mask]
    exc = np.maximum(0.0, t[mask] - cap)
    cmd_zero = first_year_after_last_active(cmd, 1e-6)
    forcing_zero = first_year_after_last_active(rf, FORCING_ZERO_THRESHOLD)
    return {
        "T2060_c": float(t[YEARS == 2060][0]),
        "T2100_c": float(t[YEARS == 2100][0]),
        "T2200_c": float(t[YEARS == 2200][0]),
        "T2300_c": float(t[YEARS == 2300][0]),
        "peak_command_wm2": float(cmd.max()),
        "peak_reservoir_wm2": float((-rf).max()),
        "max_cap_exceedance_c": float(exc.max()),
        "mean_cap_exceedance_c": float(exc.mean()),
        "first_le_1C_year": first_year_at_or_below(t, 1.0),
        "first_command_zero_after_last_active_year": cmd_zero,
        "first_physical_forcing_below_0p01_after_last_active_year": forcing_zero,
        "post_physical_zero_20y_ols_trend_c_per_yr": window_slope(t, forcing_zero, 20),
        "max_annual_warming_after_deploy_c_per_yr": float(np.max(np.diff(t[YEARS >= DEPLOY_START]))),
    }


def quantile_row(label: str, tau: float, rows: pd.DataFrame) -> dict:
    out = {"case": label, "tau_years": tau, "n": int(len(rows))}
    cols = [
        "T2060_c", "T2100_c", "T2200_c", "T2300_c",
        "peak_command_wm2", "peak_reservoir_wm2",
        "max_cap_exceedance_c", "mean_cap_exceedance_c",
        "first_le_1C_year",
        "first_command_zero_after_last_active_year",
        "first_physical_forcing_below_0p01_after_last_active_year",
        "post_physical_zero_20y_ols_trend_c_per_yr",
        "max_annual_warming_after_deploy_c_per_yr",
    ]
    for c in cols:
        vals = pd.to_numeric(rows[c], errors="coerce").dropna().to_numpy(float)
        if len(vals):
            out[f"{c}_p05"] = float(np.quantile(vals, 0.05))
            out[f"{c}_p50"] = float(np.quantile(vals, 0.50))
            out[f"{c}_p95"] = float(np.quantile(vals, 0.95))
    trends = pd.to_numeric(rows.post_physical_zero_20y_ols_trend_c_per_yr, errors="coerce").dropna().to_numpy(float)
    out["post_zero_20y_valid_n"] = int(len(trends))
    out["post_zero_20y_unclassified_n"] = int(len(rows) - len(trends))
    out["fraction_post_zero_warming_trend_gt_0_among_valid"] = float(np.mean(trends > 0)) if len(trends) else float("nan")
    out["fraction_peak_command_hits_3p5_cap"] = float(np.mean(rows.peak_command_wm2.to_numpy(float) >= 3.499))
    return out


# Central deterministic trajectories at GeoMIP-bounded aerosol lifetimes.
central_frames = []
central_metrics = []
for tau in TAUS:
    z = simulate_reservoir(tau)
    z.insert(0, "tau_years", tau)
    central_frames.append(z)
    mm = metrics(z)
    mm.update({"case": "central", "tau_years": tau, "member": "central"})
    central_metrics.append(mm)
central_ts = pd.concat(central_frames, ignore_index=True)
central_ts.to_csv(ROOT / "D32_geomip_reservoir_central_timeseries.csv", index=False)

# Re-use V72's 500 sampled climate plants exactly.  We intentionally omit the
# AR(1) observational/noise term here so actuator persistence is isolated.
ensemble_rows = []
for tau in TAUS:
    for _, r in UNC.iterrows():
        z = simulate_reservoir(
            tau_years=tau,
            C1=float(r.C1), C2=float(r.C2), lamb=float(r.lambda_wm2k), gamma=float(r.gamma),
            efficacy=float(r.sai_efficacy), alpha=float(r.sai_nonlinearity_alpha),
            logistics_lag_years=int(r.lag_years), noise=None,
        )
        mm = metrics(z)
        mm.update({
            "case": "V72_plant_ensemble_noise_free",
            "tau_years": tau,
            "member": int(r.member),
            "C1": float(r.C1), "C2": float(r.C2),
            "lambda_wm2k": float(r.lambda_wm2k), "gamma": float(r.gamma),
            "sai_efficacy": float(r.sai_efficacy),
            "sai_nonlinearity_alpha": float(r.sai_nonlinearity_alpha),
            "logistics_lag_years": int(r.lag_years),
        })
        ensemble_rows.append(mm)
ensemble = pd.DataFrame(ensemble_rows)
ensemble.to_csv(ROOT / "D33_geomip_reservoir_ensemble.csv", index=False)

summary_rows = []
for tau in TAUS:
    sub = ensemble[ensemble.tau_years == tau].copy()
    summary_rows.append(quantile_row("V72_plant_ensemble_noise_free", tau, sub))
summary = pd.DataFrame(summary_rows)
summary.to_csv(ROOT / "D34_geomip_reservoir_summary.csv", index=False)

# External response-envelope diagnostics.  These are NOT a mass-forcing
# conversion and should not be used to design injection.  They ask whether the
# amount of global cooling produced by the central reduced model lies in a
# coupled-ESM response scale reported in GeoMIP.
cool_center = 0.11
cool_low = 0.10
cool_high = 0.12
precip_center_per_c = -0.11 / 1.4
precip_low_per_c = -0.15 / 1.4
precip_high_per_c = -0.07 / 1.4
proxy_rows = []
for tau in TAUS:
    z = central_ts[central_ts.tau_years == tau]
    cooling = np.maximum(0.0, z.baseline_temperature_c.to_numpy(float) - z.controlled_temperature_c.to_numpy(float))
    for year, c in zip(z.year.to_numpy(int), cooling):
        proxy_rows.append({
            "tau_years": tau,
            "year": int(year),
            "global_cooling_c": float(c),
            "geomip_implied_injection_center_tgso2yr": float(c / cool_center),
            "geomip_implied_injection_low_tgso2yr": float(c / cool_high),
            "geomip_implied_injection_high_tgso2yr": float(c / cool_low),
            "geomip_precip_response_scale_center_mmday": float(c * precip_center_per_c),
            "geomip_precip_response_scale_low_mmday": float(c * precip_low_per_c),
            "geomip_precip_response_scale_high_mmday": float(c * precip_high_per_c),
        })
proxy = pd.DataFrame(proxy_rows)
proxy.to_csv(ROOT / "D35_geomip_response_envelope.csv", index=False)

# Compact machine-readable summary.
out = {
    "validation": "GeoMIP-constrained aerosol-reservoir robustness test",
    "controller_retuned": False,
    "controller": {"Kp": KP, "Ki": KI, "max_command_wm2": MAX_SAI_CMD, "rate_limit_wm2_per_yr": MAX_DCMD},
    "geomip_source": "Lee et al. 2026, Atmos. Chem. Phys. 26, 7463-7483",
    "geomip_doi": "10.5194/acp-26-7463-2026",
    "aerosol_tau_years_tested": TAUS,
    "ensemble_n_per_tau": int(len(UNC)),
    "central": central_metrics,
    "ensemble_summary": summary.to_dict(orient="records"),
}
(ROOT / "V73_VALIDATION_SUMMARY.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

# Human-readable results.
lines = [
    "# V73 GeoMIP-constrained aerosol-reservoir validation",
    "",
    "This test inserts a first-order aerosol/forcing reservoir into the V72 actuator and leaves the V72 controller gains unchanged.",
    "The tested 0.67, 1.0 and 1.5 yr e-folding times span the approximately 8-17 month sulfur lifetimes reported across four G6-1.5K-SAI Earth system models (Lee et al., 2026; doi:10.5194/acp-26-7463-2026).",
    "The 500-member tests re-use V72's sampled two-layer plant parameters exactly and omit AR(1) noise so the incremental effect of aerosol persistence is isolated.",
    "",
    "## Central trajectories",
]
for m in central_metrics:
    lines += [
        f"- tau={m['tau_years']:.2f} yr: T2060={m['T2060_c']:.3f} C; T2100={m['T2100_c']:.3f} C; peak command={m['peak_command_wm2']:.3f} W m-2; max cap exceedance={m['max_cap_exceedance_c']:.3f} C; command zero after last active={m['first_command_zero_after_last_active_year']}; physical forcing <0.01 W m-2 after last active={m['first_physical_forcing_below_0p01_after_last_active_year']}; post-zero 20-y OLS trend={m['post_physical_zero_20y_ols_trend_c_per_yr']:.5f} C yr-1."
    ]
lines += ["", "## 500-member V72 plant ensemble"]
for _, r in summary.iterrows():
    lines += [
        f"- tau={r.tau_years:.2f} yr: median T2060={r.T2060_c_p50:.3f} C (5-95% {r.T2060_c_p05:.3f}-{r.T2060_c_p95:.3f}); median max cap exceedance={r.max_cap_exceedance_c_p50:.3f} C; median peak command={r.peak_command_wm2_p50:.3f} W m-2; cap-hit fraction={r.fraction_peak_command_hits_3p5_cap:.3f}; complete 20-y post-zero windows={int(r.post_zero_20y_valid_n)}/{int(r.n)}; positive-trend fraction among complete windows={r.fraction_post_zero_warming_trend_gt_0_among_valid:.3f}."
    ]
lines += [
    "",
    "## External response-envelope diagnostics",
    "GeoMIP reports 0.11 +/- 0.01 C global cooling per Tg SO2 yr-1 and an approximately 0.11 +/- 0.04 mm day-1 global precipitation reduction while cooling approximately 1.4 C. We use these only to scale the central reduced-model cooling into a literature-response envelope. They are not a conversion from W m-2 to SO2 and are not a regional-impact prediction.",
]
for tau in TAUS:
    p = proxy[proxy.tau_years == tau]
    imax = p.geomip_implied_injection_center_tgso2yr.idxmax()
    rr = p.loc[imax]
    lines += [f"- tau={tau:.2f} yr: maximum cooling-equivalent injection scale={rr.geomip_implied_injection_center_tgso2yr:.2f} Tg SO2 yr-1 in {int(rr.year)}; corresponding global precipitation-response scale={rr.geomip_precip_response_scale_center_mmday:.3f} mm day-1."]
lines += [
    "",
    "## Interpretation boundary",
    "Passing this validation means the phaseout result is not an artifact of assuming zero aerosol persistence over the published GeoMIP sulfur-lifetime range. It does not validate regional precipitation, ozone chemistry, aerosol microphysics, injection logistics, governance, or safety; those require coupled Earth system experiments.",
]
(ROOT / "V73_VALIDATION_RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

# One compact validation figure: central trajectories, ensemble 2060 response,
# command requirement, and literature response-envelope scale.
fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))
ax = axes[0, 0]
ax.plot(YEARS, TBASE, label="No-SAI baseline", linewidth=1.6)
for tau in TAUS:
    z = central_ts[central_ts.tau_years == tau]
    ax.plot(z.year, z.controlled_temperature_c, label=f"tau={tau:g} yr")
ax.plot(YEARS, CAP, linestyle="--", linewidth=1.0, label="Temperature cap")
ax.set_xlim(2030, 2120); ax.set_ylabel("Global temperature anomaly (C)"); ax.set_title("a  Central controlled temperature")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

ax = axes[0, 1]
labels=[]; meds=[]; lo=[]; hi=[]
for tau in TAUS:
    s=summary[summary.tau_years==tau].iloc[0]
    labels.append(str(tau)); meds.append(s.T2060_c_p50); lo.append(s.T2060_c_p05); hi.append(s.T2060_c_p95)
x=np.arange(len(TAUS)); meds=np.array(meds); lo=np.array(lo); hi=np.array(hi)
ax.errorbar(x, meds, yerr=np.vstack([meds-lo, hi-meds]), fmt='o', capsize=4)
ax.set_xticks(x, labels); ax.set_xlabel("Aerosol e-folding time (yr)"); ax.set_ylabel("2060 temperature (C)"); ax.set_title("b  V72 plant uncertainty (5-95%)"); ax.grid(alpha=0.25)

ax = axes[1, 0]
for tau in TAUS:
    z = central_ts[central_ts.tau_years == tau]
    ax.plot(z.year, z.sai_command_magnitude_wm2, label=f"tau={tau:g} yr")
    ax.plot(z.year, z.aerosol_reservoir_magnitude_wm2, linestyle=':', linewidth=1.0)
ax.set_xlim(2030, 2300); ax.set_ylabel("Magnitude (W m-2)"); ax.set_xlabel("Year"); ax.set_title("c  Command (solid) and aerosol reservoir (dotted)"); ax.grid(alpha=0.25)

ax = axes[1, 1]
for tau in TAUS:
    p=proxy[proxy.tau_years==tau]
    ax.plot(p.year, p.geomip_implied_injection_center_tgso2yr, label=f"tau={tau:g} yr")
ax.axhspan(10, 20, alpha=0.12, label="10-20 Tg yr-1 nonlinearity context")
ax.set_xlim(2030, 2120); ax.set_ylabel("Cooling-equivalent SO2 scale (Tg yr-1)"); ax.set_xlabel("Year"); ax.set_title("d  GeoMIP response-envelope diagnostic"); ax.legend(fontsize=8); ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(ROOT / "F7_geomip_reservoir_validation.png", dpi=220, bbox_inches="tight")
plt.close(fig)

print((ROOT / "V73_VALIDATION_RESULTS.md").read_text(encoding="utf-8"))
