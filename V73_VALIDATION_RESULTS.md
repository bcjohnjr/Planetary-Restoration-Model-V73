# V73 GeoMIP-constrained aerosol-reservoir validation

This test inserts a first-order aerosol/forcing reservoir into the V72 actuator and leaves the V72 controller gains unchanged.
The tested 0.67, 1.0 and 1.5 yr e-folding times span the approximately 8-17 month sulfur lifetimes reported across four G6-1.5K-SAI Earth system models (Lee et al., 2026; doi:10.5194/acp-26-7463-2026).
The 500-member tests re-use V72's sampled two-layer plant parameters exactly and omit AR(1) noise so the incremental effect of aerosol persistence is isolated.

## Central trajectories
- tau=0.67 yr: T2060=1.164 C; T2100=1.064 C; peak command=1.147 W m-2; max cap exceedance=0.174 C; command zero after last active=2263; physical forcing <0.01 W m-2 after last active=2263; post-zero 20-y OLS trend=-0.00321 C yr-1.
- tau=1.00 yr: T2060=1.164 C; T2100=1.062 C; peak command=1.150 W m-2; max cap exceedance=0.175 C; command zero after last active=2262; physical forcing <0.01 W m-2 after last active=2263; post-zero 20-y OLS trend=-0.00319 C yr-1.
- tau=1.50 yr: T2060=1.166 C; T2100=1.058 C; peak command=1.164 W m-2; max cap exceedance=0.176 C; command zero after last active=2262; physical forcing <0.01 W m-2 after last active=2263; post-zero 20-y OLS trend=-0.00315 C yr-1.

## 500-member V72 plant ensemble
- tau=0.67 yr: median T2060=1.200 C (5-95% 1.149-1.248); median max cap exceedance=0.207 C; median peak command=1.415 W m-2; cap-hit fraction=0.000; complete 20-y post-zero windows=399/500; positive-trend fraction among complete windows=0.000.
- tau=1.00 yr: median T2060=1.206 C (5-95% 1.151-1.255); median max cap exceedance=0.208 C; median peak command=1.427 W m-2; cap-hit fraction=0.000; complete 20-y post-zero windows=398/500; positive-trend fraction among complete windows=0.000.
- tau=1.50 yr: median T2060=1.206 C (5-95% 1.155-1.255); median max cap exceedance=0.208 C; median peak command=1.433 W m-2; cap-hit fraction=0.000; complete 20-y post-zero windows=395/500; positive-trend fraction among complete windows=0.000.

## External response-envelope diagnostics
GeoMIP reports 0.11 +/- 0.01 C global cooling per Tg SO2 yr-1 and an approximately 0.11 +/- 0.04 mm day-1 global precipitation reduction while cooling approximately 1.4 C. We use these only to scale the central reduced-model cooling into a literature-response envelope. They are not a conversion from W m-2 to SO2 and are not a regional-impact prediction.
- tau=0.67 yr: maximum cooling-equivalent injection scale=4.62 Tg SO2 yr-1 in 2074; corresponding global precipitation-response scale=-0.040 mm day-1.
- tau=1.00 yr: maximum cooling-equivalent injection scale=4.63 Tg SO2 yr-1 in 2074; corresponding global precipitation-response scale=-0.040 mm day-1.
- tau=1.50 yr: maximum cooling-equivalent injection scale=4.64 Tg SO2 yr-1 in 2074; corresponding global precipitation-response scale=-0.040 mm day-1.

## Interpretation boundary
Passing this validation means the phaseout result is not an artifact of assuming zero aerosol persistence over the published GeoMIP sulfur-lifetime range. It does not validate regional precipitation, ozone chemistry, aerosol microphysics, injection logistics, governance, or safety; those require coupled Earth system experiments.
