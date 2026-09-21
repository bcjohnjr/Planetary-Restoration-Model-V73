V73 GEOMIP-CONSTRAINED AEROSOL-RESERVOIR VALIDATION

Purpose
-------
This validation leaves the V72 PI controller unchanged and inserts a persistent
first-order aerosol/forcing reservoir downstream of the existing logistics lag.
It tests e-folding times of 0.67, 1.0 and 1.5 yr, spanning the approximately
8-17 month sulfur-lifetime range reported across four G6-1.5K-SAI Earth-system
models by Lee et al. (2026), Atmos. Chem. Phys. 26, 7463-7483,
doi:10.5194/acp-26-7463-2026.

Run locally
-----------
python -m pip install -r requirements-v73.txt
python C12_geomip_reservoir_validation.py

Run on GitHub
-------------
Keep C12_geomip_reservoir_validation.py, D20_controller.csv,
D13_ctrl_uncert.csv and D31_geomip_benchmark.csv in the repository root and
copy .github/workflows/v73_geomip_validation.yml to the same repository.
Commit/push and run the workflow from GitHub Actions.

Prespecified interpretation
---------------------------
A pass means the managed-phaseout result is not an artifact of assuming zero
aerosol persistence over the tested GeoMIP-bounded lifetime range. It does not
validate regional precipitation, ozone chemistry, aerosol microphysics,
injection logistics, governance or safety. Coupled Earth-system experiments
are required for those questions.

Post-phaseout trend rule
------------------------
A 20-year OLS temperature trend is classified only if the complete 20-year
window fits within the simulation ending in 2300. Late phaseouts are explicitly
unclassified; they are not counted as stable.
