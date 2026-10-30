# Endometrial cancer & delayed childbearing in England — reproducible analysis

This bundle reproduces, from raw public data, the convergence analysis between the age
distribution of endometrial/uterine cancer diagnoses and the age distribution of
childbearing in England, including a forward projection to 2047.

## How to run

1. Install dependencies (Python 3.9+):
   ```
   pip install -r requirements.txt
   ```
2. Ensure the four raw files are in `data/` (they are included; see "Data files" below).
3. Run either form — both produce identical outputs in `outputs/`:
   - Notebook: open `ec_childbearing_convergence.ipynb` and run all cells.
   - Script:   `python convergence_analysis.py`

Outputs (9 figures + `results_all.csv`) are written to `outputs/`. A pre-generated copy
of `outputs/` is included so you can confirm your run matches.

## Data files (in `data/`)

| file | description | source |
|---|---|---|
| `Incidence_National_ICD10_byage_singleyear.csv` | Uterine cancer incidence by ICD-10, sex, 5-year age band, single years, England, 2001–2023 (age-distribution and overlap work) | NDRS CancerData dashboard export — https://digital.nhs.uk/ndrs/data/data-outputs/cancerdata |
| `Incidence_National_ICD10_allages_rolling3yr.csv` | Uterine cancer incidence, all ages, NDRS 3-year rolling age-standardised rates, England, 2001–2023 (incidence-trend figure; gap-free) | NDRS CancerData dashboard export — https://digital.nhs.uk/ndrs/data/data-outputs/cancerdata |
| `Cancer_registrations_2023_machine_readable_files.zip` | NHS cancer registration statistics, England 2023 (Table 1 = endometrial vs non-endometrial × stage) | NHS Digital / NDRS — https://files.digital.nhs.uk/16/5B8561/Cancer_registrations_2023_machine_readable_files.zip |
| `2025birthregistrations.xlsx` | ONS Births in England and Wales: birth registrations, 2025 edition (Tables 2, 4, 10, 11) | ONS — birthsinenglandandwalesbirthregistrations/2025 |
| `2022_SNPP_Births_persons.csv` | ONS 2022-based subnational projections, births by age of mother (all babies) | ONS — birthsbyageofmotherz3 (2022 SNPP, migration-category variant) |

To refresh with newer editions, replace the file(s) in `data/` with the same column
layout; the loader auto-detects ONS header rows.

## What the pipeline produces

- **figA** maternal age distribution (all vs first births), England 2025
- **figB** mean age at first birth, England & Wales, 1938–2025
- **figC** uterine cancer (C54–C55) incidence trend, England, 2001–2023
- **figD** endometrial vs non-endometrial composition and stage, 2023
- **figE** headline overlap: cancer diagnosis age vs first-birth age (shaded overlap)
- **figF** overlap trend 2001–2023 and % of EC diagnosed under 50
- **figG** absolute reproductive-age EC registrations (under 45 / under 50), 2001–2023
- **figH** forward projection: mean maternal age and % births to women 35+, to 2047
- **figI** forward projection: number of births to women 35+ / 40+, to 2047
- **results_all.csv** every reported statistic

## Key results (England)

- Uterine cancer ASR rose ~32% (2002→2022); reproductive-age burden rose more in counts:
  registrations under 45 went from 158 (2001) to 397 (2023).
- ~97% of uterine cancers are endometrial; ~61% of endometrial cancers are stage 1.
- Mean age at first birth rose from a 1969 trough of 23.7 to 29.6 in 2025; the England vs
  England-&-Wales gap is 0.0 years (validates the E&W longitudinal comparator).
- Cancer–first-birth overlap coefficient ≈ 0.05 (2023); the all-births overlap roughly
  doubled (+102%) over 2001–2023.
- Projection: mean maternal age 30.7→32.1 and births to women 35+ rising from ~25% to 34%
  of all births by 2047. The proportion-normalised overlap is pinned by EC incidence in the
  reproductive zone and is therefore reported via these demographic quantities (see Methods).

## Files

```
data/                              raw inputs
ec_childbearing_convergence.ipynb  notebook (run all cells)
convergence_analysis.py            identical pipeline as a script
methods.md                         full Methods section
README.md                          this file
requirements.txt                   Python dependencies
outputs/                           pre-generated figures + results_all.csv
```
