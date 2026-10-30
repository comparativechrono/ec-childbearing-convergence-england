# Endometrial cancer and delayed childbearing in England: a convergence analysis

Reproducible analysis accompanying the manuscript *The collision of delayed childbearing and
endometrial cancer in England: a population-level convergence analysis, 2001–2047* (in
preparation). The same material is supplied with the manuscript as Supplementary material 1.

Starting from raw public data, the pipeline quantifies the convergence between the age
distribution of uterine cancer diagnoses (ICD-10 C54–C55; 96.6% endometrial in 2023) and the
age distribution of childbearing in England, 2001–2023, with a forward projection to 2047.

Repository: https://github.com/comparativechrono/ec-childbearing-convergence-england

**Version 2 (October 2026).** Corrects an error in version 1's all-births trend comparator
and adds trend statistics and sensitivity analyses — see `CHANGELOG.md`. All other v1
figures and numbers are unchanged.

## How to run

1. Install dependencies (Python 3.9+):
   ```
   pip install -r requirements.txt
   ```
2. The five raw files are included in `data/` (see "Data files" below).
3. Run either form — both produce identical outputs in `outputs/`:
   - Script:   `python convergence_analysis.py`
   - Notebook: open `ec_childbearing_convergence.ipynb` and run all cells.

Outputs (9 figures + `results_all.csv`) are written to `outputs/`. A pre-generated copy of
`outputs/` is included so you can confirm your run matches. Version 2 outputs were generated
with Python 3.13, pandas 3.0, NumPy 2.5, Matplotlib 3.11, SciPy 1.18 and openpyxl 3.1;
`requirements.txt` gives minimum versions.

## Data files (in `data/`)

| file | description | source |
|---|---|---|
| `Incidence_National_ICD10_byage_singleyear.csv` | Uterine cancer incidence by ICD-10, sex, 5-year age band, single years, England, 2001–2023 (age-distribution, overlap and trend work; the all-cancer rows supply the female population denominators) | NDRS CancerData dashboard export — https://digital.nhs.uk/ndrs/data/data-outputs/cancerdata |
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
- **figD** endometrial vs non-endometrial composition and stage (TNM 8 / FIGO 2009), 2023
- **figE** headline overlap: cancer diagnosis age vs first-birth age (shaded overlap)
- **figF** overlap trend 2001–2023 and % of uterine cancers diagnosed under 50
- **figG** reproductive-age registrations (under 40 / 45 / 50), 2001–2023
- **figH** forward projection: mean maternal age and % births to women 35+, to 2047
- **figI** forward projection: number of births to women 35+ / 40+, to 2047
- **results_all.csv** every reported statistic

## Key results (England)

- Uterine cancer age-standardised incidence rose ~32% (3-year rolling rates, 2001–03 →
  2021–23; 1.3%/yr on single-year rates), whereas incidence below 45 rose 3.4%/yr (95% CI
  2.7–4.1) and below 40 4.8%/yr (3.8–5.7). Registrations under 40 went from 61 (2001) to 219
  (2023); under 45 from 158 to 397; under 50 from 333 to 696.
- 96.6% of uterine cancers are in the NDRS endometrial group; 61% of endometrial cancers are
  stage 1 (72% of those with known stage; TNM 8, equivalent to FIGO 2009).
- Mean age at first birth (England & Wales) rose from a 1969 trough of 23.7 to 29.6 in 2025;
  the England vs England-&-Wales gap is 0.0 years (validates the E&W longitudinal comparator).
- Cancer–first-birth overlap coefficient ≈ 0.05 (2023). The all-births overlap rose by 41%
  (0.031 → 0.044) over 2001–2023, all of the net rise after 2016. In practice the coefficient
  is the share of cancers diagnosed before 45 plus the small share of births at 45+.
- Projection: mean maternal age 30.7→32.1 and births to women 35+ rising from ~25% to 34%
  of all births by 2047 (35+: ~144,000 → 205,000 a year; 40+: ~27,000 → 58,000). Because the
  overlap is set by cancer incidence below 45, it does not respond to birth-side change alone
  and the projection is reported through these demographic quantities (see `methods.md`).
- Sensitivity: restricting to corpus uteri (C54) gives the same 2023 overlap (0.048) and a
  slightly steeper rise in younger registrations.

## Files

```
data/                              raw inputs
ec_childbearing_convergence.ipynb  notebook (run all cells)
convergence_analysis.py            identical pipeline as a script
methods.md                         full Methods description
CHANGELOG.md                       what changed between versions
README.md                          this file
requirements.txt                   Python dependencies
outputs/                           pre-generated figures + results_all.csv
```

## Data licence

The raw files in `data/` are public sector information published by the Office for National
Statistics and by NHS England (National Disease Registration Service), reused under the
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Citation

If you use this code or its outputs, please cite the associated article; full citation
details will be added here on publication.
