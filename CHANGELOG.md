# Changelog

## Version 2 — October 2026 (manuscript v4)

### Correction
- **All-births trend comparator (Fig 6 / figF).** ONS Table 10 (live births by mother's age,
  1938–2025) reports its oldest maternal age group as **"40 and over"**. Version 1 mapped only
  "40 to 44" and "45 and over", so every birth to a mother aged 40+ was silently dropped from
  the year-by-year all-births distribution. The trend overlap was then computed over ages
  15–39 only and was numerically identical to the share of uterine cancers diagnosed under 40.
  - v1 (incorrect): overlap 0.012 (2001) → 0.024 (2023), "+102%".
  - v2 (corrected): overlap **0.031 (2001) → 0.044 (2023), +41%**.
  - The headline 2023 overlap (0.048 with first births, 0.049 with all births; England 2025
    tables, which do separate 40–44 and 45+) was not affected.
  - Births at 40+ are now placed in the 40–44 band. Births at 45+ are <1% of births, so the
    corrected trend is, if anything, slightly conservative.

### Additions
- Annual percentage change (APC) in registrations and in age-specific rates below 40, 45 and
  50 (log-linear quasi-Poisson regression, 2001–2023), and APC of the all-ages
  age-standardised rate. Female population denominators are recovered from the NDRS
  all-cancer age-specific rates (count × 100,000 / rate); the 2012 rates are blank in the
  export, so the 2012 population is interpolated geometrically from 2011 and 2013.
- Under-40 registrations series (added to figG) and shares under 40/45/50 for 2001 and 2023.
- Decomposition of the headline overlap into the cancer-limited part (share of cancers
  diagnosed 15–44) and the birth-limited part (share of first births at 45+).
- Projected births to women aged 40+ (2023, 2047), projected all-ages births, and a check
  that the ONS projection allocates all births to maternal ages 15–44.
- Sensitivity analysis restricted to ICD-10 C54 (corpus uteri), and the C55 share over time.
- Share of England & Wales births to mothers aged 40+ in 2001 and 2023.
- The England & Wales standardised mean age for all births in 2025 (31.1; the manuscript
  previously quoted 31.0, which was the England banded mean) and the 2018 under-40 count
  (157), so that every number quoted in the manuscript appears in `results_all.csv`.

### Presentation
- Figure labels use "uterine cancer (C54–C55)" for registry data and state the staging system
  (TNM 8th edition / FIGO 2009) in figD.
- Figure 6 (figF): y-axis headroom so the legend no longer overlaps the data.
- Figure font pinned to DejaVu Sans (bundled with matplotlib, and the font the v1 figures
  were drawn in) so figures render identically across environments. Figures whose content
  did not change (A, B, C, H, I) differ from v1 only by sub-pixel layout differences between
  matplotlib versions.
- Stale file names in the script header corrected.
- The notebook is regenerated from the script (cell-for-cell identical code).
