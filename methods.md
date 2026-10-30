# Methods

## Study design and setting

We conducted a population-level descriptive and projection study for England, quantifying the demographic convergence between the age at which women are diagnosed with endometrial (uterine) cancer and the age at which women give birth. The analysis uses routinely collected, publicly available national data and involves no individual patient contact. England was chosen so that cancer registration (national, via the National Disease Registration Service) and the cancer denominators align; for the long-run childbearing trend, where England-only series are not published historically, England and Wales data were used and validated against England-only figures for the most recent year (see "Geography and validation").

## Data sources

**Endometrial/uterine cancer incidence by age.** Female cancer registrations for England, 2001–2023, were obtained from the National Disease Registration Service (NDRS) CancerData incidence tool, disaggregated by single diagnosis year, sex, five-year age band, and three-character ICD-10 site. The primary disease definition was the official uterine-cancer grouping ICD-10 **C54–C55** (malignant neoplasm of the corpus uteri, and uterus part unspecified); corpus uteri alone (C54) was retained for sensitivity. Counts by age band reconcile exactly to the published all-ages totals (for example, 9,048 female C54–C55 registrations in 2023). The incidence-trend figure uses the NDRS three-year rolling age-standardised rates (which are gap-free), whereas all age-resolved analyses use single-year registrations; the single-year all-ages rate for 2012 is absent from the export (the count is present), which is why the rolling series is used for the trend.

**Endometrial subtype and stage.** Because the age-banded export is provided for all stages combined, the endometrial-versus-other composition and the stage distribution were taken from the NDRS machine-readable cancer-registration tables for 2023 (Table 1, incidence by NDRS group), which separates "Uterus – endometrial" from "Uterus – non-endometrial" and gives stage at diagnosis. These data are all-ages.

**Childbearing.** Birth data were taken from the Office for National Statistics (ONS) "Births in England and Wales: birth registrations" workbook (2025 edition). Four tables were used: live births by mother's area of usual residence and age (Table 2); maternities by country, place of birth, number of previous live-born children, and mother's age (Table 4); live births and age-specific fertility rates by country and parent age, 1938–2025 (Table 10); and mothers' standardised mean age by true birth order, 1938–2025 (Table 11). The **first-birth distribution** was operationalised as maternities with zero previous live-born children (Table 4), which is the closest published proxy for births to women who have not previously given birth; we note this counts maternities rather than live births and therefore treats a multiple first birth as one event.

**Forward projection.** Projected births by single year of mother's age for English local authorities, 2023–2047, were taken from the ONS 2022-based subnational population projections (births by age of mother, "persons" variant). English local authorities (area codes beginning "E") were summed to England totals for each projection year.

## Harmonisation of age scales

The cancer and childbearing distributions were placed on a common set of five-year age bands spanning under-15 to 90-and-over. Cancer's sub-15 bands were collapsed to a single "<15" category; ONS "Under 20" was treated as 15–19 and "50 and over" as 50–54 (births above 50 are negligible). Each distribution was normalised to sum to one over its full reported support before any comparison.

## Overlap coefficient

For each comparison we computed the overlap coefficient between the age density of cancer diagnoses, *f*c(a), and the age density of births, *f*b(a), defined on the shared harmonised bands as the sum over bands of the band-wise minimum, Σ min(*f*c, *f*b). This is the standard distributional overlap measure (range 0–1; 0 = disjoint, 1 = identical) and is the metric used in the analogous breast-cancer "collision" analysis. The headline comparison was the 2023 cancer distribution against the 2025 first-birth distribution (both England). The temporal trend, 2001–2023, used the all-births age distribution available for every year (England and Wales, Table 10) as the childbearing comparator, because a single-year-of-age first-birth distribution is not published England-only for the full back-series.

## Reproductive-age burden

As assumption-free complements to the overlap coefficient, we computed the annual number and proportion of female C54–C55 registrations occurring below ages 45 and 50, and the 90th percentile of the first-birth age distribution (by linear interpolation within bands), interpreted as the upper edge of contemporary first-time childbearing.

## Forward projection metrics

Holding the cancer diagnosis-age distribution fixed at 2023, the proportion-normalised overlap coefficient is mathematically pinned by cancer incidence within the reproductive zone — where births greatly outnumber cancers in proportional terms, the band-wise minimum equals the (fixed) cancer density — and is therefore nearly invariant to further childbearing postponement. We report this insensitivity explicitly. Accordingly, the projection is summarised through directly interpretable demographic quantities derived from the ONS projection: the mean age of mother, and the share and absolute number of births occurring at ages 35 and over and 40 and over, which measure the growth of the childbearing population entering the age range of material endometrial-cancer risk.

## Geography and validation

Cancer data are England-only throughout. Current-year childbearing distributions (first births and all births, 2025) are England-only (Tables 2 and 4, area of usual residence). The longitudinal childbearing comparator (Table 10) and the long mean-age series (Table 11) are England and Wales. To justify using England and Wales as a proxy for England in those series, we compared the standardised mean age at first birth for England against England and Wales in 2025, where both are published; the difference was 0.0 years, indicating that Wales (about 5% of births) does not materially shift the maternal-age distribution.

## Software and reproducibility

All analyses were performed in Python 3 using pandas, NumPy, Matplotlib, SciPy, and openpyxl. The complete pipeline runs from the raw downloaded files to every figure and reported statistic without manual intervention; the accompanying notebook reproduces the entire analysis and writes all outputs to a single directory.

## Limitations

The age-resolved cancer data are for all stages and all uterine histologies combined; the fertility-sparing-eligible subgroup (grade 1, stage IA, endometrioid) could be characterised only in aggregate (all-ages) form, so the overlap is computed on all uterine cancer and represents an upper envelope of the strictly eligible population. The first-birth measure is a maternities-based proxy. The longitudinal childbearing comparator is England and Wales rather than England-only (validated as above). Five-year banding limits resolution at the tails. Finally, the convergence quantified here is demographic, not aetiological: the analysis does not assert that delayed childbearing causes endometrial cancer (later age at last birth is in fact protective), but characterises the growing overlap between the populations affected by each, which is what determines clinical demand for uterus-conserving treatment.
