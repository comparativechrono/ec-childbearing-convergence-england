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

### Figures (all nine redrawn)
Apart from the corrected trend series and the added under-40 series (above), the plotted data
are unchanged. Colours now mean the same thing in every figure: green for uterine cancer, navy
for all births, orange for first births, and grey for a secondary series or an unknown
category.

- **Overlap figure (Fig 5 / figE).** The shading now covers the whole area under both curves.
  In v1 the shading joined the band-wise minima at the band midpoints with straight lines;
  between 42 and 47 years, where the two curves cross, that line cut the corner and left part
  of the shared area unshaded. The shading is now the pointwise minimum of the two plotted
  lines on a fine age grid. A second panel (b) enlarges ages 20–55, where the distributions meet.
  The coefficient (computed on the five-year bands) is unchanged.
  - Points now sit at the true band midpoints: ages are in completed years, so 15–19 is plotted
    at 17.5, not 17.
  - The births curve now starts at the 15–19 band, which holds ONS "Under 20". v1 plotted a zero
    for under-15 births; ONS publishes no separate group for those ages.
  - The legend gives the coefficient to three decimals (0.048).
- **Incidence trend (Fig 1 / figC).** The NDRS rolling-rate export has no rate for the
  2012–2014 window. v1 drew a straight line from 2011–2013 to 2013–2015, and the manuscript,
  README and methods called the series "gap-free". The line is now broken and the gap is
  labelled. The x-axis shows whole years (v1 showed 2002.5, 2005.0, …).
- **No text over data.** Every figure was checked for text (titles, labels, tick labels,
  legends, annotations) touching data lines, markers, shading, bars, axis lines, other text or
  the figure edge. Two overlaps in v1 are fixed:
  - the Fig 2 / figB trough label crossed the x-axis line; it now sits below the point, clear of
    the curve and the axis;
  - the Fig 6 / figF legend crossed the dotted under-50 series; y-axis headroom was added.

  Legends now sit in fixed positions with clear space, rather than being placed
  automatically. The Fig 3 / figA legend moves above the curves.
- **Second y-axes replaced by stacked panels (Figs 6 and 8 / figF and figH).** On a two-axis
  chart each series' scale is set independently. How closely the two lines tracked (Fig 6) or
  where they crossed (Fig 8) was therefore a consequence of the scaling, not a finding. Each
  quantity now has its own panel on a shared time axis.
- **Stage (Fig 4 / figD).** v1 omitted the 1,366 endometrial cancers with stage missing (1,338)
  or unstageable (28). They are now shown as an "Unknown" bar, and the two panels share one
  y-axis, so the stage bars add up to the endometrial bar.
- **Trough (Fig 2 / figB).** The published first-birth mean age (one decimal place) is 23.7 in
  1969, 1970 and 1971. The figure marks all three years, and `results_all.csv` reports
  "First-birth SMA trough years (E&W)" as 1969-1971; v1 reported only the first year.
- **Counts (Figs 7 and 9 / figG and figI).** The y-axes start at zero, so the series can be
  compared by height.
- **Labels.** Registry data are labelled "uterine cancer (C54–C55)". figD states the staging
  system (TNM 8th edition / FIGO 2009). Count axes have thousands separators (figD, figI).
  figA labels the age bands rather than plotting them at nominal midpoints.
- Figure font pinned to DejaVu Sans (bundled with matplotlib, and the font the v1 figures
  were drawn in) so figures render identically across environments.

- **Submission format.** Following BMC's figure guidelines, the figure files carry no titles
  (titles and legends belong in the manuscript), and two-panel figures are labelled a and b.
  Each figure is written as a closely cropped 300-dpi PNG and a vector PDF with embedded
  fonts; the PDFs carry no creation date, so reruns reproduce them byte for byte.

### Other
- Stale file names in the script header corrected.
- The notebook is regenerated from the script (cell-for-cell identical code).
