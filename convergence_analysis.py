# %% [markdown]
# # Endometrial cancer and delayed childbearing in England: a convergence analysis
#
# **Fully reproducible pipeline.** Starting from the raw public data files in `DATA_DIR`,
# this notebook reproduces every figure and number in the analysis:
# maternal-age distributions, the endometrial-cancer incidence trend, the
# endometrial/stage composition, the cancer–childbearing **age-overlap** (headline,
# trend, and reproductive-age burden), and a **forward projection to 2047** using the
# ONS 2022-based subnational birth projections.
#
# **Required raw files in `DATA_DIR`:**
# | file | source |
# |---|---|
# | `Incidence_National_ICD10_byage_singleyear.csv` | NDRS CancerData dashboard — incidence by ICD-10, sex, 5-year age band, single years, England, 2001–2023 |
# | `Incidence_National_ICD10_allages_rolling3yr.csv` | NDRS CancerData dashboard — all-ages 3-year rolling age-standardised rates, England, 2001–2023 |
# | `Cancer_registrations_2023_machine_readable_files.zip` | NHS Digital / NDRS — cancer registration statistics, England 2023 (Table 1 carries the endometrial vs non-endometrial × stage split) |
# | `2025birthregistrations.xlsx` | ONS — Births in England and Wales: birth registrations, 2025 edition |
# | `2022_SNPP_Births_persons.csv` | ONS — 2022-based subnational population projections, births by age of mother (persons) |
#
# Run top to bottom. Outputs are written to `OUT_DIR`.
#
# **Version 2 (October 2026).** Corrects the all-births trend comparator (ONS Table 10 reports
# its oldest maternal age group as "40 and over"; v1 expected "40 to 44"/"45 and over" and so
# silently dropped all births at 40+, which made the trend overlap equal to the share of
# cancers diagnosed under 40). Adds annual-percentage-change (APC) trend estimates, an
# under-40 series, a decomposition of the overlap coefficient, projected births at 40+, and a
# C54-only sensitivity analysis. All figures are redrawn: no text over data, overlap shading
# that covers the whole area under both curves, and the 2012–2014 gap in the rolling incidence
# series shown rather than bridged. See CHANGELOG.md.

# %%
# ============================ CONFIG & IMPORTS ============================
from pathlib import Path
import zipfile, io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Point DATA_DIR at the folder holding the four raw files (default: ./data next to this notebook).
DATA_DIR = Path("data")
OUT_DIR  = Path("outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

INCIDENCE_BYAGE = DATA_DIR / "Incidence_National_ICD10_byage_singleyear.csv"   # File 2: counts by 5-yr age band, single years (age/overlap work)
INCIDENCE_TREND = DATA_DIR / "Incidence_National_ICD10_allages_rolling3yr.csv" # File 1: NDRS 3-year rolling age-standardised rates (incidence-trend figure)
CANCER_MRF_ZIP = DATA_DIR / "Cancer_registrations_2023_machine_readable_files.zip"
BIRTHS_XLSX   = DATA_DIR / "2025birthregistrations.xlsx"
PROJ_BIRTHS   = DATA_DIR / "2022_SNPP_Births_persons.csv"

plt.rcParams.update({"font.family": "DejaVu Sans",   # bundled with matplotlib: identical figures everywhere
                     "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 140, "axes.titlesize": 12})
INK="#1a2e4a"; ACC="#b5651d"; ACC2="#2e7d6b"; GREY="#999"
results = {}   # collects every reported number

# %% [markdown]
# ## 1. Helper functions
# Harmonised 5-year age bands shared by the cancer and births distributions, the
# band-mapping dictionaries, and the **overlap coefficient** (the shared-support sum of
# the minimum of the two age densities).

# %%
HB = ['<15','15 to 19','20 to 24','25 to 29','30 to 34','35 to 39','40 to 44','45 to 49',
      '50 to 54','55 to 59','60 to 64','65 to 69','70 to 74','75 to 79','80 to 84','85 to 89','90 and over']
MIDHB = {'<15':12,'15 to 19':17,'20 to 24':22,'25 to 29':27,'30 to 34':32,'35 to 39':37,'40 to 44':42,
         '45 to 49':47,'50 to 54':52,'55 to 59':57,'60 to 64':62,'65 to 69':67,'70 to 74':72,
         '75 to 79':77,'80 to 84':82,'85 to 89':87,'90 and over':92}
REPRO_U45 = ['15 to 19','20 to 24','25 to 29','30 to 34','35 to 39','40 to 44']
REPRO_U50 = REPRO_U45 + ['45 to 49']

CANC_MAP = {'0 (Under 1)':'<15','01 to 04':'<15','05 to 09':'<15','10 to 14':'<15'}
for _b in HB[1:]: CANC_MAP[_b] = _b
BIRTH_MAP = {'Under 20':'15 to 19','20 to 24':'20 to 24','25 to 29':'25 to 29','30 to 34':'30 to 34',
             '35 to 39':'35 to 39','40 to 44':'40 to 44','45 to 49':'45 to 49','50 and over':'50 to 54'}
# ONS Table 10 (1938–2025 back-series) publishes mothers' ages as Under 20, 20–24, ..., 35–39 and
# "40 and over". Births at 40+ are placed in the 40–44 band (births at 45+ are <1% of births), so
# the trend overlap is, if anything, slightly conservative. ('40 to 44'/'45 and over' are kept in
# case a future edition splits the open-ended group.)
ASFR_MAP = {'Under 20':'15 to 19','20 to 24':'20 to 24','25 to 29':'25 to 29','30 to 34':'30 to 34',
            '35 to 39':'35 to 39','40 and over':'40 to 44','40 to 44':'40 to 44','45 and over':'45 to 49'}
REPRO_U40 = ['15 to 19','20 to 24','25 to 29','30 to 34','35 to 39']

def poisson_apc(counts, years, population=None):
    """Annual percentage change (APC) from a log-linear (quasi-)Poisson regression of annual counts
    on calendar year, with log(population) as offset when supplied (i.e. a trend in rates).
    95% CI uses the Pearson-based dispersion (floored at 1). Returns (APC, lower, upper)."""
    y = np.asarray(counts, float); x = np.asarray(years, float)
    X = np.column_stack([np.ones(len(y)), x - x.mean()])
    off = np.zeros(len(y)) if population is None else np.log(np.asarray(population, float))
    b = np.array([np.log(y.mean()) - off.mean(), 0.0])
    for _ in range(100):                                   # IRLS
        mu = np.exp(X @ b + off); z = X @ b + (y - mu) / mu
        XtW = X.T * mu; b_new = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-12: b = b_new; break
        b = b_new
    mu = np.exp(X @ b + off)
    phi = max(1.0, float(np.sum((y - mu) ** 2 / mu) / (len(y) - 2)))
    se = float(np.sqrt(np.linalg.inv((X.T * mu) @ X)[1, 1] * phi))
    pct = lambda v: 100 * (np.exp(v) - 1)
    return pct(b[1]), pct(b[1] - 1.96 * se), pct(b[1] + 1.96 * se)

def loglinear_apc(values, years):
    """APC (95% CI) from ordinary least squares of log(rate) on year (for published ASRs)."""
    from scipy import stats
    y = np.log(np.asarray(values, float)); x = np.asarray(years, float)
    X = np.column_stack([np.ones(len(y)), x - x.mean()])
    b = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ b
    se = float(np.sqrt((r @ r) / (len(y) - 2) * np.linalg.inv(X.T @ X)[1, 1]))
    t = stats.t.ppf(0.975, len(y) - 2); pct = lambda v: 100 * (np.exp(v) - 1)
    return pct(b[1]), pct(b[1] - t * se), pct(b[1] + t * se)

def to_hb(series, mapping):
    """Aggregate a band-indexed COUNT Series onto harmonised bands; return proportions over HB."""
    out = pd.Series(0.0, index=HB)
    for lab, val in series.items():
        if lab in mapping and pd.notna(val):
            out[mapping[lab]] += float(val)
    return out / out.sum()

def OVL(pc, pb):
    """Overlap coefficient: sum over shared harmonised bands of min(p_cancer, p_birth)."""
    pc = pc.reindex(HB).fillna(0); pb = pb.reindex(HB).fillna(0)
    return float(np.minimum(pc, pb).sum())

def banded_mean(p):
    return float(sum(p[b]*MIDHB[b] for b in p.index if b in MIDHB))

def find_header_row(xlsx, sheet, sentinel, max_scan=12):
    """Locate the header row (0-indexed) by finding a row that contains the sentinel AND
    has several populated cells (ONS title rows contain only one populated cell)."""
    probe = pd.read_excel(xlsx, sheet_name=sheet, header=None, nrows=max_scan, engine="openpyxl")
    for i in range(len(probe)):
        row = probe.iloc[i]
        has_sentinel = row.astype(str).str.contains(sentinel, case=False, na=False, regex=False).any()
        if has_sentinel and row.notna().sum() >= 3:
            return i
    raise ValueError(f"Header sentinel '{sentinel}' not found in sheet {sheet}")

def read_ons_table(xlsx, sheet, sentinel):
    hdr = find_header_row(xlsx, sheet, sentinel)
    df = pd.read_excel(xlsx, sheet_name=sheet, header=hdr, engine="openpyxl")
    df.columns = [str(c).strip() for c in df.columns]
    return df

# %% [markdown]
# ## 2. Cancer age distribution (NDRS dashboard export)
# Endometrial/uterine cancer (ICD-10 **C54–C55**, females) registrations by 5-year age
# band, England, single years 2001–2023.

# %%
inc = pd.read_csv(INCIDENCE_BYAGE)
inc.columns = [c.strip() for c in inc.columns]
inc["Count"] = pd.to_numeric(inc["Count"], errors="coerce")

def cancer_counts(year, code="C54-C55"):
    """Return a band-indexed COUNT Series of female cancer registrations for one year."""
    d = inc[(inc["ICD10 code"]==code) & (inc["Gender"]=="Females") & (inc["Year"]==year) &
            (inc["Age at diagnosis"]!="All ages")].drop_duplicates(subset=["Age at diagnosis"])
    return d.set_index("Age at diagnosis")["Count"]

def cancer_dist(year, code="C54-C55"):
    return to_hb(cancer_counts(year, code), CANC_MAP)

YEARS = list(range(2001, 2024))
pc_2023 = cancer_dist(2023)
raw_2023 = cancer_counts(2023)
assert abs(raw_2023.sum() - 9048) < 1, "2023 C54-C55 female total should reconcile to all-ages count"

# Incidence trend uses NDRS 3-year rolling age-standardised rates (File 1), which smooth
# year-to-year variation. The single-year file has no 2012 rate (count present, rate blank), and
# the rolling file has no rate for the 2012-2014 window; Fig C shows that window as a gap.
roll = pd.read_csv(INCIDENCE_TREND); roll.columns = [c.strip() for c in roll.columns]
asr = roll[(roll["ICD10 code"]=="C54-C55") & (roll["Gender"]=="Females") &
           (roll["Type of rate"]=="Age-standardised")].copy()
asr["Rate"] = pd.to_numeric(asr["Rate"], errors="coerce")
asr["mid"]  = asr["Year"].astype(str).str.split("-").str[0].astype(int) + 1   # window midpoint
asr = asr.dropna(subset=["Rate"]).sort_values("mid")
results["EC ASR 2001-2003 (3-yr rolling)"] = round(float(asr[asr["mid"]==2002]["Rate"].iloc[0]), 1)
results["EC ASR 2021-2023 (3-yr rolling)"] = round(float(asr[asr["mid"]==2022]["Rate"].iloc[0]), 1)
results["EC ASR rolling change 2001-03 -> 2021-23 (%)"] = round(
    100*(asr[asr["mid"]==2022]["Rate"].iloc[0]/asr[asr["mid"]==2002]["Rate"].iloc[0]-1), 1)

# Single-year all-ages age-standardised rates (2012 is blank in the export) — used for the
# all-ages APC, and the 2018 value for comparison with published UK figures.
inc["Rate_num"] = pd.to_numeric(inc["Rate"], errors="coerce")
asr1 = (inc[(inc["ICD10 code"]=="C54-C55") & (inc["Gender"]=="Females") &
            (inc["Age at diagnosis"]=="All ages") & (inc["Type of rate"]=="Age-standardised")]
        .drop_duplicates("Year").set_index("Year")["Rate_num"].dropna().sort_index())
results["EC ASR 2018 (single year)"] = round(float(asr1.loc[2018]), 1)

# Female population by age band and year, recovered from the NDRS all-cancer age-specific
# rates (population = count x 100,000 / rate; NDRS rates use ONS mid-year estimates). The
# all-cancer counts are large, so this is precise to <0.05 per 100,000 when checked against
# the published uterine-cancer rates. 2012 rates are blank in the export, so the 2012
# population is interpolated geometrically between 2011 and 2013.
_allc = inc[(inc["ICD10 code"]=="C00-C97 excl. C44") & (inc["Gender"]=="Females") &
            (inc["Type of rate"]=="Age-specific")].copy()
_allc["pop"] = _allc["Count"] * 1e5 / _allc["Rate_num"]
POP = _allc.pivot_table(index="Year", columns="Age at diagnosis", values="pop", aggfunc="first")
POP = POP.replace([np.inf, -np.inf], np.nan)
POP.loc[2012] = np.sqrt(POP.loc[2011] * POP.loc[2013])
POP = POP.sort_index()

# %% [markdown]
# ## 3. Births: distributions and trends (ONS birth registrations workbook)
# First-birth proxy = maternities with **zero previous live-born children** (Table 4).
# All-births age distribution from Table 2; long age-specific series from Table 10;
# standardised mean age by true birth order from Table 11.

# %%
t2  = read_ons_table(BIRTHS_XLSX, "Table_2",  "Age of Mother at birth")
t4  = read_ons_table(BIRTHS_XLSX, "Table_4",  "Age of Mother at birth")
t10 = read_ons_table(BIRTHS_XLSX, "Table_10", "Age group (years)")
t11 = read_ons_table(BIRTHS_XLSX, "Table_11", "True birth order first")

# First births (England 2025): maternities, 0 previous children, any place of birth
agecol_t4 = [c for c in t4.columns if c.startswith("Age of Mother")][0]
prevcol   = [c for c in t4.columns if "previous live-born" in c][0]
matcol    = [c for c in t4.columns if "Maternities" in c][0]
fb = t4[(t4["Country"].astype(str).str.strip()=="England") &
        (t4["Place of birth"].astype(str).str.strip()=="Any place of birth") &
        (t4[prevcol].astype(str).str.strip()=="0")]
firstb = fb.set_index(agecol_t4)[matcol].astype(float)
p_first = to_hb(firstb, BIRTH_MAP)

# All births (England 2025): Table 2
agecol_t2 = [c for c in t2.columns if c.startswith("Age of Mother")][0]
namecol_t2 = [c for c in t2.columns if "Name" in c][0]
allb = (t2[t2[namecol_t2].astype(str).str.strip()=="ENGLAND"]
        .set_index(agecol_t2)["Number of live births"].astype(float))
p_all = to_hb(allb, BIRTH_MAP)

# Long all-births age distribution per year (E&W) from Table 10 (Mother rows)
t10m = t10[(t10["Parent"]=="Mother") & (t10["Country"].astype(str).str.contains("England, Wales", na=False))].copy()
t10m["Year"] = pd.to_numeric(t10m["Year"], errors="coerce")
t10m["nb"]   = pd.to_numeric(t10m["Number of live births"], errors="coerce")
def allbirths_ew(year):
    d = t10m[(t10m["Year"]==year) & (t10m["Age group (years)"].isin(ASFR_MAP))]
    return to_hb(d.set_index("Age group (years)")["nb"], ASFR_MAP)

# SMA first-birth trend (E&W) from Table 11
namecol_t11 = [c for c in t11.columns if "Name" in c][0]
fbcol  = [c for c in t11.columns if "first" in c][0]
allcol = [c for c in t11.columns if "All births" in c][0]
t11ew = t11[t11[namecol_t11].astype(str).str.contains("England, Wales", na=False)].copy()
t11ew["Year"] = pd.to_numeric(t11ew["Year"], errors="coerce")
t11ew = t11ew.dropna(subset=["Year"]).sort_values("Year")
sma_first = pd.to_numeric(t11ew[fbcol], errors="coerce")
sma_all   = pd.to_numeric(t11ew[allcol], errors="coerce")

# England-vs-E&W validation of the first-birth SMA (2025)
def _val(df, where, col):
    r = df[df[namecol_t11].astype(str).str.strip()==where]
    return float(pd.to_numeric(r[col], errors="coerce").iloc[0]) if len(r) else np.nan
t11_2025 = t11[pd.to_numeric(t11["Year"], errors="coerce")==2025]
gap = _val(t11_2025, "England", fbcol) - _val(t11_2025, "England, Wales and Elsewhere", fbcol)

results["England 2025 all-births total"]            = int(allb[allb.index.isin(BIRTH_MAP)].sum())
results["England 2025 all-births banded mean age"]  = round(banded_mean(p_all), 2)
results["England 2025 first-births total (proxy)"]  = int(firstb[firstb.index.isin(BIRTH_MAP)].sum())
results["England 2025 first-births banded mean age"]= round(banded_mean(p_first), 2)
results["England 2025 first-births % aged 35+"]     = round(100*p_first[['35 to 39','40 to 44','45 to 49']].sum(), 1)
results["England 2025 first-births % aged 30+"]     = round(100*p_first[['30 to 34','35 to 39','40 to 44','45 to 49']].sum(), 1)
results["First-birth SMA 1975 (E&W)"]               = round(float(sma_first[t11ew["Year"]==1975].iloc[0]), 2)
results["First-birth SMA 2025 (E&W)"]               = round(float(sma_first[t11ew["Year"]==2025].iloc[0]), 2)
results["All-births SMA 2025 (E&W)"]                = round(float(sma_all[t11ew["Year"]==2025].iloc[0]), 2)
results["First-birth SMA trough year (E&W)"]        = int(t11ew.loc[sma_first.idxmin(), "Year"])
results["First-birth SMA trough value (E&W)"]       = round(float(sma_first.min()), 2)
results["First-birth SMA England vs E&W gap, 2025 (yrs)"] = round(gap, 3)

# %% [markdown]
# ## 4. Endometrial vs non-endometrial composition and stage (NHS machine-readable Table 1, 2023)

# %%
with zipfile.ZipFile(CANCER_MRF_ZIP) as z:
    t1name = [n for n in z.namelist() if "Table_1" in n][0]
    with z.open(t1name) as fh:
        t1 = pd.read_csv(fh)
t1["count"] = pd.to_numeric(t1["count"], errors="coerce")
ut = t1[(t1["ndrs_detailed_group"].isin(["All Uterus","Uterus - endometrial","Uterus - non-endometrial"])) &
        (t1["gender"]=="Females")].drop_duplicates(subset=["ndrs_detailed_group","stage_at_diagnosis"])
def uval(group, stage): 
    r = ut[(ut["ndrs_detailed_group"]==group) & (ut["stage_at_diagnosis"]==stage)]["count"]
    return float(r.iloc[0]) if len(r) else 0.0
all_u   = uval("All Uterus","All stages")
endo    = uval("Uterus - endometrial","All stages")
nonendo = uval("Uterus - non-endometrial","All stages")
endo_stage = {s: uval("Uterus - endometrial", s) for s in ["1","2","3","4"]}
endo_staged = sum(endo_stage.values())
results["2023 all uterine (registrations)"]       = int(all_u)
results["2023 endometrial (registrations)"]       = int(endo)
results["2023 % endometrial of all uterine"]      = round(100*endo/all_u, 1)
results["2023 endometrial stage 1 count"]         = int(endo_stage["1"])
results["2023 % endometrial stage 1 (of staged)"] = round(100*endo_stage["1"]/endo_staged, 1)
results["2023 % endometrial stage 1 (of all endo)"]= round(100*endo_stage["1"]/endo, 1)
# Stage is as registered by NDRS: TNM 8th edition (equivalent to FIGO 2009 for the corpus
# uteri), the version used for diagnoses up to 31/12/2025.
results["2023 non-endometrial (registrations)"]   = int(nonendo)
results["2023 % non-endometrial of all uterine"]  = round(100*nonendo/all_u, 1)
results["2023 endometrial stage missing/unstageable"] = int(uval("Uterus - endometrial","Missing") +
                                                          uval("Uterus - endometrial","Unstageable"))

# %% [markdown]
# ## 5. Overlap: headline, reproductive-age burden, and trend
# Headline overlap = cancer diagnosis-age density (2023) vs first-birth age density (2025).
# Trend uses the all-births comparator available longitudinally (E&W, Table 10).

# %%
results["OVL: EC 2023 vs first-births (England 2025)"] = round(OVL(pc_2023, p_first), 3)
results["OVL: EC 2023 vs all-births (England 2025)"]   = round(OVL(pc_2023, p_all), 3)
pc_2123 = sum(cancer_dist(y) for y in [2021,2022,2023]) / 3
results["OVL: EC 2021-2023 vs first-births (2025)"]    = round(OVL(pc_2123, p_first), 3)

tot23 = raw_2023.sum()
cnt = lambda s, bands: float(sum(s.get(b, 0) for b in bands))
results["EC 2023 count age<45"] = int(cnt(raw_2023, REPRO_U45))
results["EC 2023 count age<50"] = int(cnt(raw_2023, REPRO_U50))
results["EC 2023 % age<45"]     = round(100*cnt(raw_2023, REPRO_U45)/tot23, 1)
results["EC 2023 % age<50"]     = round(100*cnt(raw_2023, REPRO_U50)/tot23, 1)

# Decomposition of the headline overlap. In each band the overlap takes the smaller of the two
# percentages. Below 45 the cancer percentage is the smaller; from 45 upwards the birth
# percentage is. So the coefficient ~= (share of cancers diagnosed 15-44) + (share of births 45+).
_cl = [b for b in HB if pc_2023[b] <= p_first[b]]          # bands limited by cancer
_bl = [b for b in HB if p_first[b] < pc_2023[b]]           # bands limited by births
results["OVL decomposition: cancer-limited bands"] = ", ".join(_cl)
results["OVL decomposition: cancer share in cancer-limited bands"] = round(float(pc_2023[_cl].sum()), 4)
results["OVL decomposition: birth share in birth-limited bands"]   = round(float(p_first[_bl].sum()), 4)
results["Cancer % aged 40-44 (2023)"]       = round(100*float(pc_2023['40 to 44']), 1)
results["First-birth % aged 40-44 (2025)"]  = round(100*float(p_first['40 to 44']), 1)
results["Cancer % aged 45-49 (2023)"]       = round(100*float(pc_2023['45 to 49']), 1)
results["First-birth % aged 45+ (2025)"]    = round(100*float(p_first[['45 to 49','50 to 54']].sum()), 1)

# 90th percentile of first-birth age (banded linear interpolation)
def band_pctile(p, q):
    edges = {'15 to 19':(15,20),'20 to 24':(20,25),'25 to 29':(25,30),'30 to 34':(30,35),
             '35 to 39':(35,40),'40 to 44':(40,45),'45 to 49':(45,50),'50 to 54':(50,55)}
    cum = 0.0
    for b in HB:
        if b not in edges:
            cum += p.get(b, 0); continue
        pb = p.get(b, 0)
        if cum + pb >= q:
            lo, hi = edges[b]; return lo + (q-cum)/pb*(hi-lo)
        cum += pb
    return 50.0
results["First-birth 90th percentile age (England 2025)"] = round(band_pctile(p_first, 0.90), 1)

# Trend 2001-2023
ovl_all=[]; ovl_firstfixed=[]; pct_u50=[]; pct_u45=[]; pct_u40=[]; cnt_u40=[]; cnt_u45=[]; cnt_u50=[]; tot_all=[]
for y in YEARS:
    pcy = cancer_dist(y); rawy = cancer_counts(y); ty = rawy.sum()
    ovl_all.append(OVL(pcy, allbirths_ew(y)))
    ovl_firstfixed.append(OVL(pcy, p_first))
    pct_u50.append(100*cnt(rawy, REPRO_U50)/ty); pct_u45.append(100*cnt(rawy, REPRO_U45)/ty)
    pct_u40.append(100*cnt(rawy, REPRO_U40)/ty)
    cnt_u40.append(cnt(rawy, REPRO_U40)); cnt_u45.append(cnt(rawy, REPRO_U45)); cnt_u50.append(cnt(rawy, REPRO_U50))
    tot_all.append(float(ty))
results["OVL(all-births) 2001"] = round(ovl_all[0], 3)
results["OVL(all-births) 2016"] = round(ovl_all[YEARS.index(2016)], 3)
results["OVL(all-births) 2023"] = round(ovl_all[-1], 3)
results["OVL(all-births) change 2001->2023 (%)"] = round(100*(ovl_all[-1]/ovl_all[0]-1), 1)
results["OVL(first-births 2025 fixed) 2001"] = round(ovl_firstfixed[0], 3)
results["OVL(first-births 2025 fixed) 2023"] = round(ovl_firstfixed[-1], 3)
results["EC registrations all ages: 2001"] = int(tot_all[0]); results["EC registrations all ages: 2023"] = int(tot_all[-1])
results["EC count age<40: 2001"] = int(cnt_u40[0]); results["EC count age<40: 2023"] = int(cnt_u40[-1])
results["EC count age<40: 2018"] = int(cnt_u40[YEARS.index(2018)])   # for comparison with Ind (2023): ~150 UK cases <40
results["EC count age<45: 2001"] = int(cnt_u45[0]); results["EC count age<45: 2023"] = int(cnt_u45[-1])
results["EC count age<50: 2001"] = int(cnt_u50[0]); results["EC count age<50: 2023"] = int(cnt_u50[-1])
results["EC % age<40: 2001"] = round(pct_u40[0], 1); results["EC % age<40: 2023"] = round(pct_u40[-1], 1)
results["EC % age<45: 2001"] = round(pct_u45[0], 1); results["EC % age<45: 2023"] = round(pct_u45[-1], 1)
results["EC % age<50: 2001"] = round(pct_u50[0], 1); results["EC % age<50: 2023"] = round(pct_u50[-1], 1)
_t10_40 = lambda y: float(t10m[(t10m["Year"]==y) & (t10m["Age group (years)"]=="40 and over")]["nb"].iloc[0])
_t10_tot = lambda y: float(t10m[(t10m["Year"]==y) & (t10m["Age group (years)"].isin(ASFR_MAP))]["nb"].sum())
results["E&W % of births to mothers 40+: 2001"] = round(100*_t10_40(2001)/_t10_tot(2001), 1)
results["E&W % of births to mothers 40+: 2023"] = round(100*_t10_40(2023)/_t10_tot(2023), 1)

# %% [markdown]
# ### 5b. Pace of change: annual percentage change (APC), 2001–2023
# Log-linear (quasi-)Poisson regression of annual registrations on calendar year, without an
# offset (counts) and with log(female population) as offset (age-specific rates for the
# combined bands). The all-ages comparator is the APC of the published single-year
# age-standardised rate (log-linear least squares; 2012 is blank in the export).

# %%
def _band_series(bands):
    c = pd.Series([cnt(cancer_counts(y), bands) for y in YEARS], index=YEARS)
    p = POP.loc[YEARS, bands].sum(axis=1)
    return c, p
apc_rows = []
for label, bands, ages in [("age<40", REPRO_U40, "15-39"), ("age<45", REPRO_U45, "15-44"), ("age<50", REPRO_U50, "15-49")]:
    c, p = _band_series(bands)
    a_cnt = poisson_apc(c.values, YEARS); a_rate = poisson_apc(c.values, YEARS, p.values)
    r01, r23 = 1e5*c[2001]/p[2001], 1e5*c[2023]/p[2023]
    results[f"EC rate {label} per 100,000 women aged {ages}: 2001"] = round(r01, 2)
    results[f"EC rate {label} per 100,000 women aged {ages}: 2023"] = round(r23, 2)
    results[f"APC counts {label} (%/yr, 95% CI)"] = f"{a_cnt[0]:.1f} ({a_cnt[1]:.1f} to {a_cnt[2]:.1f})"
    results[f"APC rate {label} (%/yr, 95% CI)"]   = f"{a_rate[0]:.1f} ({a_rate[1]:.1f} to {a_rate[2]:.1f})"
    apc_rows.append((label, a_cnt, a_rate))
a_all = poisson_apc(tot_all, YEARS)
a_asr = loglinear_apc(asr1.values, asr1.index.values)
results["APC counts all ages (%/yr, 95% CI)"] = f"{a_all[0]:.1f} ({a_all[1]:.1f} to {a_all[2]:.1f})"
results["APC age-standardised rate all ages (%/yr, 95% CI)"] = f"{a_asr[0]:.1f} ({a_asr[1]:.1f} to {a_asr[2]:.1f})"

# %% [markdown]
# ### 5c. Sensitivity: corpus uteri only (ICD-10 C54)
# C55 ("uterus, part unspecified") is excluded to check that the results are not driven by
# tumours of unspecified uterine origin.

# %%
c54_all = inc[(inc["ICD10 code"]=="C54") & (inc["Gender"]=="Females") &
              (inc["Age at diagnosis"]=="All ages")].drop_duplicates("Year").set_index("Year")["Count"]
_c5455_all = pd.Series(tot_all, index=YEARS)
results["C55 % of C54-C55 registrations: 2001"] = round(100*(1 - c54_all.loc[2001]/_c5455_all.loc[2001]), 1)
results["C55 % of C54-C55 registrations: 2023"] = round(100*(1 - c54_all.loc[2023]/_c5455_all.loc[2023]), 1)
results["C54 only - OVL vs first-births (2023 vs 2025)"] = round(OVL(cancer_dist(2023, "C54"), p_first), 3)
_ovl54 = [OVL(cancer_dist(y, "C54"), allbirths_ew(y)) for y in YEARS]
results["C54 only - OVL(all-births) 2001"] = round(_ovl54[0], 3)
results["C54 only - OVL(all-births) 2023"] = round(_ovl54[-1], 3)
for label, bands in [("age<40", REPRO_U40), ("age<45", REPRO_U45), ("age<50", REPRO_U50)]:
    results[f"C54 only - count {label}: 2001"] = int(cnt(cancer_counts(2001, "C54"), bands))
    results[f"C54 only - count {label}: 2023"] = int(cnt(cancer_counts(2023, "C54"), bands))

# %% [markdown]
# ## 6. Forward projection to 2047 (ONS 2022-based subnational birth projections)
# England births by single year of mother's age are summed across English local
# authorities for each projection year (the "persons" file = all babies; the maternal-age
# distribution is identical to the sex-specific files).
#
# **Interpretation note.** The proportion-normalised overlap coefficient is *pinned by
# endometrial-cancer incidence* in the reproductive zone (where births vastly outnumber
# cancers, the band-wise minimum equals the cancer density), so it is nearly flat under a
# birth-side shift alone. The overlap's future therefore depends chiefly on whether young
# EC incidence keeps rising — not on childbearing postponement per se. We therefore report
# the projection through the sensitive, directly interpretable demographic quantities: the
# mean maternal age and the share/number of births occurring at ages 35+ and 40+, i.e. the
# growth of the childbearing population entering the EC-risk window.
#
# **Projection age range.** The ONS projections allocate all births to mothers aged 15-44 (the
# single-age figures sum exactly to the all-ages total), so "40+" in the projection means the
# 40-44 band and no projected births fall at 45+. With the cancer distribution held at 2023,
# the projected overlap therefore equals the share of cancers diagnosed at 15-44 in every year.

# %%
proj = pd.read_csv(PROJ_BIRTHS)
proj = proj[proj["AREA_CODE"].astype(str).str.startswith("E")]          # England LAs
proj = proj[proj["AGE_GROUP"].astype(str).str.fullmatch(r"\d+")].copy() # numeric ages only
proj["age"] = proj["AGE_GROUP"].astype(int)
PROJ_YEARS = [c for c in proj.columns if str(c).isdigit() and len(str(c))==4]
proj_years_int = [int(y) for y in PROJ_YEARS]

AGE_TO_HB = {}
for a in range(10, 100):
    for b in HB:
        if b=='<15' and a<15: AGE_TO_HB[a]=b
        elif b=='90 and over' and a>=90: AGE_TO_HB[a]=b
        elif b not in ('<15','90 and over'):
            lo=int(b.split(' to ')[0])
            if lo<=a<=lo+4: AGE_TO_HB[a]=b

def proj_year_age(year):
    return proj.groupby("age")[year].sum()

def proj_birth_dist(year):
    g = proj_year_age(year); s = pd.Series(0.0, index=HB)
    for a, v in g.items():
        if a in AGE_TO_HB: s[AGE_TO_HB[a]] += float(v)
    return s/s.sum()

def proj_mean_age(year):
    g = proj_year_age(year); return float((g.index*g).sum()/g.sum())

# Demographic projection metrics (sensitive, directly interpretable)
proj_mage   = [proj_mean_age(y) for y in PROJ_YEARS]
proj_share35= [proj_year_age(y)[proj_year_age(y).index>=35].sum()/proj_year_age(y).sum() for y in PROJ_YEARS]
proj_share40= [proj_year_age(y)[proj_year_age(y).index>=40].sum()/proj_year_age(y).sum() for y in PROJ_YEARS]
proj_cnt35  = [proj_year_age(y)[proj_year_age(y).index>=35].sum() for y in PROJ_YEARS]
proj_cnt40  = [proj_year_age(y)[proj_year_age(y).index>=40].sum() for y in PROJ_YEARS]
# Overlap with cancer held fixed at 2023 (reported for completeness / rigor)
proj_ovl    = [OVL(pc_2023, proj_birth_dist(y)) for y in PROJ_YEARS]

def at(year, arr): return arr[proj_years_int.index(year)]
results["Projected mean maternal age 2023"]      = round(at(2023, proj_mage), 2)
results["Projected mean maternal age 2040"]      = round(at(2040, proj_mage), 2)
results["Projected mean maternal age 2047"]      = round(at(2047, proj_mage), 2)
results["Projected % births to age 35+ : 2023"]  = round(100*at(2023, proj_share35), 1)
results["Projected % births to age 35+ : 2040"]  = round(100*at(2040, proj_share35), 1)
results["Projected % births to age 35+ : 2047"]  = round(100*at(2047, proj_share35), 1)
results["Projected % births to age 40+ : 2023"]  = round(100*at(2023, proj_share40), 1)
results["Projected % births to age 40+ : 2047"]  = round(100*at(2047, proj_share40), 1)
results["Projected births to age 35+ : 2023"]    = int(at(2023, proj_cnt35))
results["Projected births to age 35+ : 2047"]    = int(at(2047, proj_cnt35))
results["Projected births to age 40+ : 2023"]    = int(at(2023, proj_cnt40))
results["Projected births to age 40+ : 2047"]    = int(at(2047, proj_cnt40))
results["Projected births, all ages : 2023"]     = int(proj_year_age('2023').sum())
results["Projected births, all ages : 2047"]     = int(proj_year_age('2047').sum())
_all_check = pd.read_csv(PROJ_BIRTHS)
_all_check = _all_check[_all_check["AREA_CODE"].astype(str).str.startswith("E") & (_all_check["AGE_GROUP"]=="All ages")]
assert abs(_all_check['2023'].sum() - proj_year_age('2023').sum()) < 1, "single-age births should sum to the all-ages total"
# Honest rigor note: proportion-normalised overlap is pinned by EC incidence in the
# reproductive zone, hence nearly flat under birth-side shift alone.
results["Projected overlap change 2023->2047 (cancer fixed, %)"] = round(100*(at(2047,proj_ovl)/at(2023,proj_ovl)-1), 1)

# %% [markdown]
# ## 7. Figures

# %%
# Fig A — maternal age distribution, England 2025
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle
fig, ax = plt.subplots(figsize=(7.2,4.4))
xr = [MIDHB[b] for b in REPRO_U50]
ax.plot(xr, [p_all.get(b,0) for b in REPRO_U50], "-o", color=INK, lw=2, label="All births")
ax.plot(xr, [p_first.get(b,0) for b in REPRO_U50], "-s", color=ACC, lw=2, label="First births (0 previous children)")
ax.set_xlabel("Mother's age (band midpoint, years)"); ax.set_ylabel("Proportion of births")
ax.set_ylim(bottom=-0.015, top=1.32*max(max(p_all.get(b,0) for b in REPRO_U50), max(p_first.get(b,0) for b in REPRO_U50)))  # legend headroom
ax.set_title("Maternal age distribution, England, 2025"); ax.legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(OUT_DIR/"figA_maternal_age_distribution_2025.png"); plt.close(fig)

# Fig B — first-birth SMA trend (E&W); the trough label sits below the point, clear of the curve and axis
fig, ax = plt.subplots(figsize=(7.2,4.4))
ax.plot(t11ew["Year"], sma_first, color=ACC, lw=2, label="Mean age at first birth")
ax.plot(t11ew["Year"], sma_all, color=GREY, lw=1.5, ls="--", label="Mean age, all births")
ty=int(t11ew.loc[sma_first.idxmin(),"Year"]); tv=float(sma_first.min())
ax.scatter([ty],[tv], color=INK, zorder=5)
ax.set_ylim(bottom=tv - 1.1)
ax.annotate(f"trough {ty}: {tv:.1f} years", (ty, tv), textcoords="offset points", xytext=(0, -9),
            ha="center", va="top", fontsize=9, color=INK)
ax.set_xlabel("Year"); ax.set_ylabel("Standardised mean age of mother (years)")
ax.set_title("Mean age at first birth, England & Wales, 1938–2025"); ax.legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(OUT_DIR/"figB_first_birth_age_trend.png"); plt.close(fig)

# Fig C — uterine incidence ASR trend (NDRS 3-year rolling rates). The 2012–2014 window has no
# rate in the export, so the line is broken there rather than bridged.
fig, ax = plt.subplots(figsize=(7.2,4.4))
_mids = list(range(int(asr["mid"].min()), int(asr["mid"].max()) + 1))
_roll = asr.set_index("mid")["Rate"].reindex(_mids)
ax.plot(_mids, _roll.values, "-o", color=ACC2, lw=2)
for _m in _roll.index[_roll.isna()]:
    _y = float(np.nanmean([_roll.get(_m-1, np.nan), _roll.get(_m+1, np.nan)]))
    ax.annotate(f"{_m-1}–{str(_m+1)[2:]}\nnot available", (_m, _y), textcoords="offset points", xytext=(0, -22),
                ha="center", va="top", fontsize=8, color="#666")
ax.set_xticks(range(_mids[0], _mids[-1] + 1, 4))
ax.set_xlabel("Year (midpoint of 3-year rolling window)"); ax.set_ylabel("Age-standardised rate / 100,000 women")
ax.set_title("Uterine cancer (C54–C55) incidence, England, 2001–2023\n(3-year rolling rates)")
fig.tight_layout(); fig.savefig(OUT_DIR/"figC_uterine_incidence_trend.png"); plt.close(fig)

# Fig D — endometrial composition + stage
fig, (a1,a2) = plt.subplots(1,2, figsize=(9,4.2))
a1.bar(["Endometrial","Non-endometrial"], [endo,nonendo], color=[INK,GREY])
a1.set_title(f"Uterine cancer (C54–C55) by NDRS group, 2023\n(endometrial = {100*endo/all_u:.0f}% of {int(all_u):,})")
a1.set_ylabel("Registrations"); a1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
a2.bar(["1","2","3","4"], [endo_stage[s] for s in ["1","2","3","4"]], color=ACC)
a2.set_title("Endometrial cancer by stage, 2023"); a2.set_xlabel("Stage at diagnosis (TNM 8 / FIGO 2009)"); a2.set_ylabel("Registrations")
a2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
fig.tight_layout(); fig.savefig(OUT_DIR/"figD_endometrial_type_stage_2023.png"); plt.close(fig)

# Fig E — headline overlap. The shaded area is the region under BOTH plotted curves, i.e. the
# minimum of the two piecewise-linear curves evaluated on a fine grid (including the points where
# they cross between band midpoints). The coefficient itself is computed on the five-year bands.
# Lower panel: ages 20–55 enlarged, where the two distributions meet.
fig, (ax, axz) = plt.subplots(2, 1, figsize=(8, 6.6), gridspec_kw={"height_ratios": [1.6, 1]})
_x = np.array([MIDHB[b] for b in HB], float)
_yc = np.array([pc_2023.get(b,0) for b in HB], float); _yb = np.array([p_first.get(b,0) for b in HB], float)
_xd = np.linspace(_x.min(), _x.max(), 4001)
_shared = np.minimum(np.interp(_xd, _x, _yc), np.interp(_xd, _x, _yb))
for a in (ax, axz):
    a.fill_between(_xd, _shared, color=INK, alpha=0.35, lw=0, label=f"Overlap (coefficient = {OVL(pc_2023,p_first):.2f})")
    a.plot(_x, _yc, "-o", color=ACC2, lw=2, label="Uterine cancer (C54–C55), age at diagnosis (England 2023)")
    a.plot(_x, _yb, "-s", color=ACC, lw=2, label="First births, mother's age (England 2025)")
ax.set_xlim(10, 95); ax.set_ylabel("Proportion within distribution")
ax.add_patch(Rectangle((20, 0), 35, 0.05, fill=False, ls=":", lw=1, ec="#666"))
ax.set_title("Convergence of childbearing age and uterine cancer diagnosis age, England")
_h, _l = ax.get_legend_handles_labels()
ax.legend([_h[1], _h[2], _h[0]], [_l[1], _l[2], _l[0]], frameon=False, fontsize=8.5, loc="upper right")
axz.set_xlim(20, 55); axz.set_ylim(0, 0.05)
axz.set_xlabel("Age (years)"); axz.set_ylabel("Proportion")
axz.set_title("Ages 20–55 enlarged (dotted box above)", fontsize=10, pad=12)
fig.tight_layout(); fig.savefig(OUT_DIR/"figE_overlap_headline.png"); plt.close(fig)

# Fig F — overlap trend + % under 50 (right axis), all three series in one legend
fig, ax = plt.subplots(figsize=(7.6,4.4))
_l1, = ax.plot(YEARS, ovl_all, "-o", color=INK, lw=2, label="Overlap with all births (E&W, each year)")
_l2, = ax.plot(YEARS, ovl_firstfixed, "-^", color=GREY, lw=1.6, label="Overlap with first births (England 2025, held fixed)")
ax.set_xlabel("Diagnosis year"); ax.set_ylabel("Overlap coefficient")
ax2=ax.twinx(); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_color(ACC)
_l3, = ax2.plot(YEARS, pct_u50, color=ACC, lw=1.6, ls=":", label="% of uterine cancers diagnosed under 50 (right axis)")
ax2.set_ylabel("% of uterine cancers diagnosed under 50", color=ACC); ax2.tick_params(axis="y", colors=ACC)
_lo = min(min(ovl_all), min(ovl_firstfixed)); _hi = max(max(ovl_all), max(ovl_firstfixed))
ax.set_ylim(_lo - 0.05*(_hi-_lo), _hi + 0.55*(_hi-_lo))            # headroom for the legend
ax2.set_ylim(min(pct_u50) - 0.05*(max(pct_u50)-min(pct_u50)), max(pct_u50) + 0.55*(max(pct_u50)-min(pct_u50)))
ax.legend([_l1, _l2, _l3], [_l1.get_label(), _l2.get_label(), _l3.get_label()], frameon=False, fontsize=8.5, loc="upper left")
ax.set_title("Trend in cancer–childbearing age overlap, England 2001–2023")
fig.tight_layout(); fig.savefig(OUT_DIR/"figF_overlap_trend.png"); plt.close(fig)

# Fig G — reproductive-age counts
fig, ax = plt.subplots(figsize=(7.6,4.4))
ax.plot(YEARS, cnt_u50, "-o", color=ACC2, lw=2, label="Diagnosed under 50")
ax.plot(YEARS, cnt_u45, "-s", color=ACC, lw=2, label="Diagnosed under 45")
ax.plot(YEARS, cnt_u40, "-^", color=INK, lw=2, label="Diagnosed under 40")
ax.set_xlabel("Diagnosis year"); ax.set_ylabel("Annual registrations (England)")
ax.set_title("Uterine cancer (C54–C55) in reproductive-age women, England 2001–2023")
ax.legend(frameon=False, loc="upper left"); fig.tight_layout(); fig.savefig(OUT_DIR/"figG_reproductive_age_counts.png"); plt.close(fig)

# Fig H — projection: mean maternal age (left axis) + share of births to women 35+ (right axis)
fig, ax = plt.subplots(figsize=(7.6,4.4))
_l1, = ax.plot(proj_years_int, proj_mage, color=ACC, lw=2, label="Mean age of mother (left axis)")
ax.set_xlabel("Year"); ax.set_ylabel("Projected mean age of mother (years)", color=ACC)
ax.tick_params(axis="y", colors=ACC); ax.spines["left"].set_color(ACC)
ax2=ax.twinx(); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_color(INK)
_l2, = ax2.plot(proj_years_int, [100*s for s in proj_share35], color=INK, lw=2, ls="--", label="% of births to women aged 35+ (right axis)")
ax2.set_ylabel("Projected % of births to women aged 35+", color=INK); ax2.tick_params(axis="y", colors=INK)
ax.legend([_l1, _l2], [_l1.get_label(), _l2.get_label()], frameon=False, loc="upper left")
ax.set_title("Projected continuation of delayed childbearing, England, to 2047\n(ONS 2022-based projections)")
fig.tight_layout(); fig.savefig(OUT_DIR/"figH_projection_demographic.png"); plt.close(fig)

# Fig I — projected number of births to women in the EC-risk window (35+, 40+)
fig, ax = plt.subplots(figsize=(7.6,4.4))
ax.plot(proj_years_int, proj_cnt35, "-o", color=ACC2, lw=2, label="Births to women aged 35+")
ax.plot(proj_years_int, proj_cnt40, "-s", color=ACC, lw=2, label="Births to women aged 40+")
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.set_xlabel("Year"); ax.set_ylabel("Projected annual births (England)")
ax.set_title("Projected childbearing in the endometrial-cancer risk window, England")
ax.legend(frameon=False, loc="upper left"); fig.tight_layout(); fig.savefig(OUT_DIR/"figI_projection_counts.png"); plt.close(fig)

# %% [markdown]
# ## 8. Save results table

# %%
res_df = pd.DataFrame({"metric": list(results.keys()), "value": list(results.values())})
res_df.to_csv(OUT_DIR/"results_all.csv", index=False)
print(res_df.to_string(index=False))
print("\nFigures + results_all.csv written to", OUT_DIR.resolve())
