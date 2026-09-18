# =============================================================================
# Reproducible analysis code (v2.1.2; analysis identical to the v2.1.0 data-repair release)
# "Spatial Clustering of Hepatitis A in South Korea, 2020-2024: A Nationwide
#  Bayesian Analysis of Groundwater, Land Cover, and Socioeconomic Gradients"
# Seongdae Kim, Byung Chul Chun.
# Prepared for submission to JMIR Public Health and Surveillance (Original Paper).
# Releases up to v2.0.6 accompanied an earlier version of this analysis under the
# title "Sanitation Infrastructure and Environmental Vulnerability of Hepatitis A
# Transmission in South Korea, 2020-2024".
#
# WHY v2.1: an audit of the v2.0.x compiled dataset found that several
#  administrative source files store unavailable values as 0 rather than NA and
#  that some series were redefined in particular years. Under binary/quantile
#  coding these values made covariates encode data availability by year or region
#  instead of their true level. Section [2b] repairs this, and every repaired value
#  is logged in results/data_repair_log.csv. The v2.0.x headline associations for
#  piped water-supply coverage and inpatient medical cost were artefacts of this
#  problem and are NOT reproduced here. See README.md and DATA_PROVENANCE.md.
#
# Model: Bayesian negative-binomial disease mapping with a Besag-York-Mollie
#  (BYM) convolution + first-order temporal random walk (RW1) + Knorr-Held
#  Type I space-time interaction, fitted by INLA (R-INLA). The contiguity graph
#  has 223 districts; 220 districts (1,100 district-years, 2020-2024) have
#  complete covariates and enter the likelihood. 27 final covariates.
#
# The 27-covariate specification was developed through exploratory model
#  building and is not prospectively pre-specified; all credible associations
#  are hypothesis-generating.
#
# One run regenerates every number used by the manuscript:
#  [6]  principal model M6 and Table 2 (4 credible covariates)
#  [7]  model comparison M1-M6 with DIC, WAIC and effective parameters (Table S3)
#  [8]  global Moran's I, high/low-risk districts among analysed districts (Table S4)
#  [9]  8-graph neighbourhood sensitivity (Table S6)
#  [10] Getis-Ord Gi* on all 223 graph districts (Figure S2)
#  [11] alternative specifications, incl. oyster-coverage and extreme-value sensitivity (Tables S7 and S8)
#  [12] predictive diagnostics (CPO / PIT), district effects (Table S5)
#  [13] covariate estimates under nested/alternative random-effect structures
#  [14] numerical stability: INLA variants, glmmTMB cross-fit, repeated fits (Table S8)
#
# Software: R 4.6.0 and R-INLA 25.10.19 (results/sessionInfo.txt). INLA fits of
#  this model are not bit-reproducible: 20 repeated fits of M6 with identical
#  settings differ by about 3 DIC points and in the 4th decimal of the IRRs
#  (results/repeated_fits.csv); the credible set is identical in every fit.
# Run:  HAV_DATA_DIR=/path/to/data RUN_EXTENDED=true RUN_STABILITY=true Rscript HAV_spatial_reproducible.R
#
# DATA AVAILABILITY: annual district-level HAV notifications are released by the
#  Korea Disease Control and Prevention Agency (KDCA) Infectious Disease Portal
#  (https://dportal.kdca.go.kr); covariates come from KOSIS and the open-data
#  portals of the relevant Korean ministries. Raw source extracts are NOT
#  redistributed here; the compiled 1,100-row district-year analytic table is
#  provided in results/analysis_dataset_compiled.csv. Place the raw input
#  files under ./data (or set the HAV_DATA_DIR environment variable) to rebuild
#  it from source. No personally identifiable information is used
#  (aggregated district-year counts only).
# License: MIT (see LICENSE).
# =============================================================================

# ---------------------------------------------------------------------------
# [0] Packages  (INLA is installed from its own repository, not CRAN)
# ---------------------------------------------------------------------------
local({
  rp <- getOption("repos")
  if (is.null(rp) || is.na(rp["CRAN"]) || rp["CRAN"] %in% c("@CRAN@", ""))
    options(repos = c(CRAN = "https://cloud.r-project.org"))
})
need <- c("dplyr", "tidyr", "stringr", "car", "arrow", "sf", "spdep")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("INLA", quietly = TRUE))
  install.packages("INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/stable"),
    dependencies = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(car)
  library(arrow); library(sf); library(spdep); library(INLA)
})
options(scipen = 999)
set.seed(2024)

# ---------------------------------------------------------------------------
# [1] Configuration and input paths
# ---------------------------------------------------------------------------
DISEASE_NAME <- "A형간염"          # "Hepatitis A" label in the surveillance file
YEAR_START   <- 2020
YEAR_END     <- 2024
VIF_THRESHOLD <- 10            # collinearity screen (forced confounders are never dropped)
MIN_OBS      <- 20             # minimum non-missing district-years to use a covariate

# Input directory: ./data by default, override with HAV_DATA_DIR.
BASE_IV <- Sys.getenv("HAV_DATA_DIR", unset = file.path(getwd(), "data"))
OUT_DIR <- Sys.getenv("HAV_OUTPUT_DIR", unset = file.path(getwd(), "results"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pick_input <- function(...) {
  candidates <- file.path(BASE_IV, c(...))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) candidates[1] else hit[1]
}
PATH_DISEASE   <- pick_input("foodborne_final.csv", "식중독최종.csv")
PATH_HEALTH_PQ <- pick_input("health_indicators.parquet", "국민건강결과_최종.parquet")
PATH_SHP       <- pick_input("districts.shp", "final.shp")
if (!file.exists(PATH_DISEASE))
  stop("Input data not found under '", BASE_IV,
       "'. Set HAV_DATA_DIR to the folder holding the KDCA/KOSIS inputs.")

# ---------------------------------------------------------------------------
# [2] Helper functions
# ---------------------------------------------------------------------------
# Harmonise district keys across sources. Two administrative reorganisations
# need explicit handling so that left-joins do not silently drop districts:
#  * Incheon Michuhol-gu was renamed from Nam-gu.
#  * Sejong and Gunwi-gun use inconsistent province prefixes across files
#    (Gunwi was transferred from North Gyeongsang to Daegu in 2023). Without
#    this harmonisation the single-person-elderly covariate fails to join and
#    those district-years are lost to listwise deletion.
clean_region <- function(df) df %>% mutate(
  region = str_replace_all(as.character(region), "\\s+", ""),
  region = if_else(region == "인천시미추홀구", "인천시남구", region),
  region = if_else(region == "세종시",           "세종시세종시", region),
  region = if_else(region == "경상북도군위군", "대구시군위군", region),
  year   = as.integer(year)) %>%
  filter(year >= YEAR_START, year <= YEAR_END)

# ---------------------------------------------------------------------------
# [2b] v2.1 data repair (2026-09-18)
# ---------------------------------------------------------------------------
# Audit of the v2.0.x compiled dataset found that several source extracts store
# unavailable values as 0 rather than NA:
#  * piped-water coverage is 0 for 76 metropolitan/Incheon/Ulsan/Jeju districts
#    from 2018 onward (reported only at city level after 2017);
#  * whole survey years are 0 for food security (2020-2023), basic-livelihood
#    recipients and wastewater-discharge facilities (2021-2023; 2024 absent), elderly
#    living alone (2021) and groundwater-quality tests (every year before 2023; 2024 absent);
#  * inpatient medical cost is 0 in 85 districts in 2023 (4 in each of 2020-2022) and absent for 2024 in all but 4 districts;
#  * binary coding silently turned NA into 0 (groundwater wells 2023-2024).
# Repair rules, applied to the full 2008-2024 history before restricting to the study years:
#  (R1) rate/percentage/cost covariates cannot be 0 -> any 0 is set to NA;
#  (R1b) structural zeros of in-district inpatient cost (no inpatient facility) are kept as 0;
#  (R2) count covariates: a year in which >= 90% of districts are 0 (while other
#       years are mostly non-zero) is treated as unavailable -> NA; genuine
#       district-level zero counts are kept;
#  (R3) percentages > 100 are impossible -> NA;
#  (R4) a year whose national median departs > 1.6-fold from neighbouring years is a redefined
#       series (e.g., elderly living alone in 2024: median 59.5 vs 23.6 in 2023) -> that year NA;
#  (R5) NA in a study year is filled with the same district's value from the
#       nearest available year, preferring the most recent earlier year, then
#       the nearest later year;
#  (R6) functional-form coding preserves NA (no silent NA -> 0).
RATE_VARS <- c("상수도보급률", "독거노인비율", "기초생활수급자수율", "재정자립도", "재정자주도",
  "성비", "고령인구비율", "도시지역인구비율", "관내진료비_입원", "건강생활실천율_조율",
  "우울감경험률_표준화율", "식품안정성확보율_표준화율", "공공하수처리구역인구보급률(%)")
PCT_VARS <- c("상수도보급률", "독거노인비율", "기초생활수급자수율", "재정자립도", "재정자주도", "고령인구비율",
  "도시지역인구비율", "건강생활실천율_조율", "우울감경험률_표준화율", "식품안정성확보율_표준화율", "공공하수처리구역인구보급률(%)")
STRUCT_ZERO_VARS <- c("관내진료비_입원")
REPAIR_LOG <- data.frame()
clean_region_all <- function(df) df %>% mutate(
  region = str_replace_all(as.character(region), "\\s+", ""),
  region = if_else(region == "인천시미추홀구", "인천시남구", region),
  region = if_else(region == "세종시",           "세종시세종시", region),
  region = if_else(region == "경상북도군위군", "대구시군위군", region),
  year   = as.integer(year))
nearest_fill <- function(yr, x) {
  ok <- which(!is.na(x)); out <- x; src <- rep(NA_integer_, length(x))
  idx <- which(is.na(x) & yr >= YEAR_START & yr <= YEAR_END)
  if (!length(ok)) return(list(x = out, src = src))
  for (i in idx) {
    pv <- ok[yr[ok] < yr[i]]
    j <- if (length(pv)) pv[which.max(yr[pv])] else { nx <- ok[yr[ok] > yr[i]]; if (length(nx)) nx[which.min(yr[nx])] else NA }
    if (!is.na(j)) { out[i] <- x[j]; src[i] <- yr[j] }
  }
  list(x = out, src = src)
}
repair_and_fill <- function(df, vars, label) {
  vars <- intersect(vars, names(df)); if (!length(vars)) return(df)
  for (v in vars) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
  df <- df %>% group_by(region, year) %>%
    summarise(across(all_of(vars), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    tidyr::complete(region, year = union(unique(year), YEAR_START:YEAR_END)) %>% arrange(region, year)
  for (v in vars) {
    x <- df[[v]]; st <- df$year >= YEAR_START & df$year <= YEAR_END
    if (v %in% RATE_VARS) {
      bad <- !is.na(x) & x == 0
      # (R1b) structural zeros: in-district inpatient cost is genuinely 0 where no inpatient facility
      #       exists. A district whose value is 0 in every year but at most one is kept as 0.
      if (v %in% STRUCT_ZERO_VARS) {
        npos <- tapply(!is.na(x) & x > 0, df$region, sum)
        sz <- names(npos)[npos <= 1]
        x[df$region %in% sz & !is.na(x) & x > 0] <- 0      # the single stray positive is an artefact year
        bad <- bad & !(df$region %in% sz)
      }
    } else {
      nobs <- tapply(!is.na(x), df$year, sum)
      yz <- (tapply(!is.na(x) & x == 0, df$year, sum) / pmax(nobs, 1))[nobs > 0]   # years with data only
      byr <- as.integer(names(yz)[yz >= 0.9 & !is.na(yz)])
      bad <- if (any(yz < 0.5, na.rm = TRUE)) df$year %in% byr & !is.na(x) & x == 0 else rep(FALSE, length(x))
    }
    x[bad] <- NA
    # (R3) percentages cannot exceed 100 -> NA (e.g., Nonsan elderly-living-alone 135.4 in 2020)
    n_over <- 0L
    if (v %in% PCT_VARS) { over <- !is.na(x) & x > 100; n_over <- sum(over & st); x[over] <- NA }
    # (R4) year-level comparability guard for rate/percentage covariates: if a year's national
    #      median departs from the median of its (up to 4) nearest other years by a factor > 1.6,
    #      the series was redefined in that year -> the whole year is treated as non-comparable (NA).
    n_shift <- 0L; shift_years <- integer(0)
    if (v %in% RATE_VARS) {
      ym <- tapply(x, df$year, function(z) { z <- z[!is.na(z)]; if (length(z) >= 20) median(z) else NA_real_ })
      ym <- ym[!is.na(ym)]; yy <- as.integer(names(ym))
      for (k in seq_along(ym)) {
        oth <- order(abs(yy - yy[k]))[-1]; oth <- oth[seq_len(min(4, length(oth)))]
        if (length(oth) >= 2) { ref <- median(ym[oth]); if (is.finite(ref) && ref > 0 && (ym[k] / ref > 1.6 || ym[k] / ref < 1 / 1.6)) shift_years <- c(shift_years, yy[k]) }
      }
      if (length(shift_years)) { sh <- df$year %in% shift_years & !is.na(x); n_shift <- sum(sh & st); x[sh] <- NA }
    }
    na_before <- sum(is.na(x[st]))
    res <- lapply(split(seq_along(x), df$region), function(ii) { r <- nearest_fill(df$year[ii], x[ii]); list(ii = ii, x = r$x) })
    for (r in res) x[r$ii] <- r$x
    df[[v]] <- x
    REPAIR_LOG <<- rbind(REPAIR_LOG, data.frame(source = label, variable = v,
      zeros_set_NA_study_years = sum(bad & st), over100_set_NA_study_years = n_over,
      noncomparable_year_set_NA_study_years = n_shift, noncomparable_years = paste(shift_years[shift_years >= YEAR_START & shift_years <= YEAR_END], collapse = ";"),
      NA_study_years_before_fill = na_before,
      NA_study_years_after_fill = sum(is.na(x[st]))))
  }
  df
}

read_csv_safe <- function(fp) {
  for (enc in c("UTF-8", "UTF-8-BOM", "CP949", "EUC-KR")) {
    raw <- tryCatch(read.csv(fp, fileEncoding = enc, check.names = FALSE,
                             stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) return(raw)
  }
  NULL
}

is_pct <- function(x) { xv <- x[!is.na(x) & is.finite(x)]; all(xv >= 0 & xv <= 100) & max(xv) > 1 }

# ---------------------------------------------------------------------------
# [3] Data assembly
# ---------------------------------------------------------------------------
cat("## [3] Loading data\n")
df_raw <- read.csv(PATH_DISEASE, stringsAsFactors = FALSE, check.names = FALSE)
df_target <- df_raw %>%
  filter(disease == DISEASE_NAME, year >= YEAR_START, year <= YEAR_END) %>%
  clean_region() %>%
  group_by(region, year) %>%
  summarise(cases = sum(cases, na.rm = TRUE),
            population = mean(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(rate_100k = cases / population * 1e5)
cat(sprintf("  HAV: %d district-years | %d districts | %d cases\n",
            nrow(df_target), n_distinct(df_target$region), sum(df_target$cases)))
cor_merged <- df_target

# (A) community-health indicators (parquet)
health_vars <- c("건강생활실천율_조율", "상수도보급률", "독거노인비율",
  "기초생활수급자수율", "재정자립도", "재정자주도", "성비", "고령인구비율",
  "도시지역인구비율", "순이동인구", "관내진료비_입원", "의원_가정의학과",
  "우울감경험률_표준화율", "식품안정성확보율_표준화율", "폐수배출업소수")
tryCatch({
  hpq <- read_parquet(PATH_HEALTH_PQ) %>% as.data.frame() %>% clean_region_all()
  ah <- intersect(health_vars, names(hpq))
  hagg <- repair_and_fill(hpq[, c("region", "year", ah)], ah, "health_parquet") %>%
    filter(year >= YEAR_START, year <= YEAR_END)
  cor_merged <- cor_merged %>% left_join(hagg, by = c("region", "year"))
  cat(sprintf("  parquet health indicators: %d variables\n", length(ah)))
}, error = function(e) cat(sprintf("  [warn] parquet: %s\n", e$message)))

# (B) covariates supplied as individual CSV extracts. Only the columns used by
#     the final reported model are kept; gaps are filled by repair_and_fill() (rule R5).
selected <- list(
  "groundwater_household.csv|생활용지하수이용현황_전처리.csv" = c("가정용_개소수", "간이상수도용_개소수"),
  "groundwater_quality.csv|merged_지하수수질.csv" = c("검사합계"),
  "sewer_repair.csv|merged_하수관로개보수.csv" = c("개·보수관로_부분보수(개소)_계"),
  "sewerage_coverage.csv|merged_하수도보급률.csv" = c("공공하수처리구역인구보급률(%)", "하수처리구역외_정화조인구", "총면적(㎢)"),
  "livestock.csv|가축두수_전처리.csv" = c("농가수(호)_젖소", "농가수(호)_돼지", "농가수(호)_가금"),
  "elderly_singleperson.csv|고령인구_전처리.csv" = c("1인가구_80~84세"),
  "land_use.csv|국토이용현황_전처리_수정.csv" = c("답", "임야", "대"),
  "shellfish.csv|어패류_패류_전처리.csv" = c("굴_자연채묘 생산량(kg)"))
for (fn in names(selected)) {
  aliases <- strsplit(fn, "|", fixed = TRUE)[[1]]
  fp <- pick_input(aliases); if (!file.exists(fp)) next
  raw <- read_csv_safe(fp); if (is.null(raw)) next
  raw <- raw %>% clean_region_all()
  av <- intersect(selected[[fn]], names(raw)); if (length(av) == 0) next
  agg <- repair_and_fill(raw[, c("region", "year", av)], av, basename(fp)) %>%
    filter(year >= YEAR_START, year <= YEAR_END)
  cor_merged <- cor_merged %>% left_join(agg, by = c("region", "year"))
  cat(sprintf("  csv %-30s %d variables\n", basename(fp), length(av)))
}
cat(sprintf("  merged: %d rows x %d columns\n", nrow(cor_merged), ncol(cor_merged)))
write.csv(REPAIR_LOG, file.path(OUT_DIR, "data_repair_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("  data repair log:\n"); print(REPAIR_LOG[REPAIR_LOG$zeros_set_NA_study_years > 0 | REPAIR_LOG$NA_study_years_before_fill > 0, c("variable","zeros_set_NA_study_years","over100_set_NA_study_years","noncomparable_years","NA_study_years_after_fill")], row.names = FALSE)

# ---------------------------------------------------------------------------
# [4] Final selected covariate set (27) and functional forms
#     (Table S1, Multimedia Appendix 1)
# ---------------------------------------------------------------------------
# Each covariate enters here with the final reported functional form:
#   raw    = standardised continuous
#   log1p  = log(1 + x), standardised
#   binary = above median (or non-zero for zero-inflated counts)
#   T3/Q4  = ordered tertile / quartile class, entered as a standardised score
# Four demographic/urbanisation covariates (sex_ratio, elderly_rate,
# urban_pop_rate, net_migration) are FORCED confounders: kept in every model
# and never removed by the collinearity screen.
TV <- data.frame(
  code = c("굴_자연채묘 생산량(kg)", "농가수(호)_젖소",
    "가정용_개소수", "간이상수도용_개소수", "상수도보급률", "하수처리구역외_정화조인구",
    "폐수배출업소수", "공공하수처리구역인구보급률(%)", "검사합계", "개·보수관로_부분보수(개소)_계",
    "답", "임야", "대", "건강생활실천율_조율", "식품안정성확보율_표준화율",
    "독거노인비율", "1인가구_80~84세", "기초생활수급자수율", "재정자립도", "재정자주도",
    "성비", "고령인구비율", "도시지역인구비율", "순이동인구",
    "관내진료비_입원", "의원_가정의학과", "우울감경험률_표준화율"),
  eng = c("oyster", "dairy_farm",
    "gw_household", "gw_simple", "water_supply", "septic_pop",
    "ww_facility", "pub_sewage", "test_total", "sewer_repair",
    "paddy", "forest", "residential", "health_practice", "food_safety",
    "elderly_alone", "alone_80_84", "welfare", "fiscal_indep", "fiscal_auto",
    "sex_ratio", "elderly_rate", "urban_pop_rate", "net_migration",
    "med_in", "clinic_family", "depression"),
  form = c("Q4", "raw",
    "binary", "binary", "Q4", "raw",
    "T3", "Q4", "binary", "raw",
    "raw", "T3", "Q4", "raw", "binary",
    "log1p", "raw", "binary", "raw", "raw",
    "raw", "Q4", "T3", "log1p",
    "log1p", "Q4", "binary"),
  forced = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "Y", "Y", "Y", "Y", "", "", ""),
  stringsAsFactors = FALSE)
cat(sprintf("\n## [4] Final selected covariates: %d (forced confounders: %d)\n",
            nrow(TV), sum(TV$forced == "Y")))

# District polygons and contiguity graph (principal neighbourhood structure).
# Six island districts with no land contiguity are excluded so that the spatial
# graph has 223 districts in 2 components (mainland; Jeju-si + Seogwipo-si).
shp <- st_read(PATH_SHP, quiet = TRUE) %>%
  mutate(region = str_replace_all(as.character(region), "\\s+", ""),
         region = if_else(region == "인천시미추홀구", "인천시남구", region))
islands <- c("인천시옹진군", "전라남도완도군", "전라남도진도군",
             "경상남도거제시", "경상남도남해군", "경상북도울릉군")
shp_main <- shp %>% filter(!region %in% islands)
nb_obj <- poly2nb(shp_main, snap = 0.01); iso <- which(card(nb_obj) == 0)
if (length(iso) > 0) { shp_main <- shp_main[-iso, ]; nb_obj <- poly2nb(shp_main, snap = 0.01) }
graph_file <- tempfile(fileext = ".graph")
nb2INLA(graph_file, nb_obj); g_main <- inla.read.graph(graph_file)
nb_w <- nb2listw(nb_obj, style = "W", zero.policy = TRUE)
cat(sprintf("  districts in spatial model: %d\n", nrow(shp_main)))

# ---------------------------------------------------------------------------
# [5] Build the design matrix from the final reported forms
# ---------------------------------------------------------------------------
# apply_form() materialises one covariate in its declared functional form.
# `hz` flags zero-inflated counts (>20% zeros), for which the binary/tertile
# cut-points use presence/non-zero medians rather than the overall median.
apply_form <- function(x, form, hz) {
  if (form == "raw")   return(x)
  if (form == "log1p") return(log1p(pmax(x, 0)))
  if (form == "binary") {
    if (hz) return(ifelse(is.na(x), NA_real_, as.numeric(x > 0)))
    md <- median(x, na.rm = TRUE); return(ifelse(is.na(x), NA_real_, as.numeric(x > md)))
  }
  if (form == "T3") {
    if (hz) {
      nz <- x[!is.na(x) & x > 0]; mn <- median(nz, na.rm = TRUE)
      return(dplyr::case_when(is.na(x) ~ NA_real_, x == 0 ~ 1, x <= mn ~ 2, x > mn ~ 3))
    }
    q33 <- quantile(x, c(1/3, 2/3), na.rm = TRUE)
    return(as.numeric(cut(x, unique(c(-Inf, q33[1], q33[2], Inf)),
                          labels = FALSE, include.lowest = TRUE)))
  }
  if (form == "Q4") {
    q4 <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
    return(as.numeric(cut(x, unique(c(-Inf, q4[1], q4[2], q4[3], Inf)),
                          labels = FALSE, include.lowest = TRUE)))
  }
  stop("unknown form: ", form)
}

df_w <- cor_merged %>% filter(population > 0, region %in% shp_main$region)
TV <- TV[TV$code %in% names(df_w), ]
data_ext <- df_w
zcols <- character(0); DICT <- list()
for (i in seq_len(nrow(TV))) {
  code <- TV$code[i]; x <- as.numeric(df_w[[code]])
  nv <- sum(!is.na(x) & is.finite(x)); if (nv < MIN_OBS) next
  hz <- sum(!is.na(x) & is.finite(x) & x == 0) / nv * 100 > 20
  form <- if (TV$eng[i] == "sex_ratio") "raw" else TV$form[i]   # sex ratio is symmetric -> raw
  val <- apply_form(x, form, hz)
  xs <- x[!is.na(x) & is.finite(x)]
  cut_txt <- switch(form,
    raw = "", log1p = "log(1 + max(x, 0))",
    binary = if (hz) "x > 0" else sprintf("x > %s (median)", signif(median(xs), 6)),
    T3 = if (hz) sprintf("0 | (0, %s] | > %s (median of non-zero values)", signif(median(xs[xs > 0]), 6), signif(median(xs[xs > 0]), 6))
         else paste(signif(quantile(xs, c(1/3, 2/3)), 6), collapse = " | "),
    Q4 = paste(signif(quantile(xs, c(.25, .5, .75)), 6), collapse = " | "))
  DICT[[length(DICT) + 1]] <- data.frame(covariate = TV$eng[i], source_column = code, form = form, zero_inflated = hz,
    n_levels = length(unique(val[!is.na(val)])), cut_points = cut_txt, n_negative = sum(xs < 0), stringsAsFactors = FALSE)
  s <- sd(val, na.rm = TRUE); m <- mean(val, na.rm = TRUE)
  zname <- paste0(TV$eng[i], "_z")
  data_ext[[zname]] <- if (!is.na(s) && s > 0) (val - m) / s else val
  zcols <- c(zcols, zname)
}

# Collinearity screen: drop covariates with VIF > threshold one at a time,
# but never drop a forced confounder.
forced_z <- paste0(TV$eng[TV$forced == "Y"], "_z")
vif_data <- data_ext[, c("cases", zcols), drop = FALSE]
vif_data <- vif_data[complete.cases(vif_data), ]
keep <- zcols
for (step in 1:40) {
  if (length(keep) <= 1) break
  lm_t <- tryCatch(lm(as.formula(paste("cases ~", paste0("`", keep, "`", collapse = "+"))),
                      data = vif_data), error = function(e) NULL)
  if (is.null(lm_t)) break
  vv <- tryCatch(car::vif(lm_t), error = function(e) NULL); if (is.null(vv)) break
  names(vv) <- gsub("`", "", names(vv))
  if (max(vv, na.rm = TRUE) < VIF_THRESHOLD) break
  removable <- vv[!names(vv) %in% forced_z]
  if (!length(removable) || max(removable, na.rm = TRUE) < VIF_THRESHOLD) break
  drop <- names(which.max(removable))
  keep <- setdiff(keep, drop)
}
covs <- keep
final_vif <- tryCatch({
  mm <- lm(as.formula(paste("cases ~", paste0("`", covs, "`", collapse = "+"))),
           data = vif_data)
  car::vif(mm)
}, error = function(e) NULL)
if (!is.null(final_vif))
  write.csv(data.frame(variable = names(final_vif), VIF = as.numeric(final_vif)),
            file.path(OUT_DIR, "vif_audit.csv"), row.names = FALSE)
cat(sprintf("  covariates entering INLA: %d (VIF < %d)\n", length(covs), VIF_THRESHOLD))

# Final analysis frame: complete cases on the modelled covariates, indexed by
# district (idarea) and year (idtime) with a space-time interaction index.
rmap <- data.frame(region = shp_main$region, idarea = seq_len(nrow(shp_main)))
ymap <- data.frame(year = YEAR_START:YEAR_END, idtime = seq_along(YEAR_START:YEAR_END))
ic <- data_ext[complete.cases(data_ext[, covs]), ] %>%
  left_join(rmap, by = "region") %>% left_join(ymap, by = "year") %>%
  arrange(idarea, idtime)
ic$idarea_time <- seq_len(nrow(ic))
cat(sprintf("  analysis frame: N = %d district-years | EPV = %.1f\n",
            nrow(ic), nrow(ic) / length(covs)))
write.csv(ic, file.path(OUT_DIR, "analysis_dataset_compiled.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
# Covariate dictionary: applied functional form, cut points and descriptive statistics of the raw values in the analysis frame
DICT <- do.call(rbind, DICT); DICT <- DICT[paste0(DICT$covariate, "_z") %in% covs, ]
DICT <- cbind(DICT, do.call(rbind, lapply(DICT$source_column, function(cc) { v <- as.numeric(ic[[cc]])
  data.frame(n = sum(!is.na(v)), mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), min = min(v, na.rm = TRUE),
             q1 = unname(quantile(v, .25, na.rm = TRUE)), median = median(v, na.rm = TRUE), q3 = unname(quantile(v, .75, na.rm = TRUE)),
             max = max(v, na.rm = TRUE), share_zero = mean(v == 0, na.rm = TRUE),
             districts_constant_over_years = sum(tapply(v, ic$region, function(z) length(unique(z)) == 1))) })))
write.csv(DICT, file.path(OUT_DIR, "covariate_dictionary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Priors (penalised-complexity) shared across models.
pc_bym  <- list(prec.unstruct = list(prior = "pc.prec", param = c(0.5, 0.01)),
                prec.spatial  = list(prior = "pc.prec", param = c(0.5, 0.01)))
pc_prec <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
cov_str <- paste(covs, collapse = " + ")
base_f  <- paste("cases ~", cov_str, "+ offset(log(population + 1))")
fitm <- function(fs) tryCatch(
  inla(as.formula(fs), family = "nbinomial", data = ic,
       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
       control.predictor = list(link = 1)), error = function(e) { message(e$message); NULL })

# ---------------------------------------------------------------------------
# [6] Principal model M6 (BYM + RW1 + Type I interaction) and Table 2
# ---------------------------------------------------------------------------
RE_FULL <- paste("+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)",
                 "+ f(idtime, model='rw1', hyper=pc_prec)",
                 "+ f(idarea_time, model='iid', hyper=pc_prec)")
M6 <- fitm(paste(base_f, RE_FULL))
cat(sprintf("\n## [6] Principal model M6: DIC = %.2f | WAIC = %.2f\n",
            M6$dic$dic, M6$waic$waic))
fe <- M6$summary.fixed; fe <- fe[rownames(fe) != "(Intercept)", , drop = FALSE]
table2 <- data.frame(
  covariate = gsub("_z$", "", rownames(fe)),
  IRR = round(exp(fe$mean), 6),             # 6 decimals: rounding to 3 here would double-round the 2-decimal values of the paper
  lo  = round(exp(fe$`0.025quant`), 6),
  hi  = round(exp(fe$`0.975quant`), 6),
  credible = as.integer(fe$`0.025quant` > 0 | fe$`0.975quant` < 0),
  row.names = NULL)
cat(sprintf("  Table 2: %d covariates, %d credible (95%% CrI excludes 1)\n",
            nrow(table2), sum(table2$credible)))
print(table2[table2$credible == 1, ], row.names = FALSE)

FAST_PRINCIPAL <- identical(tolower(Sys.getenv("FAST_PRINCIPAL", "false")), "true")
if (FAST_PRINCIPAL) {
  agg_fast <- ic %>% group_by(idarea) %>%
    summarise(rate = sum(cases) / sum(population) * 1e5, .groups = "drop") %>% arrange(idarea)
  rv_fast <- rep(NA, nrow(shp_main)); rv_fast[agg_fast$idarea] <- agg_fast$rate
  moran_pre_fast <- moran.test(rv_fast, nb_w, zero.policy = TRUE, na.action = na.omit)
  fit_fast <- M6$summary.fitted.values$mean[seq_len(nrow(ic))]
  res_fast <- ic %>% mutate(resid = cases - fit_fast) %>% group_by(idarea) %>%
    summarise(r = sum(resid), .groups = "drop") %>% arrange(idarea)
  rr_fast <- rep(NA, nrow(shp_main)); rr_fast[res_fast$idarea] <- res_fast$r
  moran_post_fast <- moran.test(rr_fast, nb_w, zero.policy = TRUE, na.action = na.omit)
  re_fast <- M6$summary.random$idarea; na_fast <- nrow(shp_main)
  ana_fast <- sort(unique(ic$idarea))          # districts that enter the likelihood
  n_high_fast <- sum(re_fast$`0.025quant`[ana_fast] > 0)
  n_low_fast <- sum(re_fast$`0.975quant`[ana_fast] < 0)
  write.csv(table2, file.path(OUT_DIR, "table2_principal_IRR.csv"), row.names = FALSE)
  write.csv(data.frame(region = ic$region, year = ic$year,
                       cpo = M6$cpo$cpo, pit = M6$cpo$pit,
                       numerical_failure = M6$cpo$failure),
            file.path(OUT_DIR, "cpo_pit_diagnostics.csv"), row.names = FALSE)
  write.csv(data.frame(
    metric = c("N_district_years", "districts", "cases_all_229_districts",
               "M6_DIC", "M6_WAIC", "crude_Moran_I", "crude_Moran_p",
               "residual_Moran_I", "residual_Moran_p", "high_risk", "low_risk"),
    value = c(nrow(ic), nrow(shp_main), sum(df_target$cases), M6$dic$dic, M6$waic$waic,
              moran_pre_fast$estimate[[1]], moran_pre_fast$p.value,
              moran_post_fast$estimate[[1]], moran_post_fast$p.value,
              n_high_fast, n_low_fast)),
    file.path(OUT_DIR, "core_diagnostics.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
  cat("\nFAST_PRINCIPAL run complete.\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# [7] Model comparison M1-M6 (Table S3, Multimedia Appendix 1)
# ---------------------------------------------------------------------------
cat("\n## [7] Model comparison (Table S3, Multimedia Appendix 1)\n")
M <- list(
  M1 = base_f,
  M2 = paste(base_f, "+ f(idarea, model='besag', graph=g_main, scale.model=TRUE, hyper=pc_prec)"),
  M3 = paste(base_f, "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)"),
  M4 = paste(base_f, "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)",
             "+ f(idarea_time, model='iid', hyper=pc_prec)"),
  M5 = paste(base_f, "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)",
             "+ f(idtime, model='rw1', hyper=pc_prec)"),
  M6 = paste(base_f, RE_FULL))
tableS1 <- do.call(rbind, lapply(names(M), function(n) {
  f <- if (n == "M6") M6 else fitm(M[[n]])
  data.frame(model = n, DIC = round(f$dic$dic, 2), WAIC = round(f$waic$waic, 2), pD = round(f$dic$p.eff, 1))
}))
print(tableS1, row.names = FALSE)

# ---------------------------------------------------------------------------
# [8] Global Moran's I, pre- and post-modelling (Table S4, Multimedia Appendix 1)
# ---------------------------------------------------------------------------
cat("\n## [8] Global Moran's I (Table S4, Multimedia Appendix 1)\n")
agg <- ic %>% group_by(idarea) %>%
  summarise(rate = sum(cases) / sum(population) * 1e5, .groups = "drop") %>% arrange(idarea)
rv <- rep(NA, nrow(shp_main)); rv[agg$idarea] <- agg$rate
moran_pre <- moran.test(rv, nb_w, zero.policy = TRUE, na.action = na.omit)
ic$fitv <- M6$summary.fitted.values$mean[seq_len(nrow(ic))]
res_agg <- ic %>% mutate(resid = cases - fitv) %>% group_by(idarea) %>%
  summarise(r = sum(resid), .groups = "drop") %>% arrange(idarea)
rr <- rep(NA, nrow(shp_main)); rr[res_agg$idarea] <- res_agg$r
moran_post <- moran.test(rr, nb_w, zero.policy = TRUE, na.action = na.omit)
# District-level Pearson residual: (sum y - sum mu) / sqrt(sum Var), Var = mu + mu^2 / size (negative binomial).
# Fitted values include the random effects, so this measures autocorrelation left after the whole model.
nb_size <- M6$summary.hyperpar[grep("size for the nbinomial", rownames(M6$summary.hyperpar)), "mean"]
pr_agg <- ic %>% mutate(v = fitv + fitv^2 / nb_size) %>% group_by(idarea) %>%
  summarise(pr = (sum(cases) - sum(fitv)) / sqrt(sum(v)), .groups = "drop") %>% arrange(idarea)
rp <- rep(NA, nrow(shp_main)); rp[pr_agg$idarea] <- pr_agg$pr
moran_pearson <- moran.test(rp, nb_w, zero.policy = TRUE, na.action = na.omit)
cat(sprintf("  Pearson-residual I = %+.4f (p = %.3g)\n", moran_pearson$estimate[[1]], moran_pearson$p.value))
cat(sprintf("  crude    I = %+.4f (p = %.3g)\n", moran_pre$estimate[[1]], moran_pre$p.value))
cat(sprintf("  residual I = %+.4f (p = %.3g)\n", moran_post$estimate[[1]], moran_post$p.value))

# High- and low-risk districts: combined district effect (u + v of the BYM term) whose 95% CrI
# excludes zero.
# Only districts that enter the likelihood are classified. The 3 graph districts without
# covariate data receive effects interpolated from the prior and their neighbours.
re <- M6$summary.random$idarea; na <- nrow(shp_main); ana <- sort(unique(ic$idarea))
n_high <- sum(re$`0.025quant`[ana] > 0); n_low <- sum(re$`0.975quant`[ana] < 0)
n_high_all <- sum(re$`0.025quant`[1:na] > 0); n_low_all <- sum(re$`0.975quant`[1:na] < 0)
cat(sprintf("  analysed districts = %d | high-risk = %d | low-risk = %d (all %d graph nodes: %d / %d)\n",
            length(ana), n_high, n_low, na, n_high_all, n_low_all))
# Hyperparameters on the SD scale (exact transform of the precision quantiles). The posterior MEAN of a
# precision is dominated by its right tail, so the median SD is the quantity to read.
hp6 <- M6$summary.hyperpar; hp6 <- hp6[grepl("^Precision", rownames(hp6)), ]
write.csv(data.frame(component = rownames(hp6), sd_median = 1 / sqrt(hp6$`0.5quant`), sd_lo = 1 / sqrt(hp6$`0.975quant`),
                     sd_hi = 1 / sqrt(hp6$`0.025quant`), precision_mean = hp6$mean, precision_mode = hp6$mode),
          file.path(OUT_DIR, "m6_hyperparameter_sd.csv"), row.names = FALSE)

# Persist the independently regenerated core results before optional extended
# analyses. This makes a failed sensitivity fit unable to erase the audit trail.
write.csv(table2, file.path(OUT_DIR, "table2_principal_IRR.csv"), row.names = FALSE)
write.csv(tableS1, file.path(OUT_DIR, "model_comparison.csv"), row.names = FALSE)
write.csv(data.frame(region = ic$region, year = ic$year,
                     cpo = M6$cpo$cpo, pit = M6$cpo$pit,
                     numerical_failure = M6$cpo$failure),
          file.path(OUT_DIR, "cpo_pit_diagnostics.csv"), row.names = FALSE)
write.csv(data.frame(
  metric = c("N_district_years", "districts", "cases_all_229_districts",
             "M6_DIC", "M6_WAIC", "crude_Moran_I", "crude_Moran_p",
             "residual_Moran_I", "residual_Moran_p", "high_risk", "low_risk",
             "analysed_districts", "M6_pD", "high_risk_all_graph_nodes", "low_risk_all_graph_nodes",
             "pearson_residual_Moran_I", "pearson_residual_Moran_p"),
  value = c(nrow(ic), nrow(shp_main), sum(df_target$cases), M6$dic$dic,
            M6$waic$waic, moran_pre$estimate[[1]], moran_pre$p.value,
            moran_post$estimate[[1]], moran_post$p.value, n_high, n_low,
            length(ana), M6$dic$p.eff, n_high_all, n_low_all,
            moran_pearson$estimate[[1]], moran_pearson$p.value)),
  file.path(OUT_DIR, "core_diagnostics.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

RUN_EXTENDED <- identical(tolower(Sys.getenv("RUN_EXTENDED", "false")), "true")
if (!RUN_EXTENDED) {
  cat("\nCore reproducibility run complete. Set RUN_EXTENDED=true for sections 9-14.\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# [9] Eight-graph neighbourhood sensitivity (Table S6, Multimedia Appendix 1)
# ---------------------------------------------------------------------------
cat("\n## [9] 8-graph neighbourhood sensitivity (Table S6, Multimedia Appendix 1)\n")
cz <- st_coordinates(st_centroid(st_geometry(shp_main)))
mkgraph <- function(nb) { f <- tempfile(); nb2INLA(f, nb); inla.read.graph(f) }
graphs <- list(Queen = poly2nb(shp_main, queen = TRUE,  snap = 0.01),
               Rook  = poly2nb(shp_main, queen = FALSE, snap = 0.01))
for (k in 2:7) graphs[[paste0("knn", k)]] <- make.sym.nb(knn2nb(knearneigh(cz, k = k)))
cred <- paste0(table2$covariate[table2$credible == 1], "_z")   # credible set of the refitted M6
cred <- intersect(cred, covs)
graph_irr <- list(); graph_cnt <- list()
graph_cred <- setNames(integer(length(cred)), cred)
for (gn in names(graphs)) {
  gg <- tryCatch(mkgraph(graphs[[gn]]), error = function(e) NULL); if (is.null(gg)) next
  ff <- fitm(paste(base_f,
    "+ f(idarea, model='bym', graph=gg, scale.model=TRUE, hyper=pc_bym)",
    "+ f(idtime, model='rw1', hyper=pc_prec)",
    "+ f(idarea_time, model='iid', hyper=pc_prec)"))
  if (is.null(ff)) next
  fe2 <- ff$summary.fixed
  graph_irr[[gn]] <- data.frame(graph = gn, covariate = rownames(fe2), IRR = exp(fe2$mean),
    lo = exp(fe2$`0.025quant`), hi = exp(fe2$`0.975quant`), DIC = ff$dic$dic,
    n_high = sum(ff$summary.random$idarea$`0.025quant`[1:nrow(shp_main)] > 0),
    n_low = sum(ff$summary.random$idarea$`0.975quant`[1:nrow(shp_main)] < 0))
  graph_cnt[[gn]] <- data.frame(graph = gn, n_nodes = nrow(shp_main), n_analysed = length(ana),
    n_high_all = sum(ff$summary.random$idarea$`0.025quant`[1:nrow(shp_main)] > 0),
    n_low_all = sum(ff$summary.random$idarea$`0.975quant`[1:nrow(shp_main)] < 0),
    n_high_220 = sum(ff$summary.random$idarea$`0.025quant`[ana] > 0),
    n_low_220 = sum(ff$summary.random$idarea$`0.975quant`[ana] < 0),
    mean_links = mean(card(graphs[[gn]])), DIC = ff$dic$dic)
  for (c in cred) if (c %in% rownames(fe2))
    graph_cred[c] <- graph_cred[c] +
      as.integer(fe2[c, "0.025quant"] > 0 | fe2[c, "0.975quant"] < 0)
  cat(sprintf("  %-6s DIC = %.2f\n", gn, ff$dic$dic))
}
cat("  credible across graphs (out of 8):\n")
for (c in cred) cat(sprintf("    %-18s %d/8\n", c, graph_cred[c]))
write.csv(do.call(rbind, graph_irr), file.path(OUT_DIR, "graph_sensitivity.csv"), row.names = FALSE)
# counts restricted to the analysed districts (column suffix _220 = the 220 districts with complete covariates)
write.csv(do.call(rbind, graph_cnt), file.path(OUT_DIR, "graph_counts_analysed220.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# [10] Getis-Ord Gi* local clustering (Figure S2, Multimedia Appendix 3)
# ---------------------------------------------------------------------------
cat("\n## [10] Getis-Ord Gi* (Figure S2, Multimedia Appendix 3)\n")
nb_self <- include.self(nb_obj)
lw_self <- nb2listw(nb_self, style = "B", zero.policy = TRUE)
# Descriptive Gi* uses the observed crude rate of ALL graph districts. Outcome data exist for every
# district; only covariates are missing for the 3 complete-case exclusions. (v2.1 fix: these 3 districts
# were previously entered as rate = 0, which produced spurious cold spots.)
rv_full <- data_ext %>% group_by(region) %>%
  summarise(rate = sum(cases) / sum(population) * 1e5, cases = sum(cases), .groups = "drop")
cases_223 <- sum(rv_full$cases[rv_full$region %in% shp_main$region])
rv_full <- rv_full$rate[match(shp_main$region, rv_full$region)]
stopifnot(!anyNA(rv_full))
moran_pre_full <- moran.test(rv_full, nb_w, zero.policy = TRUE)
write.csv(data.frame(metric = c("crude_Moran_I_223", "crude_Moran_p_223", "min_rate_223", "max_rate_223", "cases_223"),
                     value = c(moran_pre_full$estimate[[1]], moran_pre_full$p.value, min(rv_full), max(rv_full), cases_223)),
          file.path(OUT_DIR, "descriptive_223.csv"), row.names = FALSE)
gi <- localG(rv_full, lw_self, zero.policy = TRUE)
gi_z <- as.numeric(gi)
cat(sprintf("  Gi* z-scores: hot spots (z > 1.96) = %d | cold spots (z < -1.96) = %d\n",
            sum(gi_z > 1.96, na.rm = TRUE), sum(gi_z < -1.96, na.rm = TRUE)))
getis <- data.frame(region = shp_main$region, Gi_z = round(gi_z, 3))
write.csv(getis, file.path(OUT_DIR, "getis_ord_gi.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# [11] Alternative-specification robustness checks (Table S7, Multimedia Appendix 1)
# ---------------------------------------------------------------------------
cat("\n## [11] Robustness checks (Table S7, Multimedia Appendix 1)\n")
irr_txt <- function(fit, name) {
  if (!name %in% rownames(fit$summary.fixed)) return("-")
  r <- fit$summary.fixed[name, ]
  sprintf("%.6f (%.6f-%.6f)%s", exp(r$mean), exp(r$`0.025quant`), exp(r$`0.975quant`),
          ifelse(r$`0.025quant` > 0 | r$`0.975quant` < 0, "*", ""))
}
fit_alt <- function(dat, cv) tryCatch(
  inla(as.formula(paste("cases ~", paste(cv, collapse = " + "),
       "+ offset(log(population + 1))", RE_FULL)),
       family = "nbinomial", data = dat,
       control.compute = list(dic = TRUE, waic = TRUE),
       control.predictor = list(link = 1)), error = function(e) NULL)

# (a) exclude the COVID-19 years 2020-2021
fa <- fit_alt(ic %>% filter(!year %in% c(2020, 2021)), covs)
# (b) drop inpatient medical cost (guard against over-adjustment / ascertainment)
fb <- fit_alt(ic, setdiff(covs, "med_in_z"))
# (c) add swine and poultry farm counts (tests specificity of the dairy signal).
#     Source columns are 농가수(호) = number of farms, not animal headcount or area density.
add_cov <- function(dat, raw_col, newname) {
  key <- cor_merged %>% transmute(region, year, val = suppressWarnings(as.numeric(.data[[raw_col]])))
  m <- dat %>% left_join(key, by = c("region", "year"))
  v <- m$val; v[!is.finite(v)] <- NA
  z <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  if (anyNA(z)) stop("add_cov: ", sum(is.na(z)), " missing values in ", raw_col, " - refusing to impute silently")
  dat[[newname]] <- z; dat
}
swine_col <- grep("농가수.*돼지", names(cor_merged), value = TRUE)[1]
poul_col  <- grep("농가수.*가금", names(cor_merged), value = TRUE)[1]
ic_fc <- ic
if (!is.na(swine_col)) ic_fc <- add_cov(ic_fc, swine_col, "swine_farm_z")
if (!is.na(poul_col))  ic_fc <- add_cov(ic_fc, poul_col,  "poultry_farm_z")
extra <- intersect(c("swine_farm_z", "poultry_farm_z"), names(ic_fc))
fc <- fit_alt(ic_fc, c(covs, extra))
# (d) add a binary COVID-era indicator (2020-2021 vs 2022-2024)
ic_cv <- ic %>% mutate(covid_era = as.numeric(year %in% c(2020, 2021)))
fd <- fit_alt(ic_cv, c(covs, "covid_era"))
tS6 <- do.call(rbind, lapply(c(covs, "covid_era", "swine_farm_z", "poultry_farm_z"), function(cv) data.frame(
  covariate = cv, principal = irr_txt(M6, cv),
  excl_2020_2021 = if (!is.null(fa)) irr_txt(fa, cv) else "-",
  covid_indicator = if (!is.null(fd)) irr_txt(fd, cv) else "-",
  drop_inpatient = if (!is.null(fb)) irr_txt(fb, cv) else "-",
  swine_poultry = if (!is.null(fc)) irr_txt(fc, cv) else "-")))
tS6 <- rbind(tS6, data.frame(covariate = "N_district_years", principal = nrow(ic),
  excl_2020_2021 = sum(!ic$year %in% c(2020, 2021)), covid_indicator = nrow(ic), drop_inpatient = nrow(ic), swine_poultry = nrow(ic)))
write.csv(tS6, file.path(OUT_DIR, "alternative_specifications.csv"), row.names = FALSE)
# (e) oyster production: the source has only 2020, 2022 and 2023, and the 2020 file omits the major producing
#     districts Tongyeong, Goseong and Yeosu. Replace the annual score by a time-invariant producer indicator, or drop it.
KEY_OY <- covs                                  # all 27 covariates
rows_oy <- function(f, nm) { k <- intersect(KEY_OY, rownames(f$summary.fixed)); fe <- f$summary.fixed[k, ]
  data.frame(spec = nm, covariate = k, IRR = round(exp(fe$mean), 6), lo = round(exp(fe$`0.025quant`), 6), hi = round(exp(fe$`0.975quant`), 6),
             credible = as.integer(fe$`0.025quant` > 0 | fe$`0.975quant` < 0), DIC = round(f$dic$dic, 2), WAIC = round(f$waic$waic, 2)) }
ic_oy <- ic; ever <- tapply(ic$oyster_z, ic$region, max); ic_oy$oyster_z <- as.numeric(ever[as.character(ic$region)])
fe1 <- fit_alt(ic_oy, covs); fe2_ <- fit_alt(ic, setdiff(covs, "oyster_z"))
write.csv(rbind(rows_oy(M6, "principal"), rows_oy(fe1, "oyster_ever_producer"), rows_oy(fe2_, "oyster_dropped")),
          file.path(OUT_DIR, "sens_oyster_coverage.csv"), row.names = FALSE)
print(tS6, row.names = FALSE)
# (f) extreme-value screen (post hoc). Among rate/percentage/cost covariates, values beyond 5 interquartile ranges
#     from the quartiles are flagged. Several of them (sex ratio of 49 and 51 males per 100 females; basic-livelihood
#     recipients of 35-38%) look like aggregation errors of the same kind as rule R3; all but one of these occur in
#     cities with non-autonomous wards. Sewer-pipe repair sites enters the model as a raw count whose maximum exceeds 500 times the third
#     quartile. Three refits: (f1) flagged values set to NA, filled within the study years by rule R5, covariate
#     re-transformed, unfillable district-years dropped; (f2) sewer-pipe repair sites as log(1 + x); (f3) both.
mk_extreme <- function(recode, sewer_log) {
  d <- ic; fl <- list()
  for (i in seq_len(nrow(TV))) {
    code <- TV$code[i]; zn <- paste0(TV$eng[i], "_z"); if (!zn %in% covs || !code %in% names(d)) next
    x <- as.numeric(d[[code]]); form <- if (TV$eng[i] == "sex_ratio") "raw" else TV$form[i]; touched <- FALSE
    if (recode && code %in% RATE_VARS) {
      q <- quantile(x, c(.25, .75), na.rm = TRUE); k <- 5 * diff(q); bad <- !is.na(x) & (x < q[1] - k | x > q[2] + k)
      if (any(bad)) {
        fl[[code]] <- data.frame(covariate = TV$eng[i], region = d$region[bad], year = d$year[bad], value = x[bad])
        x[bad] <- NA
        for (g in unique(d$region[bad])) { j <- which(d$region == g); x[j] <- nearest_fill(d$year[j], x[j])$x }
        touched <- TRUE
      }
    }
    if (sewer_log && TV$eng[i] == "sewer_repair") { form <- "log1p"; touched <- TRUE }
    if (touched) { nv <- sum(!is.na(x)); hz <- sum(!is.na(x) & x == 0) / nv * 100 > 20
      val <- apply_form(x, form, hz); d[[zn]] <- (val - mean(val, na.rm = TRUE)) / sd(val, na.rm = TRUE) }
  }
  list(data = d[complete.cases(d[, covs]), ], flags = do.call(rbind, fl))
}
ext <- list(extreme_values_recoded = mk_extreme(TRUE, FALSE), sewer_repair_log = mk_extreme(FALSE, TRUE), both = mk_extreme(TRUE, TRUE))
write.csv(ext$extreme_values_recoded$flags, file.path(OUT_DIR, "extreme_value_flags.csv"), row.names = FALSE, fileEncoding = "UTF-8")
sx <- do.call(rbind, lapply(names(ext), function(nm) { dd <- ext[[nm]]$data; f <- fit_alt(dd, covs); fe_ <- f$summary.fixed[covs, ]
  cat(sprintf("  (f) %-24s N = %d | credible: %s\n", nm, nrow(dd), paste(sub("_z$", "", covs[fe_$`0.025quant` > 0 | fe_$`0.975quant` < 0]), collapse = ",")))
  data.frame(spec = nm, covariate = covs, IRR = round(exp(fe_$mean), 6), lo = round(exp(fe_$`0.025quant`), 6), hi = round(exp(fe_$`0.975quant`), 6),
             credible = as.integer(fe_$`0.025quant` > 0 | fe_$`0.975quant` < 0), N = nrow(dd), districts = length(unique(dd$region)),
             DIC = round(f$dic$dic, 2), WAIC = round(f$waic$waic, 2)) }))
write.csv(sx, file.path(OUT_DIR, "sens_extreme_values.csv"), row.names = FALSE)
cat(sprintf("  (a) exclude 2020-2021:  water_supply IRR %s | dairy %s\n",
            if (!is.null(fa)) irr_txt(fa, "water_supply_z") else "-",
            if (!is.null(fa)) irr_txt(fa, "dairy_farm_z") else "-"))
cat(sprintf("  (b) drop inpatient cost: water_supply IRR %s | dairy %s\n",
            if (!is.null(fb)) irr_txt(fb, "water_supply_z") else "-",
            if (!is.null(fb)) irr_txt(fb, "dairy_farm_z") else "-"))
cat(sprintf("  (c) + swine/poultry:     dairy IRR %s | swine %s | poultry %s\n",
            if (!is.null(fc)) irr_txt(fc, "dairy_farm_z") else "-",
            if (!is.null(fc)) irr_txt(fc, "swine_farm_z") else "-",
            if (!is.null(fc)) irr_txt(fc, "poultry_farm_z") else "-"))

# ---------------------------------------------------------------------------
# [12] Predictive diagnostics (CPO / PIT)
# ---------------------------------------------------------------------------
cat("\n## [12] Predictive diagnostics\n")
cpo_fail <- sum(M6$cpo$failure > 0, na.rm = TRUE)
cat(sprintf("  CPO failures = %d/%d | mean PIT = %.3f (well-calibrated ~ 0.5)\n",
            cpo_fail, length(M6$cpo$cpo), mean(M6$cpo$pit, na.rm = TRUE)))

cat(sprintf("  CPO < 0.001: %d\n", sum(M6$cpo$cpo < 0.001, na.rm = TRUE)))
reS <- M6$summary.random$idarea[1:nrow(shp_main), ]
tS5 <- data.frame(region = shp_main$region, mean = reS$mean, lo = reS$`0.025quant`, hi = reS$`0.975quant`) %>%
  arrange(desc(mean))
write.csv(tS5, file.path(OUT_DIR, "district_effects.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(M6_fitted = M6$summary.fitted.values$mean[seq_len(nrow(ic))], ic = ic,
             re = reS, shp_region = shp_main$region, cpo = M6$cpo), file.path(OUT_DIR, "m6_outputs.rds"))

# ---------------------------------------------------------------------------
# [13] Covariate estimates under nested / alternative random-effect structures
# ---------------------------------------------------------------------------
cat("\n## [13] Covariate estimates by random-effect structure\n")
BYM_T <- "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)"
SPEC13 <- list(M3 = paste(base_f, BYM_T),
               M4 = paste(base_f, BYM_T, "+ f(idarea_time, model='iid', hyper=pc_prec)"),
               M5 = paste(base_f, BYM_T, "+ f(idtime, model='rw1', hyper=pc_prec)"),
               M6 = paste(base_f, RE_FULL),
               M7_yearFE = paste(sub("cases ~", "cases ~ factor(year) +", base_f), BYM_T, "+ f(idarea_time, model='iid', hyper=pc_prec)"))
nest <- list()
for (nm in names(SPEC13)) {
  m <- if (nm == "M6") M6 else fitm(SPEC13[[nm]]); if (is.null(m)) next
  fe13 <- m$summary.fixed; fe13 <- fe13[grepl("_z$", rownames(fe13)), ]
  nest[[nm]] <- data.frame(model = nm, covariate = rownames(fe13), IRR = exp(fe13$mean), lo = exp(fe13$`0.025quant`), hi = exp(fe13$`0.975quant`),
    credible = as.integer(fe13$`0.025quant` > 0 | fe13$`0.975quant` < 0), DIC = m$dic$dic, WAIC = m$waic$waic, pD = m$dic$p.eff)
  cat(sprintf("  %-10s DIC %.2f WAIC %.2f pD %.1f | credible: %s\n", nm, m$dic$dic, m$waic$waic, m$dic$p.eff,
              paste(sub("_z$", "", rownames(fe13)[fe13$`0.025quant` > 0 | fe13$`0.975quant` < 0]), collapse = ",")))
}
write.csv(do.call(rbind, nest), file.path(OUT_DIR, "nested_model_covariates.csv"), row.names = FALSE)
# M6 without covariates: how much of the spatial structure do the covariates take up?
M6_null <- fitm(paste("cases ~ 1 + offset(log(population + 1))", RE_FULL))
sp_row <- function(m, nm) { hp <- m$summary.hyperpar; g <- function(k) hp[grep(k, rownames(hp), fixed = TRUE), ]
  sp <- g("spatial component"); ii <- g("iid component"); rw <- g("idtime"); it <- g("idarea_time"); rr_ <- m$summary.random$idarea
  data.frame(model = nm, sd_spatial = 1 / sqrt(sp$`0.5quant`), sd_spatial_lo = 1 / sqrt(sp$`0.975quant`), sd_spatial_hi = 1 / sqrt(sp$`0.025quant`),
    sd_unstructured = 1 / sqrt(ii$`0.5quant`), sd_rw1 = 1 / sqrt(rw$`0.5quant`), sd_interaction = 1 / sqrt(it$`0.5quant`),
    sd_district_effect = sd(rr_$mean[ana]), n_high = sum(rr_$`0.025quant`[ana] > 0), n_low = sum(rr_$`0.975quant`[ana] < 0),
    DIC = m$dic$dic, WAIC = m$waic$waic, pD = m$dic$p.eff) }
write.csv(rbind(sp_row(M6_null, "M6_without_covariates"), sp_row(M6, "M6")), file.path(OUT_DIR, "spatial_structure_with_without_covariates.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# [14] Numerical stability (set RUN_STABILITY=true)
# ---------------------------------------------------------------------------
RUN_STABILITY <- identical(tolower(Sys.getenv("RUN_STABILITY", "false")), "true")
if (RUN_STABILITY) {
  cat("\n## [14] Numerical stability\n")
  KEY <- covs                                   # all 27 covariates are monitored
  fit_v <- function(ctrl = list(), hyp_b = pc_bym, hyp_p = pc_prec, bym2 = FALSE) {
    sp <- if (bym2) "f(idarea, model='bym2', graph=g_main, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec',param=c(0.5,0.01)), phi=list(prior='pc', param=c(0.5,0.5))))"
          else "f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=hyp_b)"
    fs <- paste(base_f, "+", sp, "+ f(idtime, model='rw1', hyper=hyp_p) + f(idarea_time, model='iid', hyper=hyp_p)")
    tryCatch(inla(as.formula(fs), family = "nbinomial", data = ic, control.compute = list(dic = TRUE, waic = TRUE),
                  control.predictor = list(link = 1), control.inla = ctrl),
             error = function(e) { message("FAIL ", conditionMessage(e)); NULL })
  }
  pcw <- function(u) list(hb = list(prec.unstruct = list(prior = "pc.prec", param = c(u, 0.01)), prec.spatial = list(prior = "pc.prec", param = c(u, 0.01))),
                          hp = list(prec = list(prior = "pc.prec", param = c(u, 0.01))))
  variants <- list(default = list(), vb_off = list(control.vb = list(enable = FALSE)), laplace = list(strategy = "laplace"),
                   grid_int = list(int.strategy = "grid"), eb_int = list(int.strategy = "eb"),
                   prior_wide = pcw(1), prior_tight = pcw(0.25), bym2 = list(bym2 = TRUE))
  rows <- list(); hyp <- list()
  for (nm in names(variants)) {
    v <- variants[[nm]]
    f <- if (!is.null(v$hb)) fit_v(hyp_b = v$hb, hyp_p = v$hp) else if (isTRUE(v$bym2)) fit_v(bym2 = TRUE) else fit_v(ctrl = v)
    if (is.null(f)) { cat(sprintf("  %-11s FAILED\n", nm)); next }
    fe <- f$summary.fixed[KEY, ]
    rows[[nm]] <- data.frame(variant = nm, covariate = KEY, IRR = exp(fe$mean), lo = exp(fe$`0.025quant`), hi = exp(fe$`0.975quant`),
      credible = as.integer(fe$`0.025quant` > 0 | fe$`0.975quant` < 0), DIC = f$dic$dic, WAIC = f$waic$waic)
    hp <- f$summary.hyperpar
    hyp[[nm]] <- data.frame(variant = nm, hyper = rownames(hp), mean = hp$mean, lo = hp$`0.025quant`, hi = hp$`0.975quant`)
    cat(sprintf("  %-11s DIC %.2f WAIC %.2f | credible: %s\n", nm, f$dic$dic, f$waic$waic,
                paste(KEY[fe$`0.025quant` > 0 | fe$`0.975quant` < 0], collapse = ",")))
  }
  write.csv(do.call(rbind, rows), file.path(OUT_DIR, "stability_inla_variants.csv"), row.names = FALSE)
  write.csv(do.call(rbind, hyp), file.path(OUT_DIR, "stability_hyperparameters.csv"), row.names = FALSE)

  # glmmTMB cross-check: NB2 with district and year random intercepts, no spatial structure
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    icg <- ic; icg$fregion <- factor(icg$region); icg$fyear <- factor(icg$year)
    gl <- glmmTMB::glmmTMB(as.formula(paste("cases ~", cov_str, "+ offset(log(population + 1)) + (1 | fregion) + (1 | fyear)")),
                           family = glmmTMB::nbinom2, data = icg)
    cf <- summary(gl)$coefficients$cond
    g <- data.frame(covariate = rownames(cf), IRR = exp(cf[, 1]), lo = exp(cf[, 1] - 1.96 * cf[, 2]), hi = exp(cf[, 1] + 1.96 * cf[, 2]))
    g <- g[g$covariate %in% KEY, ]; g$credible <- as.integer(g$lo > 1 | g$hi < 1)   # "credible" = 95% CI excludes 1
    write.csv(g, file.path(OUT_DIR, "stability_glmmTMB.csv"), row.names = FALSE)
  } else cat("  glmmTMB not installed - cross-check skipped\n")

  # Repeated fits with identical settings: run-to-run numerical variation of each nested model.
  N_REP <- as.integer(Sys.getenv("HAV_N_REPEATS", "10"))
  reps <- list()
  for (nm in names(M)) for (k in seq_len(if (nm == "M6") 2 * N_REP else N_REP)) {
    m <- fitm(M[[nm]]); if (is.null(m)) next
    fe <- m$summary.fixed; fe <- fe[grepl("_z$", rownames(fe)), ]
    cr <- sort(sub("_z$", "", rownames(fe)[fe$`0.025quant` > 0 | fe$`0.975quant` < 0]))
    g4 <- function(v) if (v %in% rownames(fe)) exp(fe[v, "mean"]) else NA
    reps[[paste(nm, k)]] <- data.frame(model = nm, run = k, DIC = m$dic$dic, WAIC = m$waic$waic, pD = m$dic$p.eff,
      n_credible = length(cr), credible_set = paste(cr, collapse = ";"),
      gw_household = g4("gw_household_z"), dairy_farm = g4("dairy_farm_z"), forest = g4("forest_z"), fiscal_indep = g4("fiscal_indep_z"))
  }
  reps <- do.call(rbind, reps); write.csv(reps, file.path(OUT_DIR, "repeated_fits.csv"), row.names = FALSE)
  print(aggregate(cbind(DIC, WAIC, pD) ~ model, reps, function(x) round(range(x), 2)))
}

cat("\n===== DONE =====\n")
cat(sprintf("N = %d | M6 DIC = %.2f | WAIC = %.2f | residual Moran's I = %+.4f (p = %.3g)\n",
            nrow(ic), M6$dic$dic, M6$waic$waic,
            moran_post$estimate[[1]], moran_post$p.value))
cat(sprintf("credible covariates = %d | high/low-risk = %d/%d\n",
            sum(table2$credible), n_high, n_low))
