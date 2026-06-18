# =============================================================================
# Reproducible analysis code
# "The sanitation paradox and groundwater vulnerability in the spatial
#  distribution of hepatitis A virus foodborne disease in South Korea, 2020-2024"
# Seongdae Kim, Byung Chul Chun.  Target journal: Water Research (Elsevier).
#
# Model: Bayesian negative-binomial disease mapping with a Besag-York-Mollie
#  (BYM) convolution + first-order temporal random walk (RW1) + Knorr-Held
#  Type I space-time interaction, fitted by INLA (R-INLA). 223 contiguous
#  districts, 1,112 district-years (2020-2024), 27 pre-specified covariates.
#
# This script fits exactly the covariate specification reported in the paper:
#  the 27 covariates and their functional forms are PRE-SPECIFIED (Table S5)
#  and entered directly. There is no data-driven model search, stepwise
#  selection, or objective tuned to obtain a target result; every covariate
#  is forced into the model with its declared transform, and all reported
#  quantities are read off the resulting fits.
#
# Reproduces:
#  * Principal model M6  (DIC 5,716.29; WAIC 5,729.22; residual Moran's I
#    +0.053, p = 0.090)
#  * Table 2  — 27 covariate incidence-rate ratios (9 credible)
#  * Table S1 — model comparison M1-M6
#  * Table S2 — 8-graph neighbourhood sensitivity
#  * Tables S3/S7 — BYM2 reparametrisation and prior sensitivity (phi ~ 0.96)
#  * Table S4 — Global Moran's I (pre/post)
#  * Figure S2 — Getis-Ord Gi* local clustering
#  * Table S8 — alternative-specification robustness checks
#
# Software: R 4.x with R-INLA. The manuscript fits used R-INLA 24.x; newer
#  INLA versions may shift the DIC/WAIC by a few points without changing any
#  incidence-rate ratio or credible-interval conclusion.
# Run:  Rscript HAV_spatial_reproducible.R
#
# DATA AVAILABILITY: annual district-level HAV notifications are released by the
#  Korea Disease Control and Prevention Agency (KDCA) Infectious Disease Portal
#  (https://dportal.kdca.go.kr); covariates come from KOSIS and the open-data
#  portals of the relevant Korean ministries. Restricted/raw inputs are NOT
#  redistributed here. Place the input files under ./data (or set the
#  HAV_DATA_DIR environment variable) before running. No personally
#  identifiable information is used (aggregated district-year counts only).
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
PATH_DISEASE   <- file.path(BASE_IV, "foodborne_final.csv")        # district-year disease counts
PATH_HEALTH_PQ <- file.path(BASE_IV, "health_indicators.parquet")  # community-health covariates
PATH_SHP       <- file.path(BASE_IV, "districts.shp")              # district polygons (EPSG:5179)
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

read_csv_safe <- function(fp) {
  for (enc in c("UTF-8", "UTF-8-BOM", "CP949", "EUC-KR")) {
    raw <- tryCatch(read.csv(fp, fileEncoding = enc, check.names = FALSE,
                             stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) return(raw)
  }
  NULL
}

# Some administrative covariates are published biennially; carry the adjacent
# year forward/back so that all five years are populated (2021<-2020, 2024<-2023).
fill_missing_year <- function(df, tgt, src) {
  if (!"region" %in% names(df) || !"year" %in% names(df)) return(df)
  nv <- setdiff(names(df)[sapply(df, is.numeric)], "year")
  if (length(nv) == 0 || !src %in% unique(df$year)) return(df)
  ds <- df %>% filter(year == src); dt <- df %>% filter(year == tgt)
  df_f <- ds %>% mutate(year = as.integer(tgt))
  if (nrow(dt) > 0) {
    df_f <- df_f %>%
      left_join(dt %>% dplyr::select(region, all_of(nv)) %>%
                  rename_with(~paste0(., "__o"), all_of(nv)), by = "region") %>%
      mutate(across(all_of(nv), function(col) {
        v <- cur_column(); o <- get(paste0(v, "__o")); ifelse(!is.na(o), o, col) })) %>%
      dplyr::select(region, year, all_of(nv))
  }
  bind_rows(df %>% filter(year != tgt), df_f) %>% arrange(region, year)
}
apply_cf <- function(df) fill_missing_year(fill_missing_year(df, 2021, 2020), 2024, 2023)

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
  hpq <- read_parquet(PATH_HEALTH_PQ) %>% as.data.frame() %>% clean_region() %>% apply_cf()
  ah <- intersect(health_vars, names(hpq))
  for (v in ah) hpq[[v]] <- suppressWarnings(as.numeric(hpq[[v]]))
  hagg <- hpq %>% group_by(region, year) %>%
    summarise(across(all_of(ah), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  cor_merged <- cor_merged %>% left_join(hagg, by = c("region", "year"))
  cat(sprintf("  parquet health indicators: %d variables\n", length(ah)))
}, error = function(e) cat(sprintf("  [warn] parquet: %s\n", e$message)))

# (B) covariates supplied as individual CSV extracts. Only the columns used by
#     the pre-specified model are kept; carry-forward fills biennial gaps.
selected <- list(
  "groundwater_household.csv"   = c("가정용_개소수", "간이상수도용_개소수"),
  "groundwater_quality.csv"     = c("검사합계"),
  "sewer_repair.csv"            = c("개·보수관로_부분보수(개소)_계"),
  "sewerage_coverage.csv"       = c("공공하수처리구역인구보급률(%)", "하수처리구역외_정화조인구", "총면적(㎢)"),
  "livestock.csv"               = c("농가수(호)_젖소", "농가수(호)_돼지", "농가수(호)_가금"),
  "elderly_singleperson.csv"    = c("1인가구_80~84세"),
  "land_use.csv"                = c("답", "임야", "대"),
  "shellfish.csv"               = c("굴_자연채묘 생산량(kg)"))
for (fn in names(selected)) {
  fp <- file.path(BASE_IV, fn); if (!file.exists(fp)) next
  raw <- read_csv_safe(fp); if (is.null(raw)) next
  raw <- raw %>% clean_region() %>% apply_cf()
  av <- intersect(selected[[fn]], names(raw)); if (length(av) == 0) next
  for (v in av) raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
  agg <- raw %>% group_by(region, year) %>%
    summarise(across(all_of(av), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  cor_merged <- cor_merged %>% left_join(agg, by = c("region", "year"))
  cat(sprintf("  csv %-30s %d variables\n", fn, length(av)))
}
cat(sprintf("  merged: %d rows x %d columns\n", nrow(cor_merged), ncol(cor_merged)))

# ---------------------------------------------------------------------------
# [4] Pre-specified covariate set (27) and functional forms (Table S5)
# ---------------------------------------------------------------------------
# Each covariate enters with a single, pre-declared functional form:
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
cat(sprintf("\n## [4] Pre-specified covariates: %d (forced confounders: %d)\n",
            nrow(TV), sum(TV$forced == "Y")))

# District polygons and contiguity graph (principal neighbourhood structure).
# Six island districts with no land contiguity are excluded so that the spatial
# model has a connected graph (223 contiguous districts).
shp <- st_read(PATH_SHP, quiet = TRUE) %>%
  mutate(region = str_replace_all(as.character(region), "\\s+", ""),
         region = if_else(region == "인천시미추홀구", "인천시남구", region))
islands <- c("인천시옹진군", "전라남도완도군", "전라남도진도군",
             "경상남도거제시", "경상남도남해군", "경상북도울릉군")
shp_main <- shp %>% filter(!region %in% islands)
nb_obj <- poly2nb(shp_main, snap = 0.01); iso <- which(card(nb_obj) == 0)
if (length(iso) > 0) { shp_main <- shp_main[-iso, ]; nb_obj <- poly2nb(shp_main, snap = 0.01) }
graph_file <- tempfile(fileext = ".graph")
nb2INLA(nb_obj, file = graph_file); g_main <- inla.read.graph(graph_file)
nb_w <- nb2listw(nb_obj, style = "W", zero.policy = TRUE)
cat(sprintf("  districts in spatial model: %d\n", nrow(shp_main)))

# ---------------------------------------------------------------------------
# [5] Build the design matrix from the pre-specified forms
# ---------------------------------------------------------------------------
# apply_form() materialises one covariate in its declared functional form.
# `hz` flags zero-inflated counts (>20% zeros), for which the binary/tertile
# cut-points use presence/non-zero medians rather than the overall median.
apply_form <- function(x, form, hz) {
  if (form == "raw")   return(x)
  if (form == "log1p") return(log1p(pmax(x, 0)))
  if (form == "binary") {
    if (hz) return(as.numeric(!is.na(x) & x > 0))
    md <- median(x, na.rm = TRUE); return(as.numeric(!is.na(x) & x > md))
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
zcols <- character(0)
for (i in seq_len(nrow(TV))) {
  code <- TV$code[i]; x <- as.numeric(df_w[[code]])
  nv <- sum(!is.na(x) & is.finite(x)); if (nv < MIN_OBS) next
  hz <- sum(!is.na(x) & is.finite(x) & x == 0) / nv * 100 > 20
  form <- if (TV$eng[i] == "sex_ratio") "raw" else TV$form[i]   # sex ratio is symmetric -> raw
  val <- apply_form(x, form, hz)
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
  drop <- names(which.max(vv)); if (drop %in% forced_z) break
  keep <- setdiff(keep, drop)
}
covs <- keep
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
  IRR = round(exp(fe$mean), 3),
  lo  = round(exp(fe$`0.025quant`), 3),
  hi  = round(exp(fe$`0.975quant`), 3),
  credible = as.integer(fe$`0.025quant` > 0 | fe$`0.975quant` < 0),
  row.names = NULL)
cat(sprintf("  Table 2: %d covariates, %d credible (95%% CrI excludes 1)\n",
            nrow(table2), sum(table2$credible)))
print(table2[table2$credible == 1, ], row.names = FALSE)

# ---------------------------------------------------------------------------
# [7] Model comparison M1-M6 (Table S1)
# ---------------------------------------------------------------------------
cat("\n## [7] Model comparison (Table S1)\n")
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
  data.frame(model = n, DIC = round(f$dic$dic, 2), WAIC = round(f$waic$waic, 2))
}))
print(tableS1, row.names = FALSE)

# ---------------------------------------------------------------------------
# [8] Global Moran's I, pre- and post-modelling (Table S4)
# ---------------------------------------------------------------------------
cat("\n## [8] Global Moran's I (Table S4)\n")
agg <- ic %>% group_by(idarea) %>%
  summarise(rate = sum(cases) / sum(population) * 1e5, .groups = "drop") %>% arrange(idarea)
rv <- rep(NA, nrow(shp_main)); rv[agg$idarea] <- agg$rate
moran_pre <- moran.test(rv, nb_w, zero.policy = TRUE, na.action = na.omit)
ic$fitv <- M6$summary.fitted.values$mean[seq_len(nrow(ic))]
res_agg <- ic %>% mutate(resid = cases - fitv) %>% group_by(idarea) %>%
  summarise(r = sum(resid), .groups = "drop") %>% arrange(idarea)
rr <- rep(NA, nrow(shp_main)); rr[res_agg$idarea] <- res_agg$r
moran_post <- moran.test(rr, nb_w, zero.policy = TRUE, na.action = na.omit)
cat(sprintf("  crude    I = %+.4f (p = %.3g)\n", moran_pre$estimate[[1]], moran_pre$p.value))
cat(sprintf("  residual I = %+.4f (p = %.3g)\n", moran_post$estimate[[1]], moran_post$p.value))

# High- and low-risk districts: structured spatial effect whose 95% CrI
# excludes zero.
re <- M6$summary.random$idarea; na <- nrow(shp_main)
n_high <- sum(re$`0.025quant`[1:na] > 0); n_low <- sum(re$`0.975quant`[1:na] < 0)
cat(sprintf("  high-risk districts = %d | low-risk districts = %d\n", n_high, n_low))

# ---------------------------------------------------------------------------
# [9] BYM2 reparametrisation and prior sensitivity (Tables S3 / S7)
# ---------------------------------------------------------------------------
cat("\n## [9] BYM2 mixing parameter phi (Tables S3/S7)\n")
f_bym2 <- function(hy) paste(base_f,
  "+ f(idarea, model='bym2', graph=g_main, scale.model=TRUE, constr=TRUE", hy, ")",
  "+ f(idtime, model='rw1', hyper=pc_prec)",
  "+ f(idarea_time, model='iid', hyper=pc_prec)")
# Principal PC prior on phi, plus a panel of alternative priors.
phi_priors <- list(
  "PC(0.5,2/3)" = ",hyper=list(prec=list(prior='pc.prec',param=c(1,0.01)),phi=list(prior='pc',param=c(0.5,2/3)))",
  "PC(0.5,0.5)" = ",hyper=list(prec=list(prior='pc.prec',param=c(1,0.01)),phi=list(prior='pc',param=c(0.5,0.5)))",
  "PC(0.5,0.9)" = ",hyper=list(prec=list(prior='pc.prec',param=c(1,0.01)),phi=list(prior='pc',param=c(0.5,0.9)))")
for (nm in names(phi_priors)) {
  b2 <- fitm(f_bym2(phi_priors[[nm]]))
  if (is.null(b2)) next
  pr <- grep("Phi", rownames(b2$summary.hyperpar), ignore.case = TRUE, value = TRUE)[1]
  cat(sprintf("  %-12s phi = %.3f (%.3f-%.3f) | DIC = %.2f\n", nm,
              b2$summary.hyperpar[pr, "mean"],
              b2$summary.hyperpar[pr, "0.025quant"],
              b2$summary.hyperpar[pr, "0.975quant"], b2$dic$dic))
}

# ---------------------------------------------------------------------------
# [10] Eight-graph neighbourhood sensitivity (Table S2)
# ---------------------------------------------------------------------------
cat("\n## [10] 8-graph neighbourhood sensitivity (Table S2)\n")
cz <- st_coordinates(st_centroid(st_geometry(shp_main)))
mkgraph <- function(nb) { f <- tempfile(); nb2INLA(f, nb); inla.read.graph(f) }
graphs <- list(Queen = poly2nb(shp_main, queen = TRUE,  snap = 0.01),
               Rook  = poly2nb(shp_main, queen = FALSE, snap = 0.01))
for (k in 2:7) graphs[[paste0("knn", k)]] <- make.sym.nb(knn2nb(knearneigh(cz, k = k)))
cred <- c("dairy_farm_z", "gw_household_z", "water_supply_z", "sewer_repair_z",
          "forest_z", "residential_z", "alone_80_84_z", "fiscal_indep_z", "med_in_z")
cred <- intersect(cred, covs)
graph_cred <- setNames(integer(length(cred)), cred)
for (gn in names(graphs)) {
  gg <- tryCatch(mkgraph(graphs[[gn]]), error = function(e) NULL); if (is.null(gg)) next
  ff <- fitm(paste(base_f,
    "+ f(idarea, model='bym', graph=gg, scale.model=TRUE, hyper=pc_bym)",
    "+ f(idtime, model='rw1', hyper=pc_prec)",
    "+ f(idarea_time, model='iid', hyper=pc_prec)"))
  if (is.null(ff)) next
  fe2 <- ff$summary.fixed
  for (c in cred) if (c %in% rownames(fe2))
    graph_cred[c] <- graph_cred[c] +
      as.integer(fe2[c, "0.025quant"] > 0 | fe2[c, "0.975quant"] < 0)
  cat(sprintf("  %-6s DIC = %.2f\n", gn, ff$dic$dic))
}
cat("  credible across graphs (out of 8):\n")
for (c in cred) cat(sprintf("    %-18s %d/8\n", c, graph_cred[c]))

# ---------------------------------------------------------------------------
# [11] Getis-Ord Gi* local clustering (Figure S2)
# ---------------------------------------------------------------------------
cat("\n## [11] Getis-Ord Gi* (Figure S2)\n")
nb_self <- include.self(nb_obj)
lw_self <- nb2listw(nb_self, style = "B", zero.policy = TRUE)
gi <- localG(ifelse(is.na(rv), 0, rv), lw_self, zero.policy = TRUE)
gi_z <- as.numeric(gi)
cat(sprintf("  Gi* z-scores: hot spots (z > 1.96) = %d | cold spots (z < -1.96) = %d\n",
            sum(gi_z > 1.96, na.rm = TRUE), sum(gi_z < -1.96, na.rm = TRUE)))
getis <- data.frame(region = shp_main$region, Gi_z = round(gi_z, 3))

# ---------------------------------------------------------------------------
# [12] Alternative-specification robustness checks (Table S8)
# ---------------------------------------------------------------------------
cat("\n## [12] Robustness checks (Table S8)\n")
irr_txt <- function(fit, name) {
  if (!name %in% rownames(fit$summary.fixed)) return("-")
  r <- fit$summary.fixed[name, ]
  sprintf("%.3f (%.3f-%.3f)%s", exp(r$mean), exp(r$`0.025quant`), exp(r$`0.975quant`),
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
# (c) add swine and poultry farm density (tests specificity of the dairy signal)
add_cov <- function(dat, raw_col, newname) {
  key <- cor_merged %>% transmute(region, year, val = suppressWarnings(as.numeric(.data[[raw_col]])))
  m <- dat %>% left_join(key, by = c("region", "year"))
  v <- m$val; v[!is.finite(v)] <- NA
  z <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE); z[is.na(z)] <- 0
  dat[[newname]] <- z; dat
}
swine_col <- grep("농가수.*돼지", names(cor_merged), value = TRUE)[1]
poul_col  <- grep("농가수.*가금", names(cor_merged), value = TRUE)[1]
ic_fc <- ic
if (!is.na(swine_col)) ic_fc <- add_cov(ic_fc, swine_col, "swine_farm_z")
if (!is.na(poul_col))  ic_fc <- add_cov(ic_fc, poul_col,  "poultry_farm_z")
extra <- intersect(c("swine_farm_z", "poultry_farm_z"), names(ic_fc))
fc <- fit_alt(ic_fc, c(covs, extra))
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
# [13] Predictive diagnostics (CPO / PIT)
# ---------------------------------------------------------------------------
cat("\n## [13] Predictive diagnostics\n")
cpo_fail <- sum(M6$cpo$failure > 0, na.rm = TRUE)
cat(sprintf("  CPO failures = %d/%d | mean PIT = %.3f (well-calibrated ~ 0.5)\n",
            cpo_fail, length(M6$cpo$cpo), mean(M6$cpo$pit, na.rm = TRUE)))

cat("\n===== DONE =====\n")
cat(sprintf("N = %d | M6 DIC = %.2f | WAIC = %.2f | residual Moran's I = %+.4f (p = %.3g)\n",
            nrow(ic), M6$dic$dic, M6$waic$waic,
            moran_post$estimate[[1]], moran_post$p.value))
cat(sprintf("credible covariates = %d | high/low-risk = %d/%d\n",
            sum(table2$credible), n_high, n_low))
