# =============================================================================
# HAV JMIR PHS #112448 revision R1 - additional analyses (v3, 2026-10-31, after GPT audit rounds 1-2)
# v3: sign-preserving log in every log1p form (net migration can be negative; v2 clipped it in the mixed-form arms);
#   HAV_FRAME / HAV_ARMS environment variables allow the revised principal model to be refitted on alternative input frames.
# v2 changes: input frame = pipeline rerun with corrected Busanjin-gu land-use/well files (G4_thresholds/out2_base);
#   district index re-derived from the boundary file by code+geometry (the submitted frame had Incheon Nam-gu and Namdong-gu interchanged);
#   signed log for the uniform log form; SVC local slopes as INLA lincombs (beta + b_i); Moran tests two-sided with a fixed queen reference;
#   (y-mu)/sqrt(mu) named Poisson-scaled; repeated fits of A0 and A5 for numerical stability.
# Self-contained: reads the archived analysis frame (results/analysis_dataset_compiled.csv, v2.3.3),
# the district boundary file, and the KDCA district x age-group notifications and age-group populations
# used for indirect age standardisation. Every arm writes a coefficient table, fit indices and a
# residual Moran's I; arms are wrapped in tryCatch so one failure does not stop the run.
# Arms (editor/reviewer items in brackets):
#  A0  baseline M6 as submitted (sanity check: 5 credible covariates)
#  A1  offset log(population) instead of log(population + 1)                     [E-m1, AB-9]
#  A2  indirectly age-standardised expected counts as offset                      [E-M1, AB-1]
#  A3  population offset + proportion aged 20-49 years as covariate               [E-M1, AB-1]
#  A4  rescaled covariates (shares / per 1000 / per km2) + population density     [E-M2, AB-2, AR-2]
#  A5  A4 + age-standardised offset (candidate revised principal model)           [E-M1, E-M2]
#  A6  uniform form: all covariates continuous untransformed (rescaled), A5 offset [E-M3, AB-3]
#  A7  uniform form: all covariates continuous signed log(1+|x|) (rescaled), A5 offset [E-M3, AB-3]
#  A8  parsimonious models P1 (credible set) and P2 (credible set + demographic)  [E-M9, BZ-3]
#  A9  restricted spatial regression version of A5                                [E-M9, AB-8]
#  A10 zero-inflated negative binomial (type 1) for A0 and A5                     [E-M4]
#  A11 nested structures M1-M6 with residual Moran's I (A0 and A5 covariates)     [E-M5, BZ-1]
#  A12 eight neighbourhood graphs for A5 with residual Moran's I                  [E-M5]
#  A13 spatially varying coefficients (besag slopes) for forest share and wells   [E-M6, AB-4]
#  A14 ascertainment proxy: physicians per 1000 population added to A5           [E-M7, AB-5]
#  A15 A5 excluding 2020-2021                                                     [E-M8, AB-6]
# Run (from the repository root): HAV_DATA_DIR=/path/to/data Rscript revision_R1/HAV_revision_R1.R   (about 60-90 min; log in revision_R1/results/run.log)
# HAV_FRAME=<compiled frame> HAV_ARMS=A5 HAV_OUT=<dir> refits the revised principal model on an alternative input frame (run_g4_a5.sh).
# =============================================================================
suppressMessages({library(dplyr); library(sf); library(spdep); library(INLA)})
ROOT <- Sys.getenv("HAV_R1_DIR", file.path(getwd(), "revision_R1")); OUT <- Sys.getenv("HAV_OUT", file.path(ROOT, "results")); dir.create(OUT, showWarnings = FALSE)
COMPILED <- Sys.getenv("HAV_FRAME", file.path(ROOT, "data/analysis_dataset_compiled_v2.4.0.csv"))   # pipeline rerun with corrected inputs (Busanjin-gu)
SHP  <- file.path(Sys.getenv("HAV_DATA_DIR", file.path(getwd(), "data")), "final.shp")   # district boundary file (not redistributed; see DATA_PROVENANCE.md)
CACHE <- file.path(ROOT, "data")
DICT <- read.csv(file.path(ROOT, "data/covariate_dictionary_v2.4.0.csv"), check.names = FALSE, stringsAsFactors = FALSE)
T2   <- read.csv(file.path(getwd(), "results/table2_principal_IRR.csv"), stringsAsFactors = FALSE)   # submitted credible set (archived)
set.seed(20261031)
sink(file.path(OUT, "run.log"), split = TRUE)
cat("HAV revision R1 additional analyses -", format(Sys.time()), "\n")

# ---------------------------------------------------------------- 1. data
C <- read.csv(COMPILED, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
stopifnot(all(c("region","year","cases","population","idarea","idtime") %in% names(C)))
cat(sprintf("analysis frame: %d rows, %d districts, %d cases\n", nrow(C), n_distinct(C$region), sum(C$cases)))
shp <- st_read(SHP, quiet = TRUE); st_crs(shp) <- 5179   # no .prj; coordinates are Korea 2000 Unified (EPSG:5179) as in the submitted script
shp$region <- gsub("\\s+", "", as.character(shp$region)); shp$region[shp$region == "인천시미추홀구"] <- "인천시남구"
shp$area_km2 <- as.numeric(st_area(shp)) / 1e6
islands <- c("인천시옹진군", "전라남도완도군", "전라남도진도군", "경상남도거제시", "경상남도남해군", "경상북도울릉군")   # as submitted
shp_main <- shp[!shp$region %in% islands, ]; nb <- poly2nb(shp_main, snap = 0.01); iso <- which(card(nb) == 0)
if (length(iso) > 0) { shp_main <- shp_main[-iso, ]; nb <- poly2nb(shp_main, snap = 0.01) }
g <- nb2mat(nb, style = "B", zero.policy = TRUE); g_main <- inla.read.graph(g); NA_ <- nrow(shp_main)
cat(sprintf("graph: %d districts with neighbours (of %d)\n", NA_, nrow(shp)))
# idarea in the archived frame must equal the row position in shp_main
rm_chk <- data.frame(region = shp_main$region, idarea_shp = seq_len(NA_))
chk <- C %>% distinct(region, idarea) %>% left_join(rm_chk, by = "region")
stopifnot(!any(is.na(chk$idarea_shp)))
cat(sprintf("idarea check: %d of %d districts differ from the archived index (index re-derived from this boundary order)\n", sum(chk$idarea != chk$idarea_shp), nrow(chk)))
if (any(chk$idarea != chk$idarea_shp)) write.csv(chk[chk$idarea != chk$idarea_shp, ], file.path(OUT, "idarea_remap.csv"), row.names = FALSE, fileEncoding = "UTF-8")
# v2: the index is RE-DERIVED from the boundary file. Audit (01:20): rows 52 (code 23030, 24.6 km2, centroid 126.665E 37.453N)
# and 54 (code 23050, 57.8 km2, 126.727E 37.431N) are the true Nam-gu (Michuhol-gu) and Namdong-gu polygons; the archived frame
# carried them interchanged. Cases, covariates and expected counts join by district name; polygons by this corrected index.
C <- C %>% dplyr::select(-idarea) %>% left_join(rm_chk %>% rename(idarea = idarea_shp), by = "region") %>% arrange(idarea, idtime) %>% mutate(idarea_time = row_number())
stopifnot(!anyDuplicated(C[c("region","year")]), !anyDuplicated(shp_main$region), all(C$population > 0), all(C$cases >= 0), all(shp$area_km2 > 0))
C <- C %>% left_join(data.frame(region = shp$region, area_km2 = shp$area_km2), by = "region")
stopifnot(!any(is.na(C$area_km2)))
# note: the sewerage-statistics column 총면적(㎢) in the frame is ~2x the cadastral area in most districts; not used.
C$area_ratio_check <- C[["총면적(㎢)"]] / C$area_km2
write.csv(C %>% group_by(region) %>% summarise(sewer_area = mean(`총면적(㎢)`), shp_area = mean(area_km2), ratio = mean(area_ratio_check)),
          file.path(OUT, "area_check_sewer_vs_shapefile.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# age-standardised expected counts (indirect; national age-specific rates by year from all 229 districts)
ac <- read.csv(file.path(CACHE, "ST_HepatitisA_agecase.csv"), fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
ap <- read.csv(file.path(CACHE, "ST_agepop.csv"), fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
ac <- ac %>% filter(year >= 2020, year <= 2024); ap <- ap %>% filter(year >= 2020, year <= 2024)
names(ap)[names(ap) == "age10"] <- "age_group"
nat <- ac %>% group_by(year, age_group) %>% summarise(cases = sum(cases), .groups = "drop") %>%
  inner_join(ap %>% group_by(year, age_group) %>% summarise(pop = sum(pop), .groups = "drop"), by = c("year","age_group")) %>%
  mutate(rate = cases / pop)
write.csv(nat, file.path(OUT, "national_age_rates_2020_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")
E <- ap %>% inner_join(nat %>% dplyr::select(year, age_group, rate), by = c("year","age_group")) %>%
  group_by(region, year) %>% summarise(E_age = sum(pop * rate), pop_age = sum(pop),
    p2049 = sum(pop[age_group %in% c("20~29세","30~39세","40~49세")]) / sum(pop), .groups = "drop")
C <- C %>% left_join(E, by = c("region","year"))
stopifnot(!any(is.na(C$E_age)), all(C$E_age > 0))
cat(sprintf("E_age: sum %.1f vs cases %d (SIR %.3f); population check ratio %.3f-%.3f\n",
            sum(C$E_age), sum(C$cases), sum(C$cases)/sum(C$E_age), min(C$population/C$pop_age), max(C$population/C$pop_age)))
# physicians per 1000 (ascertainment proxy) from the FBD canonical panel
cp <- tryCatch(read.csv(file.path(CACHE, "physicians_count.csv"), check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"), error = function(e) NULL)
if (!is.null(cp) && "physicians_count" %in% names(cp)) {
  C <- C %>% left_join(cp %>% dplyr::select(region, year, physicians_count) , by = c("region","year"))
  C$phys_p1k <- C$physicians_count / C$population * 1000
  cat(sprintf("physicians: %d rows merged, %d NA\n", sum(!is.na(C$phys_p1k)), sum(is.na(C$phys_p1k))))
}

# ---------------------------------------------------------------- 2. covariate forms
zs <- function(v) { s <- sd(v, na.rm = TRUE); if (is.na(s) || s == 0) return(v); (v - mean(v, na.rm = TRUE)) / s }
apply_form <- function(x, form, hz) {
  if (form == "raw")   return(x)
  if (form == "log1p") return(sign(x) * log1p(abs(x)))   # v3: sign-preserving
  if (form == "binary") { if (hz) return(as.numeric(x > 0)); return(as.numeric(x > median(x, na.rm = TRUE))) }
  if (form == "T3") { if (hz) { nz <- x[x > 0]; mn <- median(nz); return(dplyr::case_when(x == 0 ~ 1, x <= mn ~ 2, TRUE ~ 3)) }
    q <- quantile(x, c(1/3, 2/3), na.rm = TRUE); return(as.numeric(cut(x, unique(c(-Inf, q, Inf)), labels = FALSE, include.lowest = TRUE))) }
  if (form == "Q4") { q <- quantile(x, c(.25,.5,.75), na.rm = TRUE); return(as.numeric(cut(x, unique(c(-Inf, q, Inf)), labels = FALSE, include.lowest = TRUE))) }
  stop(form)
}
# original 27 covariates: source column and form from the archived dictionary
ORIG <- DICT$covariate; names(ORIG) <- DICT$covariate
Z_ORIG <- paste0(DICT$covariate, "_z"); stopifnot(all(Z_ORIG %in% names(C)))
# rescaled raw values (same names as original covariates; rescaling only where the source was a count or an area)
pop <- C$population; km2 <- C$area_km2
RES <- list(
  oyster        = C[["굴_자연채묘 생산량(kg)"]],                    # kept (production, zero-inflated)
  dairy_farm    = C[["농가수(호)_젖소"]] / km2,                      # farms per km2
  gw_household  = C[["가정용_개소수"]] / pop * 1000,                 # wells per 1000 population
  gw_simple     = C[["간이상수도용_개소수"]] / pop * 1000,
  water_supply  = C[["상수도보급률"]],
  septic_pop    = C[["하수처리구역외_정화조인구"]] / pop * 100,       # % of population
  ww_facility   = C[["폐수배출업소수"]] / km2,                       # facilities per km2
  pub_sewage    = C[["공공하수처리구역인구보급률(%)"]],
  test_total    = C[["검사합계"]] / pop * 1000,                      # groundwater-quality tests per 1000 population
  sewer_repair  = C[["개·보수관로_부분보수(개소)_계"]] / pop * 1000,  # repair sites per 1000 population
  paddy         = C[["답"]] / (km2 * 1e6) * 100,                    # % of district area
  forest        = C[["임야"]] / (km2 * 1e6) * 100,
  residential   = C[["대"]] / (km2 * 1e6) * 100,
  health_practice = C[["건강생활실천율_조율"]],
  food_safety   = C[["식품안정성확보율_표준화율"]],
  elderly_alone = C[["독거노인비율"]],
  alone_80_84   = C[["1인가구_80~84세"]] / pop * 1000,               # per 1000 population
  welfare       = C[["기초생활수급자수율"]],
  fiscal_indep  = C[["재정자립도"]],
  fiscal_auto   = C[["재정자주도"]],
  sex_ratio     = C[["성비"]],
  elderly_rate  = C[["고령인구비율"]],
  urban_pop_rate = C[["도시지역인구비율"]],
  net_migration = C[["순이동인구"]] / pop * 1000,                    # per 1000 population
  med_in        = C[["관내진료비_입원"]] / pop,                       # 1000 KRW per capita
  clinic_family = C[["의원_가정의학과"]] / pop * 1e5,                 # clinics per 100 000 population
  depression    = C[["우울감경험률_표준화율"]])
stopifnot(setequal(names(RES), DICT$covariate))
RESCALED <- c("dairy_farm","gw_household","gw_simple","septic_pop","ww_facility","test_total","sewer_repair",
              "paddy","forest","residential","alone_80_84","net_migration","med_in","clinic_family")
C$pop_density <- pop / km2
# A4-type design: original forms applied to rescaled values (+ log population density, continuous)
for (v in names(RES)) {
  x <- RES[[v]]; form <- DICT$form[DICT$covariate == v]; if (v == "sex_ratio") form <- "raw"
  hz <- mean(x == 0, na.rm = TRUE) > 0.2
  C[[paste0(v, "_r")]] <- zs(apply_form(x, form, hz))
  C[[paste0(v, "_u")]] <- zs(x)                      # uniform: continuous untransformed
  C[[paste0(v, "_l")]] <- zs(sign(x) * log1p(abs(x)))  # uniform: continuous signed log(1 + |x|) (net migration can be negative)
}
C$popdens_z  <- zs(log(C$pop_density)); C$p2049_z <- zs(C$p2049 * 100)
if ("phys_p1k" %in% names(C)) C$phys_z <- zs(C$phys_p1k)
Z_R <- paste0(names(RES), "_r"); Z_U <- paste0(names(RES), "_u"); Z_L <- paste0(names(RES), "_l")
desc <- do.call(rbind, lapply(names(RES), function(v) { x <- RES[[v]]; data.frame(covariate = v, rescaled = v %in% RESCALED,
  form = DICT$form[DICT$covariate == v], mean = mean(x), sd = sd(x), min = min(x), q1 = quantile(x,.25), median = median(x), q3 = quantile(x,.75), max = max(x), share_zero = mean(x == 0)) }))
desc <- rbind(desc, data.frame(covariate = "pop_density", rescaled = TRUE, form = "log", mean = mean(C$pop_density), sd = sd(C$pop_density), min = min(C$pop_density),
  q1 = quantile(C$pop_density,.25), median = median(C$pop_density), q3 = quantile(C$pop_density,.75), max = max(C$pop_density), share_zero = 0))
write.csv(desc, file.path(OUT, "rescaled_covariate_dictionary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(C, file.path(OUT, "analysis_frame_R1.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------- 3. model infrastructure (as submitted)
pc_bym <- list(prec.unstruct = list(prior = "pc.prec", param = c(0.5, 0.01)), prec.spatial = list(prior = "pc.prec", param = c(0.5, 0.01)))
pc_prec <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
RE_M6 <- "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym) + f(idtime, model='rw1', hyper=pc_prec) + f(idarea_time, model='iid', hyper=pc_prec)"
nb_w <- nb2listw(nb, style = "W", zero.policy = TRUE)
fit <- function(covs, offset = "log(population + 1)", re = RE_M6, dat = C, family = "nbinomial", lincomb = NULL) {
  fs <- paste("cases ~", if (length(covs)) paste(covs, collapse = " + ") else "1", "+ offset(", offset, ")", re)
  tryCatch(inla(as.formula(fs), family = family, data = dat, lincomb = lincomb, control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor = list(link = 1)), error = function(e) { message("fit error: ", e$message); NULL })
}
coeftab <- function(m, arm) { if (is.null(m)) return(NULL); fe <- m$summary.fixed; fe <- fe[rownames(fe) != "(Intercept)", , drop = FALSE]; if (nrow(fe) == 0) return(NULL)
  data.frame(arm = arm, covariate = sub("_(z|r|u|l)$", "", rownames(fe)), term = rownames(fe), IRR = exp(fe$mean), lo = exp(fe$`0.025quant`), hi = exp(fe$`0.975quant`),
             credible = as.integer(fe$`0.025quant` > 0 | fe$`0.975quant` < 0), row.names = NULL) }
nb_w_ref <- nb_w   # fixed queen weights for like-for-like diagnostics across graph arms
resid_moran <- function(m, dat = C, lw = nb_w) { if (is.null(m)) return(rep(NA, 6))
  fv <- m$summary.fitted.values$mean[seq_len(nrow(dat))]
  d <- dat %>% mutate(res = cases - fv, pr = (cases - fv) / sqrt(fv)) %>% group_by(idarea) %>% summarise(r = sum(res), p = mean(pr), .groups = "drop")
  rr <- rep(NA, NA_); rr[d$idarea] <- d$r; pp <- rep(NA, NA_); pp[d$idarea] <- d$p
  m1 <- moran.test(rr, lw, zero.policy = TRUE, na.action = na.omit, alternative = "two.sided")
  m2 <- moran.test(pp, lw, zero.policy = TRUE, na.action = na.omit, alternative = "two.sided")
  m3 <- moran.test(rr, nb_w_ref, zero.policy = TRUE, na.action = na.omit, alternative = "two.sided")
  c(m1$estimate[[1]], m1$p.value, m2$estimate[[1]], m2$p.value, m3$estimate[[1]], m3$p.value) }
fitrow <- function(m, arm, dat = C, lw = nb_w, note = "") { if (is.null(m)) return(data.frame(arm = arm, DIC = NA, WAIC = NA, pD = NA, moran_raw = NA, p_raw = NA, moran_poisson_scaled = NA, p_poisson_scaled = NA, moran_raw_queen_ref = NA, p_raw_queen_ref = NA, n_high = NA, n_low = NA, note = "fit failed"))
  mi <- resid_moran(m, dat, lw); re <- m$summary.random$idarea; ana <- sort(unique(dat$idarea))
  nh <- if (!is.null(re)) sum(re$`0.025quant`[ana] > 0) else NA; nl <- if (!is.null(re)) sum(re$`0.975quant`[ana] < 0) else NA
  data.frame(arm = arm, DIC = m$dic$dic, WAIC = m$waic$waic, pD = m$dic$p.eff, moran_raw = mi[1], p_raw = mi[2], moran_poisson_scaled = mi[3], p_poisson_scaled = mi[4], moran_raw_queen_ref = mi[5], p_raw_queen_ref = mi[6], n_high = nh, n_low = nl, note = note) }
COEF <- list(); FITS <- list(); SAVE <- list()
run_arm <- function(arm, covs, offset = "log(population + 1)", re = RE_M6, dat = C, family = "nbinomial", lw = nb_w, note = "", lincomb = NULL) {
  cat(sprintf("\n[%s] %s  (%d covariates, offset %s)\n", arm, note, length(covs), offset)); t0 <- Sys.time()
  m <- fit(covs, offset, re, dat, family, lincomb); cat(sprintf("   %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ct <- coeftab(m, arm); fr <- fitrow(m, arm, dat, lw, note)
  if (!is.null(ct)) { COEF[[arm]] <<- ct; print(ct[ct$credible == 1, c("covariate","IRR","lo","hi")], row.names = FALSE, digits = 4) }
  FITS[[arm]] <<- fr; print(fr[, c("DIC","WAIC","moran_raw","p_raw","moran_raw_queen_ref","n_high","n_low")], row.names = FALSE, digits = 4)
  write.csv(do.call(rbind, COEF), file.path(OUT, "coefficients_all_arms.csv"), row.names = FALSE); write.csv(do.call(rbind, FITS), file.path(OUT, "fit_all_arms.csv"), row.names = FALSE)
  invisible(m) }
CRED5 <- T2$covariate[T2$credible == 1]; cat("credible set as submitted:", CRED5, "\n")

# ---------------------------------------------------------------- 4. arms
ARMS <- Sys.getenv("HAV_ARMS", "all")
m0 <- run_arm("A0_baseline", Z_ORIG, note = "submitted specification refitted on this input frame")
m1 <- run_arm("A1_offset_pop", Z_ORIG, offset = "log(population)", note = "offset population without +1")
m2 <- run_arm("A2_agestd_offset", Z_ORIG, offset = "log(E_age)", note = "indirect age-standardised expected counts as offset")
m3 <- run_arm("A3_p2049", c(Z_ORIG, "p2049_z"), offset = "log(population)", note = "proportion aged 20-49 added")
m4 <- run_arm("A4_rescaled", c(Z_R, "popdens_z"), offset = "log(population)", note = "rescaled covariates (original forms) + log population density")
m5 <- run_arm("A5_rescaled_agestd", c(Z_R, "popdens_z"), offset = "log(E_age)", note = "rescaled + age-standardised offset (candidate revised principal)")
if (ARMS == "A5") { CO <- do.call(rbind, COEF); FI <- do.call(rbind, FITS); write.csv(CO, file.path(OUT, "coefficients_all_arms.csv"), row.names = FALSE); write.csv(FI, file.path(OUT, "fit_all_arms.csv"), row.names = FALSE); cat("\nHAV_ARMS=A5: stopped after A5\n"); sink(); quit(save = "no") }
m6 <- run_arm("A6_uniform_raw", c(Z_U, "popdens_z"), offset = "log(E_age)", note = "uniform form: continuous untransformed")
m7 <- run_arm("A7_uniform_log1p", c(Z_L, "popdens_z"), offset = "log(E_age)", note = "uniform form: continuous log(x+1)")
P1 <- paste0(CRED5, "_r"); P2 <- c(P1, "sex_ratio_r", "elderly_rate_r", "urban_pop_rate_r", "popdens_z")
m8a <- run_arm("A8_parsimonious_P1", P1, offset = "log(E_age)", note = "credible set only (rescaled)")
m8b <- run_arm("A8_parsimonious_P2", P2, offset = "log(E_age)", note = "credible set + sex ratio, elderly, urban share, density")
if (exists("m5") && !is.null(m5)) saveRDS(m5, file.path(OUT, "A5_fit.rds"))
# A9 restricted spatial regression: besag effect constrained orthogonal to district means of the covariates
rsr <- function(arm, covs, offset = "log(E_age)", note = "RSR") {
  Xd <- sapply(covs, function(v) tapply(C[[v]], C$idarea, mean)[as.character(1:NA_)]); Xd[is.na(Xd)] <- 0; A <- t(Xd)
  re <- paste0("+ f(idarea, model='besag', graph=g_main, scale.model=TRUE, constr=TRUE, extraconstr=list(A=A, e=rep(0, nrow(A))), hyper=pc_prec)",
               " + f(idtime, model='rw1', hyper=pc_prec) + f(idarea_time, model='iid', hyper=pc_prec)")
  assign("A", A, envir = .GlobalEnv); run_arm(arm, covs, offset, re, note = note) }
m9 <- rsr("A9_RSR_rescaled_agestd", c(Z_R, "popdens_z"), note = "restricted spatial regression (A5 covariates)")
m9b <- rsr("A9_RSR_baseline", Z_ORIG, offset = "log(population + 1)", note = "restricted spatial regression (submitted covariates)")
# A10 zero-inflated negative binomial
m10a <- run_arm("A10_ZINB_baseline", Z_ORIG, family = "zeroinflatednbinomial1", note = "ZINB type 1, submitted specification")
m10b <- run_arm("A10_ZINB_rescaled_agestd", c(Z_R, "popdens_z"), offset = "log(E_age)", family = "zeroinflatednbinomial1", note = "ZINB type 1, A5 specification")
cat(sprintf("\nshare of zero district-years: %.1f%%\n", 100 * mean(C$cases == 0)))
# A11 nested random-effect structures with residual Moran's I
RE_LIST <- list(M1 = "", M2 = "+ f(idarea, model='besag', graph=g_main, scale.model=TRUE, hyper=pc_prec)",
  M3 = "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym)",
  M4 = "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym) + f(idarea_time, model='iid', hyper=pc_prec)",
  M5 = "+ f(idarea, model='bym', graph=g_main, scale.model=TRUE, hyper=pc_bym) + f(idtime, model='rw1', hyper=pc_prec)", M6 = RE_M6)
for (nm in names(RE_LIST)) run_arm(paste0("A11_", nm, "_submitted"), Z_ORIG, re = RE_LIST[[nm]], note = paste("nested", nm, "submitted covariates"))
for (nm in names(RE_LIST)) run_arm(paste0("A11_", nm, "_A5"), c(Z_R, "popdens_z"), offset = "log(E_age)", re = RE_LIST[[nm]], note = paste("nested", nm, "A5 covariates"))
run_arm("A11_M0_nocov_submitted", character(0), re = RE_M6, note = "M6 random effects without covariates")
# A12 eight neighbourhood graphs (queen, rook, k-nearest 2-7) for A5
cents <- st_coordinates(st_centroid(st_geometry(shp_main)))
GR <- list(queen = nb, rook = poly2nb(shp_main, queen = FALSE, snap = 0.01)); for (k in 2:7) GR[[paste0("knn", k)]] <- knn2nb(knearneigh(cents, k = k), sym = TRUE)
for (gn in names(GR)) { nbg <- GR[[gn]]; g_main <- inla.read.graph(nb2mat(nbg, style = "B", zero.policy = TRUE)); lwg <- nb2listw(nbg, style = "W", zero.policy = TRUE)
  run_arm(paste0("A12_graph_", gn), c(Z_R, "popdens_z"), offset = "log(E_age)", lw = lwg, note = paste("graph", gn, "A5 covariates")) }
g_main <- inla.read.graph(g)
# A13 spatially varying coefficients: besag-structured local slopes for forest share and household wells
svc <- function(arm, var) { C$idarea_s <- C$idarea
  re <- paste0(RE_M6, " + f(idarea_s, ", var, ", model='besag', graph=g_main, scale.model=TRUE, hyper=pc_prec)")
  lc_args <- list(rep(1, NA_), diag(NA_)); names(lc_args) <- c(var, "idarea_s"); lc <- do.call(inla.make.lincombs, lc_args); names(lc) <- sprintf("d%03d", seq_len(NA_))
  m <- run_arm(arm, c(Z_R, "popdens_z"), offset = "log(E_age)", re = re, dat = C, note = paste("SVC on", var), lincomb = lc)
  if (!is.null(m)) { sl <- m$summary.lincomb.derived; sl <- sl[sprintf("d%03d", seq_len(NA_)), , drop = FALSE]
    loc <- data.frame(idarea = seq_len(NA_), local_slope = sl$mean, lo = sl$`0.025quant`, hi = sl$`0.975quant`, region = shp_main$region)
    loc$analysed <- loc$idarea %in% unique(C$idarea)
    write.csv(loc, file.path(OUT, paste0(arm, "_local_slopes.csv")), row.names = FALSE, fileEncoding = "UTF-8")
    hy <- m$summary.hyperpar; cat("   hyperparameters (slope field):\n"); print(hy[grep("idarea_s", rownames(hy)), , drop = FALSE], digits = 3)
    la <- loc[loc$analysed, ]; cat(sprintf("   total local slope (beta+b_i) IRR range %.3f-%.3f; analysed districts with CrI excluding 0: negative %d, positive %d; global %.3f (%.3f-%.3f)\n",
      exp(min(la$local_slope)), exp(max(la$local_slope)), sum(la$hi < 0), sum(la$lo > 0), exp(m$summary.fixed[var,"mean"]), exp(m$summary.fixed[var,"0.025quant"]), exp(m$summary.fixed[var,"0.975quant"]))) }
  invisible(m) }
svc("A13_SVC_forest", "forest_r"); svc("A13_SVC_gw_household", "gw_household_r")
# A14 ascertainment proxy
if ("phys_z" %in% names(C) && !any(is.na(C$phys_z))) run_arm("A14_physicians", c(Z_R, "popdens_z", "phys_z"), offset = "log(E_age)", note = "A5 + physicians per 1000")
# A15 exclude 2020-2021
C2 <- C %>% filter(year >= 2022) %>% mutate(idtime = year - 2021, idarea_time = row_number())
run_arm("A15_excl_2020_2021", c(Z_R, "popdens_z"), offset = "log(E_age)", dat = C2, note = "A5 specification, 2022-2024 only")

# A16 repeated fits with identical settings (numerical stability of the credible sets)
rep_tab <- list()
for (k in 1:5) { for (nm in c("A0","A5")) { covs <- if (nm == "A0") Z_ORIG else c(Z_R, "popdens_z"); off <- if (nm == "A0") "log(population + 1)" else "log(E_age)"
  m <- fit(covs, off); ct <- coeftab(m, paste0("A16_rep", k, "_", nm)); if (!is.null(ct)) { ct$DIC <- m$dic$dic; rep_tab[[length(rep_tab) + 1]] <- ct } } }
RT <- do.call(rbind, rep_tab); write.csv(RT, file.path(OUT, "repeated_fits_A0_A5.csv"), row.names = FALSE)
cat("\n[A16] repeated fits: credible-set identity across 5 fits\n")
for (nm in c("A0","A5")) { sets <- lapply(1:5, function(k) sort(RT$covariate[RT$arm == paste0("A16_rep", k, "_", nm) & RT$credible == 1]))
  cat(sprintf("  %s: identical in all fits = %s; DIC range %.2f-%.2f; set = %s\n", nm, all(sapply(sets, identical, sets[[1]])), min(RT$DIC[grepl(nm, RT$arm)]), max(RT$DIC[grepl(nm, RT$arm)]), paste(sets[[1]], collapse = ", "))) }

# ---------------------------------------------------------------- 5. summary
CO <- do.call(rbind, COEF); FI <- do.call(rbind, FITS)
write.csv(CO, file.path(OUT, "coefficients_all_arms.csv"), row.names = FALSE); write.csv(FI, file.path(OUT, "fit_all_arms.csv"), row.names = FALSE)
pers <- CO %>% filter(covariate %in% CRED5) %>% dplyr::select(arm, covariate, IRR, lo, hi, credible) %>% arrange(covariate, arm)
write.csv(pers, file.path(OUT, "credible_set_persistence.csv"), row.names = FALSE)
cat("\n==== persistence of the five submitted credible associations ====\n"); print(as.data.frame(pers), digits = 4, row.names = FALSE)
cat("\n==== fit indices ====\n"); print(FI[, c("arm","DIC","WAIC","moran_raw","p_raw","moran_raw_queen_ref","n_high","n_low")], digits = 4, row.names = FALSE)
cat("\ndone -", format(Sys.time()), "\n"); sink()
