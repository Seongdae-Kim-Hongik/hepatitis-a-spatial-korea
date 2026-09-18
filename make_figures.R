# =============================================================================
# Figures for "Spatial Clustering of Hepatitis A in South Korea, 2020-2024" (v2.1.0)
# Draws Figure 2 and Figures S2-S4 from the outputs of HAV_spatial_reproducible.R.
# (Figure 1 and Figure S1 are descriptive maps of the public KDCA counts.)
# Needs results/m6_outputs.rds (written by a RUN_EXTENDED=true run; not committed
# because it is a binary file), results/getis_ord_gi.csv, results/cpo_pit_diagnostics.csv
# and the district shapefile under HAV_DATA_DIR.
# Run:  HAV_DATA_DIR=/path/to/data Rscript make_figures.R
# =============================================================================
suppressMessages({library(sf); library(dplyr); library(ggplot2); library(patchwork); library(ggspatial); library(classInt)})
DATA <- Sys.getenv("HAV_DATA_DIR", unset = file.path(getwd(), "data"))
RES  <- Sys.getenv("HAV_OUTPUT_DIR", unset = file.path(getwd(), "results"))
FIG  <- Sys.getenv("HAV_FIGURE_DIR", unset = file.path(getwd(), "figures")); dir.create(FIG, showWarnings = FALSE)
x <- readRDS(file.path(RES, "m6_outputs.rds")); ic <- x$ic
shp_file <- c(file.path(DATA, "districts.shp"), file.path(DATA, "final.shp")); shp_file <- shp_file[file.exists(shp_file)][1]
shp <- st_read(shp_file, quiet = TRUE) %>% mutate(region = gsub("\\s+", "", region), region = ifelse(region == "인천시미추홀구", "인천시남구", region))
if (is.na(st_crs(shp))) st_crs(shp) <- 5179          # the boundary file carries no CRS; coordinates are Korea 2000 Unified CS
main <- shp %>% filter(region %in% x$shp_region)
ana <- unique(ic$region)                              # districts that enter the likelihood
th <- theme_void(base_size = 9) + theme(legend.position = "right", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9, hjust = 0))
minus <- function(s) gsub("-", "−", s)
north <- function(h = 0.8, w = 0.6, ts = 6) annotation_north_arrow(location = "tr", height = unit(h, "cm"), width = unit(w, "cm"), style = north_arrow_orienteering(text_size = ts))

# ---- Figure 2: combined district effect and risk classification (3 districts without covariate data are not classified)
re <- data.frame(region = x$shp_region, eff = x$re$mean, lo = x$re$`0.025quant`, hi = x$re$`0.975quant`) %>%
  mutate(class = case_when(!(region %in% ana) ~ "Not analyzed", lo > 0 ~ "High", hi < 0 ~ "Low", TRUE ~ "Indeterminate"), eff = ifelse(region %in% ana, eff, NA))
m2 <- main %>% left_join(re, by = "region") %>% mutate(class = factor(class, levels = c("High", "Low", "Indeterminate", "Not analyzed")))
br <- classIntervals(m2$eff[!is.na(m2$eff)], n = 6, style = "quantile")$brks; br[1] <- floor(min(m2$eff, na.rm = TRUE) * 100) / 100; br[7] <- ceiling(max(m2$eff, na.rm = TRUE) * 100) / 100
lab <- minus(sprintf("%.2f to %.2f", head(br, -1), tail(br, -1)))
m2$effc <- as.character(cut(m2$eff, br, include.lowest = TRUE, labels = lab)); m2$effc[is.na(m2$effc)] <- "Not analyzed"; m2$effc <- factor(m2$effc, levels = c(lab, "Not analyzed"))
pal <- setNames(c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b", "white"), levels(m2$effc))
pA <- ggplot(m2) + geom_sf(aes(fill = effc), colour = "grey40", linewidth = 0.1) + scale_fill_manual(values = pal, name = "Combined district\neffect (log RR)", drop = FALSE) + th + ggtitle("(A)")
pB <- ggplot(m2) + geom_sf(aes(fill = class), colour = "grey40", linewidth = 0.1) +
  scale_fill_manual(values = c(High = "#b2182b", Low = "#2166ac", Indeterminate = "grey85", `Not analyzed` = "white"), name = "Residual risk", drop = FALSE) + th + ggtitle("(B)") +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) + north()
ggsave(file.path(FIG, "Figure2.png"), pA + pB, width = 9, height = 4.5, dpi = 600, bg = "white")

# ---- Figure S3: observed / predicted / residual incidence per 100,000 person-years
ic$fit <- x$M6_fitted
ag <- ic %>% group_by(region) %>% summarise(obs = sum(cases) / sum(population) * 1e5, pred = sum(fit) / sum(population) * 1e5, .groups = "drop") %>% mutate(res = obs - pred)
s1 <- main %>% left_join(ag, by = "region")
b5 <- unique(quantile(s1$obs, seq(0, 1, .2), na.rm = TRUE)); b5[1] <- min(b5[1], min(s1$pred, na.rm = TRUE)); b5[6] <- max(b5[6], max(s1$pred, na.rm = TRUE))
l5 <- sprintf("%.1f to %.1f", head(b5, -1), tail(b5, -1))
fz <- function(v, b, l) { o <- as.character(cut(v, b, include.lowest = TRUE, labels = l)); o[is.na(o)] <- "Not analyzed"; factor(o, levels = c(l, "Not analyzed")) }
rb <- c(-Inf, -1.7, -0.8, 0, 0.8, 1.7, Inf); rl <- minus(c("< -1.7", "-1.7 to -0.8", "-0.8 to 0.0", "0.0 to 0.8", "0.8 to 1.7", "> 1.7"))
s1$obsc <- fz(s1$obs, b5, l5); s1$predc <- fz(s1$pred, b5, l5); s1$resc <- fz(s1$res, rb, rl)
yv <- setNames(c(RColorBrewer::brewer.pal(5, "YlOrRd"), "white"), c(l5, "Not analyzed")); rv <- setNames(c(rev(RColorBrewer::brewer.pal(6, "RdBu")), "white"), c(rl, "Not analyzed"))
q1 <- ggplot(s1) + geom_sf(aes(fill = obsc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = yv, name = "per 100,000", drop = FALSE) + th + ggtitle("(A) Observed")
q2 <- ggplot(s1) + geom_sf(aes(fill = predc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = yv, name = "per 100,000", drop = FALSE) + th + ggtitle("(B) Predicted")
q3 <- ggplot(s1) + geom_sf(aes(fill = resc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = rv, name = "Residual\n(per 100,000)", drop = FALSE) + th + ggtitle("(C) Residual") +
  annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.5) + north(0.6, 0.45, 5)
ggsave(file.path(FIG, "FigureS3.png"), q1 + q2 + q3, width = 12, height = 4.2, dpi = 600, bg = "white")

# ---- Figure S2: Getis-Ord Gi* on the observed incidence of all 223 graph districts
gi <- read.csv(file.path(RES, "getis_ord_gi.csv"), fileEncoding = "UTF-8")
s2 <- main %>% left_join(gi, by = "region") %>% mutate(gc = cut(Gi_z, c(-Inf, -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, Inf),
  labels = c("≤ −2.58", "−2.58 to −1.96", "−1.96 to −1.65", "−1.65 to 1.65", "1.65 to 1.96", "1.96 to 2.58", "≥ 2.58")))
pS2 <- ggplot(s2) + geom_sf(aes(fill = gc), colour = "grey40", linewidth = 0.1) +
  scale_fill_manual(values = c("#2166ac", "#67a9cf", "#d1e5f0", "#f0f0f0", "#fddbc7", "#ef8a62", "#b2182b"), name = "Getis-Ord Gi*\n(z score)", drop = FALSE) + th +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) + north()
ggsave(file.path(FIG, "FigureS2.png"), pS2, width = 6, height = 6.5, dpi = 600, bg = "white")

# ---- Figure S4: PIT histogram
cp <- read.csv(file.path(RES, "cpo_pit_diagnostics.csv"))
pS3 <- ggplot(data.frame(pit = cp$pit), aes(pit)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#3182bd", colour = "white") +
  geom_hline(yintercept = nrow(cp) / 20, linetype = 2) + labs(x = "Probability integral transform (PIT)", y = "Count") + theme_classic(base_size = 11)
ggsave(file.path(FIG, "FigureS4.png"), pS3, width = 6, height = 4, dpi = 600, bg = "white")
cat("figures written to", FIG, "\n"); print(table(m2$class)); cat("Gi* hot/cold (95%):", sum(gi$Gi_z > 1.96), sum(gi$Gi_z < -1.96), "\n")
