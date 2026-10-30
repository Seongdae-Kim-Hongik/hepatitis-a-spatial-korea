# Figures for the R1 revision from the revised principal model (A5): Figure 2 (zero-boundary classes), S3, S4; district effects table.
suppressMessages({library(sf); library(dplyr); library(ggplot2); library(patchwork); library(ggspatial); library(INLA)})
R <- Sys.getenv("HAV_R1_DIR", file.path(getwd(), "revision_R1")); RES <- file.path(R, "results"); FIG <- file.path(R, "figures"); dir.create(FIG, showWarnings = FALSE)
SHP <- file.path(Sys.getenv("HAV_DATA_DIR", file.path(getwd(), "data")), "final.shp")
m <- readRDS(file.path(RES, "A5_fit.rds")); C <- read.csv(file.path(RES, "analysis_frame_R1.csv"), check.names = FALSE, fileEncoding = "UTF-8")
shp <- st_read(SHP, quiet = TRUE); st_crs(shp) <- 5179
shp$region <- gsub("[[:space:]]+", "", as.character(shp$region)); shp$region[shp$region == "인천시미추홀구"] <- "인천시남구"
islands <- c("인천시옹진군","전라남도완도군","전라남도진도군","경상남도거제시","경상남도남해군","경상북도울릉군")
main <- shp[!shp$region %in% islands, ]; nb <- spdep::poly2nb(main, snap = 0.01); iso <- which(spdep::card(nb) == 0); if (length(iso)) main <- main[-iso, ]
NA_ <- nrow(main); stopifnot(NA_ == 223)
re <- m$summary.random$idarea[1:NA_, ]            # BYM: first NA_ rows = u + v
ana <- sort(unique(C$idarea))
eff <- data.frame(idarea = 1:NA_, region = main$region, mean = re$mean, lo = re$`0.025quant`, hi = re$`0.975quant`) %>%
  mutate(analysed = idarea %in% ana, class = case_when(!analysed ~ "Not analyzed", lo > 0 ~ "High", hi < 0 ~ "Low", TRUE ~ "Indeterminate"), mean = ifelse(analysed, mean, NA))
write.csv(eff %>% arrange(desc(mean)), file.path(RES, "A5_district_effects.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("classes:\n"); print(table(eff$class))
m2 <- main %>% left_join(eff, by = "region") %>% mutate(class = factor(class, levels = c("High","Low","Indeterminate","Not analyzed")))
# zero-boundary classes (reviewer AR-8 / editor m6)
br <- c(-Inf, -0.5, -0.25, 0, 0.25, 0.5, Inf); lab <- c("< −0.50", "−0.50 to −0.25", "−0.25 to 0", "0 to 0.25", "0.25 to 0.50", "> 0.50")
m2$effc <- as.character(cut(m2$mean, br, labels = lab)); m2$effc[is.na(m2$effc)] <- "Not analyzed"; m2$effc <- factor(m2$effc, levels = c(lab, "Not analyzed"))
cat("effect range:", range(m2$mean, na.rm = TRUE), "\n"); print(table(m2$effc))
th <- theme_void(base_size = 9) + theme(legend.position = "right", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9, hjust = 0))
north <- function(h = 0.8, w = 0.6, ts = 6) annotation_north_arrow(location = "tr", height = unit(h, "cm"), width = unit(w, "cm"), style = north_arrow_orienteering(text_size = ts))
pal <- setNames(c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b", "white"), levels(m2$effc))
pA <- ggplot(m2) + geom_sf(aes(fill = effc), colour = "grey40", linewidth = 0.1) + scale_fill_manual(values = pal, name = "Combined district\neffect (log RR)", drop = FALSE) + th + ggtitle("(A)")
pB <- ggplot(m2) + geom_sf(aes(fill = class), colour = "grey40", linewidth = 0.1) +
  scale_fill_manual(values = c(High = "#b2182b", Low = "#2166ac", Indeterminate = "grey85", `Not analyzed` = "white"), name = "Residual risk", drop = FALSE) + th + ggtitle("(B)") +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) + north()
ggsave(file.path(FIG, "Figure2.png"), pA + pB, width = 9, height = 4.5, dpi = 600, bg = "white")
# S3 observed / predicted / residual (A5 fitted includes random effects)
C$fit <- m$summary.fitted.values$mean[seq_len(nrow(C))]
ag <- C %>% group_by(region) %>% summarise(obs = sum(cases) / sum(population) * 1e5, pred = sum(fit) / sum(population) * 1e5, .groups = "drop") %>% mutate(res = obs - pred)
s1 <- main %>% left_join(ag, by = "region")
b5 <- unique(quantile(s1$obs, seq(0, 1, .2), na.rm = TRUE)); b5[1] <- min(b5[1], min(s1$pred, na.rm = TRUE)); b5[length(b5)] <- max(b5[length(b5)], max(s1$pred, na.rm = TRUE))
l5 <- sprintf("%.1f to %.1f", head(b5, -1), tail(b5, -1))
fz <- function(v, b, l) { o <- as.character(cut(v, b, include.lowest = TRUE, labels = l)); o[is.na(o)] <- "Not analyzed"; factor(o, levels = c(l, "Not analyzed")) }
rb <- c(-Inf, -1.7, -0.8, 0, 0.8, 1.7, Inf); rl <- c("< −1.7", "−1.7 to −0.8", "−0.8 to 0.0", "0.0 to 0.8", "0.8 to 1.7", "> 1.7")
s1$obsc <- fz(s1$obs, b5, l5); s1$predc <- fz(s1$pred, b5, l5); s1$resc <- fz(s1$res, rb, rl)
yv <- setNames(c(RColorBrewer::brewer.pal(5, "YlOrRd"), "white"), c(l5, "Not analyzed")); rv <- setNames(c(rev(RColorBrewer::brewer.pal(6, "RdBu")), "white"), c(rl, "Not analyzed"))
q1 <- ggplot(s1) + geom_sf(aes(fill = obsc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = yv, name = "per 100,000", drop = FALSE) + th + ggtitle("(A) Observed")
q2 <- ggplot(s1) + geom_sf(aes(fill = predc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = yv, name = "per 100,000", drop = FALSE) + th + ggtitle("(B) Predicted")
q3 <- ggplot(s1) + geom_sf(aes(fill = resc), colour = "grey40", linewidth = 0.08) + scale_fill_manual(values = rv, name = "Residual\n(per 100,000)", drop = FALSE) + th + ggtitle("(C) Residual") +
  annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.5) + north(0.6, 0.45, 5)
ggsave(file.path(FIG, "FigureS3.png"), q1 + q2 + q3, width = 12, height = 4.2, dpi = 600, bg = "white")
# S4 PIT
cp <- data.frame(region = C$region, year = C$year, cpo = m$cpo$cpo, pit = m$cpo$pit, failure = m$cpo$failure); write.csv(cp, file.path(RES, "A5_cpo_pit.csv"), row.names = FALSE)
cat(sprintf("CPO failures %d; CPO<0.001: %d; mean PIT %.3f\n", sum(cp$failure > 0), sum(cp$cpo < 0.001), mean(cp$pit)))
pS4 <- ggplot(cp, aes(pit)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#3182bd", colour = "white") + geom_hline(yintercept = nrow(cp) / 20, linetype = 2) + labs(x = "Probability integral transform (PIT)", y = "Count") + theme_classic(base_size = 11)
ggsave(file.path(FIG, "FigureS4.png"), pS4, width = 6, height = 4, dpi = 600, bg = "white")
hy <- m$summary.hyperpar; write.csv(hy, file.path(RES, "A5_hyperparameters.csv")); print(round(hy[, c("0.5quant","0.025quant","0.975quant")], 4))
cat("done\n")
