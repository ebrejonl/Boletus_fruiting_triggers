# Code accompaning the manuscript "Predicting porcini: a decade of sporocarp monitoring reveals the meteorological triggers of Boletus edulis fruiting in central European beech forests"

library(data.table)
library(tidyverse)
library(zoo)
library(climateExtract)
library(scales)
library(patchwork)
library(cowplot)
library(glmmTMB)
library(parallel)
library(ggeffects)
renv::status()


# I extracted the environmental data from the GPS location this way: 
# Note: Unfortunately, since the collection site is a privately own site, public release of the GPS data will not be permitted. 
# However, the environmental and daily sporocarp counts data needed to fully reproduce all the analysis and the figures is shared.
# For transparency, I still share the code used to retrieve the environmental data below:

mydata=fread(file = "Data/Note_data_no_env_Nov21.csv") 
# centroid position 
meanlon=mean(as.numeric(as.character(mydata$lon)), na.rm = T)
meanlat=mean(as.numeric(as.character(mydata$lat)), na.rm = T)

#==========================================================================================================================================================#

# Environmental data

# One grid cell for the study area
gc()
#climate_data_temp <- extract_nc_value(first_year = 2015, last_year = 2024, local_file = FALSE, file_path = NULL, sml_chunk = "2011-2025", 
#                                      spatial_extent = c(meanlon - 0.001, meanlat - 0.001, meanlon + 0.001, meanlat + 0.001),
#                                      clim_variable = "mean temp",statistic = "mean",grid_size = 0.1,ecad_v = NULL,write_raster = TRUE,out = "Data/raster_mean_temp.tiff",return_data =TRUE)
#
#climate_data_precip <- extract_nc_value(first_year = 2015,last_year = 2024,local_file = FALSE,file_path = NULL,sml_chunk = "2011-2025",
#                                        spatial_extent = c(meanlon - 0.001, meanlat - 0.001, meanlon + 0.001, meanlat + 0.001),
#                                        clim_variable = "precipitation", statistic = "mean", grid_size = 0.1, ecad_v = NULL, write_raster = TRUE, out = "Data/raster_precip.tiff", return_data =TRUE)
#==========================================================================================================================================================#

# Format the data

gps <- data.frame(lon=meanlon, lat=meanlat)
gc()

# climate_data_temp = terra::rast("Data/raster_mean_temp.tiff")
# climate_data_precip = terra::rast("Data/raster_precip.tiff")
# 
# temp_small <- raster::extract(x = climate_data_temp, y = gps, method = "bilinear") # method bilinear is to get an estinated "better" depending on proximity to nearby grid cells and their values
# Precip_small <- raster::extract(x = climate_data_precip, y = gps, method = "bilinear") 
# 
# # Here all data is from our point, so we do not need coordinates. All we need is long format dates temp and precip
# 
# Climat_data=temp_small  |>  t() |> as_tibble() |> mutate(Date=colnames(temp_small)) |> rename("Temp"="V1") |> 
#   left_join(Precip_small |>  t() |> as_tibble() |> mutate(Date=colnames(temp_small)) |> rename("Precip"="V1")) |> filter(!Date=="ID") |> mutate(date=as.Date(Date))
# 
# Climat_data
#fwrite(file="Data/Climate.csv",Climat_data)
Climat_data <- fread(file="Data/Climate.csv") |> as_tibble() |> mutate(date=as.Date(Date))
gc()

Climat_roll <- Climat_data |> arrange(date) |> mutate(year = substr(date, 1, 4)) |> group_by(year)

# Formatting the sporocarp data per site and full database together
#==========================================================================================================================================================#



windows <- seq(2, 35, by = 3)  # 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35 

for (w in windows) Climat_roll <- Climat_roll |>
  mutate(!!paste0("rollT_", w) := rollmean(Temp,   w, align = "right", fill = NA),
         !!paste0("rollP_", w) := rollmean(Precip, w, align = "right", fill = NA))
Climat_roll <- ungroup(Climat_roll)

# Attach sporocarps per site, then subset to fruiting season
# full database
sporo=mydata |> group_by(date) |> summarise(Mnb=n()) |> mutate(date=as.Date(date))
# Site_A only (merged with below Site_A, a tiny small just next to it)
sporoSite_A=mydata |> mutate(Sitenum=parse_number(Site)) |> filter(Sitenum==23| Sitenum==3) |> group_by(date) |> summarise(Mnb=n()) |> mutate(date=as.Date(date))
# Site_B Site
sporoSite_B=mydata |> mutate(Sitenum=parse_number(Site)) |> filter(Sitenum==6) |> group_by(date) |> summarise(Mnb=n()) |> mutate(date=as.Date(date))

# add env
Full_ready    <- Climat_roll |> left_join(sporo,        by = "date") |> 
  mutate(month = substr(date, 6, 7)) |> filter(month %in% c("08","09","10","11")) |> mutate(Mnb = ifelse(is.na(Mnb), 0, Mnb))

Site_A_ready <- Climat_roll |> left_join(sporoSite_A, by = "date") |> 
  mutate(month = substr(date, 6, 7)) |> filter(month %in% c("08","09","10","11")) |> mutate(Mnb = ifelse(is.na(Mnb), 0, Mnb))

Site_B_ready  <- Climat_roll |> left_join(sporoSite_B,  by = "date") |> 
  mutate(month = substr(date, 6, 7)) |> filter(month %in% c("08","09","10","11")) |> mutate(Mnb = ifelse(is.na(Mnb), 0, Mnb))

all_data <- bind_rows(Full_ready    %>% mutate(site = "All"), Site_A_ready %>% mutate(site = "Site_A"),  Site_B_ready  %>% mutate(site = "Site_B"))

#fwrite(all_data, file = "Data/Note_data_Jun6_climate_counts_windows.csv")


# Then:
windows <- seq(2, 35, by = 3)  # 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35 
sites_to_run <- unique(all_data$site)
# All 64 combinations: temp window × precip window
combos <- expand.grid(win_T = windows, win_P = windows, stringsAsFactors = FALSE)

fit_one_combo <- function(i, data, combos) {
  wT <- combos$win_T[i]
  wP <- combos$win_P[i]
  d <- data.frame(Mnb  = data$Mnb,rT   = data[[paste0("rollT_", wT)]],rP   = data[[paste0("rollP_", wP)]],year = as.factor(data$year))
  d <- d[complete.cases(d), ]

  res_re <- tryCatch({
    m <- glmmTMB(Mnb ~ rP + poly(rT, 2, raw = TRUE) + (1 | year),data = d, family = nbinom1())
    coefs <- summary(m)$coefficients$cond
    coef_P  <- coefs["rP", "Estimate"];  se_P  <- coefs["rP", "Std. Error"];  pval_P  <- coefs["rP", "Pr(>|z|)"]
    coef_T1 <- coefs["poly(rT, 2, raw = TRUE)1", "Estimate"]; se_T1 <- coefs["poly(rT, 2, raw = TRUE)1", "Std. Error"]; pval_T1 <- coefs["poly(rT, 2, raw = TRUE)1", "Pr(>|z|)"]
    coef_T2 <- coefs["poly(rT, 2, raw = TRUE)2", "Estimate"]; se_T2 <- coefs["poly(rT, 2, raw = TRUE)2", "Std. Error"]; pval_T2 <- coefs["poly(rT, 2, raw = TRUE)2", "Pr(>|z|)"]
    peak_T  <- ifelse(coef_T2 < 0, -coef_T1 / (2 * coef_T2), NA)
    peak_ht <- NA
    if (!is.na(peak_T)) {
      nd      <- data.frame(rT = peak_T, rP = mean(d$rP, na.rm = TRUE), year = d$year[1])
      peak_ht <- as.numeric(predict(m, newdata = nd, type = "response", re.form = NA))
    }

    data.frame(aic = AIC(m),
               coef_P = coef_P, se_P = se_P, pval_P = pval_P,
               coef_T1 = coef_T1, se_T1 = se_T1, pval_T1 = pval_T1,
               coef_T2 = coef_T2, se_T2 = se_T2, pval_T2 = pval_T2,
               peak_T = peak_T, peak_ht = peak_ht, converged = TRUE)
  }, error = function(e) {
    data.frame(aic = NA,
               coef_P = NA, se_P = NA, pval_P = NA,
               coef_T1 = NA, se_T1 = NA, pval_T1 = NA,
               coef_T2 = NA, se_T2 = NA, pval_T2 = NA,
               peak_T = NA, peak_ht = NA, converged = FALSE)
  })

  cbind(combos[i, ], setNames(res_re, paste0("RE_", names(res_re))))
}
all_results <- list()

for (s in sites_to_run) {
  site_data <- all_data |> filter(site == s)
  res_list  <- mclapply(1:nrow(combos), fit_one_combo, data = site_data, combos = combos, mc.cores = 20)
  res  <- bind_rows(res_list)
  res$site  <- s
  all_results[[s]] <- res
}

results_all <- bind_rows(all_results) |>
  mutate(RE_dAIC = RE_aic - min(RE_aic, na.rm = TRUE), .by = site) |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))

#fwrite(results_all, "Data/window_screen_all_sites_grouped_year_rounded_pvalTemp.csv") # saved for supplementary Table 1


#~~~~~~~~~~~~~~~ Compare model performance for all window length combinations ~~~~~~~~~~~~~#

#  Panel A: AIC vs Precipitation window, temperature fixed at 20 
results_all= results_all |> mutate(Site=case_when(site=="Site_A"~"Site A", site=="Site_B"~"Site B", site=="All"~"All sites", .default = site))
results_all <- results_all |> mutate(Site = factor(Site, levels = c("All sites", "Site A", "Site B")))
p_precip <- results_all |>
  filter(win_T == 20) |>
  ggplot(aes(x = win_P, y = RE_dAIC, color = Site)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  labs(x = "Precipitation window length (days)", y = "ΔAIC",
       title = "Temperature window fixed at 20 days") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 16), plot.title = element_text(size=16))

#  Panel B: AIC vs Temperature window, precipitation fixed at 26 
p_temp <- results_all |>
  filter(win_P == 26) |>
  ggplot(aes(x = win_T, y = RE_dAIC, color = Site)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  labs(x = "Temperature window length (days)",
       title = "Precipitation window fixed at 26 days") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
axis.title = element_text(size = 16),
axis.title.y.left = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 16), plot.title = element_text(size=16))


# combine
aic_curves <-  p_temp +p_precip +
  plot_layout(ncol = 2, guides = "collect") +  plot_annotation( tag_levels = "A") & theme(plot.tag = element_text(size = 20, face = "plain"));aic_curves
#ggsave(aic_curves, filename = "Figures/AIC_vs_window_curves_Jun6.pdf", device=cairo_pdf, width = 14, height = 6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Figure 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Heatmap function 
# shared color range across all sites so visually comparable
results_all <- results_all |>
  group_by(site) |>
  mutate(RE_dAIC_scaled = scales::rescale(-RE_dAIC, to = c(0, 1))) |>
  ungroup()

aic_heatmap <- function(df, site_name, show_y = FALSE, show_x_title = FALSE) {
  p <- ggplot(df, aes(x = factor(win_P), y = factor(win_T), fill = RE_dAIC_scaled)) +
    geom_raster(alpha=1) +
    scale_fill_viridis_c(option = "H",direction=1,begin = 0.0, end = 1,
      name = "Relative model fit\n(scaled -AIC)") +
    labs(x = NULL, y = NULL, title = site_name) +
    theme_classic() +
    theme(
      axis.line   = element_blank(),
      plot.title = element_text(size = 17, color = "black", hjust = 0.5),
      axis.text.x= element_text(color = "black", size = 14),
      axis.ticks = element_blank(),
      legend.title= element_text(size = 16, color = "black", vjust = 0.5),
      legend.text = element_text(size = 16, color = "black"))
  if (show_y) {
    p <- p + ylab("Temperature window length (days)") +
      theme(axis.text.y = element_text(color = "black", size = 14),
            axis.title.y.left = element_text(size = 16, color = "black"))
  } else {
    p <- p + theme(axis.text.y = element_blank())
  }
  if (show_x_title) {
    p <- p + xlab("Precipitation window length (days)") +
      theme(axis.title.x.bottom = element_text(size = 16, color = "black"))
  }
  return(p)
}

h_all <- aic_heatmap(results_all |> filter(site == "All"), "All sites", show_y = TRUE)
h_Site_A <- aic_heatmap(results_all |> filter(site == "Site_A"), "Site A",   show_x_title = TRUE)
h_Site_B  <- aic_heatmap(results_all |> filter(site == "Site_B"),  "Site B")

aic_panel <- h_all + plot_spacer() + h_Site_A + plot_spacer() + h_Site_B + plot_layout(widths = c(0.25, 0.001, 0.25, 0.001, 0.25), guides = "collect");aic_panel
#ggsave(aic_panel, filename = "Figures/AIC_heatmaps_RE_Jun6.pdf", width = 16, height = 5)
gc()

#  Fit the best model for each site 
# Using window: 20 days for temperature, 26 days for precipitation
# With year as random effect, nbinom1 family, like in the window exploration (identical as the best window one from before)
all_sub <- all_data |> filter(site == "All", !is.na(rollT_20), !is.na(rollP_26))
global_mean_T <- mean(all_sub$rollT_20, na.rm = TRUE)
global_mean_P <- mean(all_sub$rollP_26, na.rm = TRUE)


fit_best <- function(df) {
  d <- df |> filter(!is.na(rollT_20), !is.na(rollP_26)) |>
    mutate(rT = rollT_20, rP = rollP_26, year = as.factor(year))
  glmmTMB(Mnb ~ rP + poly(rT, 2, raw = TRUE) + (1 | year), data = d, family = nbinom1())
}

m_all     <- fit_best(all_data |> filter(site == "All"))
m_Site_A <- fit_best(all_data |> filter(site == "Site_A"))
m_Site_B  <- fit_best(all_data |> filter(site == "Site_B"))

# Get optimal temperature
get_peak <- function(m) {
  coefs <- summary(m)$coefficients$cond
  b1 <- coefs["poly(rT, 2, raw = TRUE)1", "Estimate"]
  b2 <- coefs["poly(rT, 2, raw = TRUE)2", "Estimate"]
  -b1 / (2 * b2)  # vertex of the quadratic
}

peak_all     <- get_peak(m_all)
peak_Site_A <- get_peak(m_Site_A)
peak_Site_B  <- get_peak(m_Site_B)

# Get prediction for temperature to showcase the curvature/quadratic nature of the effect
pred_T_all <- ggpredict(m_all,     terms = "rT [all]", condition = c(rP = global_mean_P)) |> as_tibble() |> mutate(site = "All sites")
pred_T_Site_A <- ggpredict(m_Site_A, terms = "rT [all]", condition = c(rP = global_mean_P)) |> as_tibble() |> mutate(site = "Site A")
pred_T_Site_B  <- ggpredict(m_Site_B,  terms = "rT [all]", condition = c(rP = global_mean_P)) |> as_tibble() |> mutate(site = "Site B")

pred_T <- bind_rows(pred_T_all, pred_T_Site_A, pred_T_Site_B) |>
  group_by(site) |>
  mutate(max_pred= max(predicted, na.rm = TRUE),
  scaled = predicted / max_pred, scaled_lo = (predicted - std.error) / max_pred,scaled_hi  = (predicted + std.error) / max_pred) |>
  ungroup()

# Get prediction for temperature precipitation curves (at the mean study temperature)
pred_P_all     <- ggpredict(m_all,     terms = "rP [all]", condition = c(rT = global_mean_T)) |> as_tibble() |> mutate(site = "All sites")
pred_P_Site_A <- ggpredict(m_Site_A, terms = "rP [all]", condition = c(rT = global_mean_T)) |> as_tibble() |> mutate(site = "Site A")
pred_P_Site_B  <- ggpredict(m_Site_B,  terms = "rP [all]", condition = c(rT = global_mean_T)) |> as_tibble() |> mutate(site = "Site B")

pred_P <- bind_rows(pred_P_all, pred_P_Site_A, pred_P_Site_B) |>
  group_by(site) |>
  mutate(max_pred= max(predicted, na.rm = TRUE),
    scaled = predicted / max_pred, scaled_lo = (predicted - std.error) / max_pred, scaled_hi = (predicted + std.error) / max_pred) |>
  ungroup()


p_temp <- ggplot(pred_T, aes(x = x, y = scaled, color = site, fill = site)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = peak_all, linetype = "dashed", color = "grey40") +
  labs(x = "Temperature (°C)", y = NULL) +
  theme_classic() +
  theme( axis.text    = element_text(color = "black", size = 14),
    axis.title.y.left   = element_text(size = 16, color = "black"),
    axis.title.x.bottom = element_text(size = 16, color = "black"),
    legend.title        = element_blank(),
    legend.text         = element_text(size = 12))

#  Precipitation plot 
p_precip <- ggplot(pred_P, aes(x = x, y = scaled, color = site, fill = site)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Precipitation (mm/day)", y ="Relative sporocarp production") +
  theme_classic() +
  theme( axis.text = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(size = 16, color = "black"),
    axis.title.x.bottom = element_text(size = 16, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12))

# Combine 
p_temp <- p_temp + labs(x = "Temperature, 20-day mean (°C)")
p_precip <- p_precip + labs(x = "Precipitation, 26-day mean (mm/day)")

gc()



## Heat map conditions fruiting vs conditions occurence for all years together, three different sites. 

kde_input <- all_data |>
  filter(!is.na(rollT_20), !is.na(rollP_26)) |>
  mutate(rT = rollT_20,rP = rollP_26,site = factor(site, levels = c("All", "Site_A", "Site_B")))

# Observed climate 
obs_range <- kde_input |>
  group_by(site) |>
  summarise(T_min = min(rT), T_max = max(rT),P_min = min(rP), P_max = max(rP))

# Fruiting days only
fruiting <- kde_input |> filter(Mnb > 0)

# KDE on a common grid
T_lims <- range(kde_input$rT)
P_lims <- range(kde_input$rP)

compute_kde <- function(df, n = 300) {
  # Tiny jitter to avoid ties that can trip up KDE
  rP_j <- abs(jitter(df$rP, amount = 0.01))
  rT_j <- abs(jitter(df$rT, amount = 0.01))
  kde <- MASS::kde2d(rP_j, rT_j, n = n, lims = c(P_lims, T_lims))
  tidyr::crossing(T_val = kde$y, P_val = kde$x) |> mutate(density = as.vector(kde$z))
}

### Combined

# Observed conditions (all days) only "All" site
kde_obs <- kde_input |>
  filter(site == "All") |>
  reframe(compute_kde(pick(everything()))) |>
  mutate(panel = "All sites\n(all days)")

# Fruiting conditions, all three sites
kde_fruit <- fruiting |>
  group_by(site) |>
  reframe(compute_kde(pick(everything()))) |>
  ungroup() |>
  mutate(panel = case_when(site == "All"     ~  "All sites\n(fruiting days only)",
    site == "Site_A" ~ "Site A\n(fruiting days only)",
    site == "Site_B"  ~ "Site B\n(fruiting days only)"))

kde_combined <- bind_rows(kde_obs, kde_fruit) |>
  mutate(panel = factor(panel, levels = c( "All sites\n(all days)", "All sites\n(fruiting days only)",
      "Site A\n(fruiting days only)", "Site B\n(fruiting days only)" )), scaled_density = scales::rescale(density, to = c(0, 1)))

#  Plot
p_combined <- ggplot(kde_combined, aes(x = P_val, y = T_val, fill = scaled_density)) +
  geom_raster() +
  scale_fill_viridis_c(option = "H", 
    na.value = 0, begin = 0, end = 1,
    name = "Scaled\ndensity") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Precipitation, 26-day mean (mm/day)",
       y = "Temperature, 20-day mean (°C)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme( axis.line            = element_blank(),
    axis.text            = element_text(color = "black", size = 14),
    axis.title.x.bottom  = element_text(size = 16, color = "black"),
    axis.title.y.left    = element_text(size = 16, color = "black"),
    strip.text.x.top     = element_text(color = "black", size = 16),
    strip.background     = element_blank(),
    legend.title         = element_text(color = "black", size = 16),
    legend.text          = element_text(color = "black", size = 16))

p_combined

# with the models
model_curves_clean <- p_precip + theme(legend.text = element_text(size=16))+ p_temp + theme(legend.text = element_text(size=16))+
  plot_layout(ncol = 2, guides = "collect")

pfull <- p_combined / plot_spacer() / model_curves_clean + plot_layout(heights = c(0.4, 0.01, 0.5)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))
gc()

### panel A to F version :
kde_plot <- function(lv) {
  ggplot(filter(kde_combined, panel == lv), aes(P_val, T_val, fill = scaled_density)) +
    geom_raster() +
    scale_fill_viridis_c(option = "H", limits = c(0, 1),       # one comparable scale -> one legend
                         na.value = "transparent", begin = 0, end = 1, name = "Scaled\ndensity") +
    labs(title = lv, x = "Precipitation, 26-day mean (mm/day)",
         y = "Temperature, 20-day mean (°C)") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text  = element_text(color = "black", size = 14),
          axis.title = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, color = "black", hjust = 0.5),
          legend.title = element_text(size = 16, color = "black"),
          legend.text  = element_text(size = 16, color = "black"))
}
kde_list <- Map(\(lv, tg) kde_plot(lv) + labs(tag = tg),
                levels(kde_combined$panel), c("A", "B", "C", "D"))

kde_row <- wrap_plots(kde_list, nrow = 1,
                      guides = "collect",
                      axes = "collect_y", 
                      axis_titles = "collect")      

curve_row <- wrap_plots(p_precip + labs(tag = "E") + theme(legend.text = element_text(size = 16)),plot_spacer(),p_temp   + labs(tag = "F") + theme(legend.text = element_text(size = 16)),ncol = 3, widths = c(1, 0.05, 1), guides = "collect")

pfull <- ((kde_row / plot_spacer() / curve_row) +
          plot_layout(heights = c(0.4, 0.05, 0.5))) &
  theme(plot.tag = element_text(size = 22));pfull

ggsave(pfull, filename = "Figures/KDE_and_models_A_to_F_Jun8.pdf", width = 16, height = 10)

gc()




#~~~~~~~~~~~~~~~~~ Figure 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Full_ready


# fruiting season
rect_data <- Full_ready  |> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  dplyr::filter(!is.na(Mnb)) |> 
  dplyr::group_by(year) |> arrange(date) |> filter(Mnb>0) |> 
  dplyr::summarise(xmin = min(Dm), xmax = max(Dm), .groups = "drop")

rect_data |> mutate(length=xmax-xmin+1) |> pull(length) |> summary()


# Get a symmetric y-range so both sides are comparable
max_precip <- max(Full_ready$Precip, na.rm = TRUE)
max_M  <- max(Full_ready$Mnb,    na.rm = TRUE)

P=Full_ready |> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  ggplot(aes(x = Dm)) +
  geom_col(aes(y = Precip), fill = "#49a2f5ff") +
  geom_col(aes(y = -Mnb), fill = "black") +
  facet_wrap( ~ year, strip.position = "left", ncol = 1 ) +  
 scale_y_continuous(limits = c(-max_M, max_precip),                
  breaks = c(0, 40), name   = "Precipitation (mm)", sec.axis = sec_axis(trans = ~ ., breaks = c(0, -50),  
  labels = function(x) abs(round(x)), name = "Sporocarp number"))+
#sec.axis = sec_axis(trans   = ~ ., breaks = -pretty(c(0, max_M)), labels  = function(x) abs(x), name= "Sporocarp number"))+
  theme_classic() + theme(legend.position = "none",strip.placement = "outside",strip.text.x.bottom = element_text(margin = margin(t = 6)),
    axis.text.x   = element_text(margin = margin(t = 4)),axis.title.x  = element_text(margin = margin(t = 8)), 
    plot.margin = margin(5, 5, 10, 5),
    axis.title.y.left = element_text(color = "#49a2f5ff",size=16),axis.text.y.left  = element_text(color = "#49a2f5ff"),
    axis.line.y.left = element_line(color = "#49a2f5ff"),
    axis.ticks.y.left = element_line(color="#49a2f5ff"), 
    axis.title.x.bottom = element_blank(),
    axis.title.y.right = element_text(color = "black", size=16),
   strip.text.y.left = element_blank(),
    axis.text.y.right  = element_text(color = "black"),
     strip.background.y = element_blank()) #;P

P = P + 
  geom_rect(data = rect_data,
          aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = -max_M, ymax = 0),
          inherit.aes = FALSE, fill = "grey60", alpha = 0.2)

gc()
max_temp <- max(Full_ready$Temp, na.rm = TRUE)
max_M    <- max(Full_ready$Mnb,na.rm = TRUE)

T= Full_ready|> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  ggplot(aes(x = Dm)) +
  geom_col(aes(y = Temp), fill = "#de450f") +
  geom_col(aes(y = -Mnb), fill = "black") +
    facet_wrap( ~ year, strip.position = "left", ncol = 1) +  
 scale_y_continuous(limits = c(-max_M, max_precip), breaks = c(0, 30), name = "Temperature (C)",
  sec.axis = sec_axis(trans = ~ ., breaks = c(0, -50),
                      labels = function(x) abs(round(x)),
                      name = "Sporocarp number"))  +   
  theme_classic() + theme(legend.position = "none",strip.placement = "outside",strip.text.x.bottom = element_text(margin = margin(t = 6)),
    axis.text.x   = element_text(margin = margin(t = 4)),axis.title.x  = element_text(margin = margin(t = 8)), plot.margin   = margin(5, 5, 10, 5),
    axis.title.y.left = element_text(color = "#de450f", size=16),axis.text.y.left  = element_text(color = "#de450f"),
    axis.line.y.left = element_line(color = "#de450f"),
    axis.ticks.y.left = element_line(color="#de450f"), 
    axis.title.y.right = element_text(color = "black", size=16),axis.title.x.bottom = element_blank(),
    axis.text.y.right  = element_text(color = "black"), strip.text.y.left = element_text(size=16, color="black"), strip.background.y = element_blank())#;T
#adding fruiting season
T = T + 
  geom_rect(data = rect_data,
          aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = -max_M, ymax = 0),
          inherit.aes = FALSE, fill = "grey60", alpha = 0.2)
gc()

# Density and boxplot of observed values during the study period every year
densT=Full_ready |> ggplot()+
  geom_density(mapping=aes(x=Temp), fill="#de450f", color="#de450f",alpha=0.8, trim=TRUE )+
  geom_boxplot(mapping=aes(x=Temp, y=0.15), fill="#de450f", 
color="#de450f",alpha=0.4, outliers = F, width=0.05)+
  theme_void()+
  facet_wrap(~year, strip.position = "left", ncol = 1)+
  xlab("Temperature (C)")+
  theme(axis.title.x.bottom = element_text(color="#de450f", size=16),
        axis.text.x.bottom =element_text(color="#de450f"),
        axis.line.x.bottom = element_line(color="#de450f"),
      strip.background = element_blank(), strip.text = element_blank() );densT


densP=Full_ready|> ggplot()+
  geom_density(mapping=aes(x=Precip), fill="#49a2f5ff", color="#49a2f5ff",alpha=0.8, trim=TRUE )+
  geom_boxplot(mapping=aes(x=Precip, y=0.8), fill="#49a2f5ff", 
color="#49a2f5ff",alpha=0.2, outliers = F, width=0.35)+
  facet_wrap(~year, strip.position = "left", ncol = 1)+
  xlab("Precipitation (mm)")+
    theme_void()+
  theme(axis.title.x.bottom = element_text(color="#49a2f5ff", size=16),
        axis.text.x.bottom =element_text(color="#49a2f5ff"),
        axis.line.x.bottom = element_line(color="#49a2f5ff"),
      strip.background = element_blank(), strip.text = element_blank() );densP

bump_main <- theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16))
bump_dens <- theme(axis.text.x.bottom  = element_text(size = 14),axis.title.x.bottom = element_text(size = 16), plot.margin = margin(5.5, 16, 5.5, 5.5))

full_long3 <- plot_grid(T + bump_main, NULL, densT + bump_dens, NULL,
                        P + bump_main, NULL, densP + bump_dens,
                        ncol = 7,
                        rel_widths = c(0.4, -0.001, 0.1, 0.01, 0.4, 0.001, 0.1),
                        align = "v", axis = "b", scale = 1,
                        labels = c("A", "", "", "", "B", ""),
                        label_size = 22, label_fontface = "plain")
full_long3
#ggsave(full_long3, filename = "Figures/Figure1_Jun4.pdf", width = 17, height = 8.5)