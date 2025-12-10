# Code accompaning the manuscript "Predicting porcini: a decade of sporocarp monitoring reveals the meteorological triggers of Boletus edulis fruiting in central European beech forests"


library(data.table)
library(tidyverse)
library(zoo)
library(jtools)
library(scales)
library(viridis)
library(patchwork)
library(cowplot)


climate <- fread(file="Data_ready.csv") # Mnb = Mushroom number, Temp=Temperature, Precip=Precipitations
climate=climate |> mutate(rollP=abs(jitter(rollP, amount=0.05)), rollT=abs(jitter(rollT, amount=0.05))) # jittering values by 0.01, meaningless but improve kde compute

# The right align rolling temperature and precipitation 5 days averages rollT and rollPT were computed with the full annual using the function:
#climate <- climate |> mutate(rollT=rollmean(Temp, 5, align = "right", fill = "NA"),rollP=rollmean(Precip, 5, align = "right", fill="NA")) 
# then the data was subsetted for the first to last sporocarp sighting days.

# -------------------------------------------------------------- General -----------------------------------------------------------------#
# first # 08/05 and last date 11/21 
climate |>  group_by(year) |>  arrange(Date) |> group_split()
climate |>  group_by(year) |>  arrange(desc(Date)) |> group_split()
# total sporocarp nb per year
climate |> group_by(year) |> summarize(S=sum(Mnb))
# % of the data with 0 sporocarps
pct0 <- sum(climate$Mnb == 0) / length(climate$Mnb) * 100 ;pct0

#-------------------------------------------------------------- Figure 1  -----------------------------------------------------------------#

# fruiting season
rect_data <- climate  |> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  dplyr::filter(!is.na(Mnb)) |> 
  dplyr::group_by(year) |> arrange(date) |> filter(Mnb>0) |> 
  dplyr::summarise(xmin = min(Dm), xmax = max(Dm), .groups = "drop")

rect_data |> mutate(length=xmax-xmin+1) |> pull(length) |> summary()


# Get a symmetric y-range so both sides are comparable
max_precip <- max(climate$Precip, na.rm = TRUE)
max_M      <- max(climate$Mnb,    na.rm = TRUE)

P=climate |> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  ggplot(aes(x = Dm)) +
  geom_col(aes(y = Precip), fill = "#49a2f5ff") +
  geom_col(aes(y = -Mnb), fill = "black") +
  facet_wrap( ~ year, strip.position = "left", ncol = 1 ) +  
 scale_y_continuous(limits = c(-max_M, max_precip), breaks = c(0, 20, 40),  name   = "Precipitation (mm)",
   sec.axis = sec_axis(trans   = ~ ., breaks = -pretty(c(0, max_M)), labels  = function(x) abs(x), name= "Sporocarp number"))+
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
# adding season
P = P + 
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = -max_M, ymax = 0),
            inherit.aes = FALSE, fill = "grey60", alpha = 0.2)
gc()
max_temp <- max(climate$Temp, na.rm = TRUE)
max_M    <- max(climate$Mnb,na.rm = TRUE)

T= climate |> mutate(Dm=as.Date(format(Date,"%m-%d"),"%m-%d" )) |> 
  ggplot(aes(x = Dm)) +
  geom_col(aes(y = Temp), fill = "#de450f") +
  geom_col(aes(y = -Mnb), fill = "black") +
    facet_wrap( ~ year, strip.position = "left", ncol = 1) +  
 scale_y_continuous(limits = c(-max_M, max_temp), breaks = c(0, 30),  name   = "Temperature (C)",
    sec.axis = sec_axis(trans   = ~ ., breaks = -pretty(c(0, max_M)),labels  = function(x) abs(x), name = "Sporocarp number"))+
  theme_classic() + theme(legend.position = "none",strip.placement = "outside",strip.text.x.bottom = element_text(margin = margin(t = 6)),
    axis.text.x   = element_text(margin = margin(t = 4)),axis.title.x  = element_text(margin = margin(t = 8)), plot.margin   = margin(5, 5, 10, 5),
    axis.title.y.left = element_text(color = "#de450f", size=16),axis.text.y.left  = element_text(color = "#de450f"),
    axis.line.y.left = element_line(color = "#de450f"),
    axis.ticks.y.left = element_line(color="#de450f"), 
    axis.title.y.right = element_text(color = "black", size=16),axis.title.x.bottom = element_blank(),
    axis.text.y.right  = element_text(color = "black"), strip.text.y.left = element_text(size=16, color="black"), strip.background.y = element_blank())#;T
#adding season
T = T + 
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = -max_M, ymax = 0),
            inherit.aes = FALSE, fill = "grey60", alpha = 0.2)

gc()

# Density and boxplot of observed values during the study period every year
densT=climate |> ggplot()+
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


densP=climate |> ggplot()+
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


# Merging all
full_long3=plot_grid(T,NULL, densT,NULL, P,NULL, densP,
                      ncol=7, 
                      rel_widths = c(0.4,-0.001,0.1,0.01,0.4,0.001,0.1),
                      align = "v", 
                      axis = "b",
                      scale = 0.95,
                      labels = c("A", "", "", "","B", "" ), label_size = 26, label_fontface = "plain");full_long3

ggsave(full_long3, filename = "Figures/Figure1.pdf", width=18, height = 11)







#-------------------------------------------------------------- KDE + FIG 2  -----------------------------------------------------------------#
# remove no sporocarp values
expanded_data <- climate |>
  filter(Mnb > 0 & !is.na(Mnb))

h_rollP <- bw.nrd0(expanded_data$rollP)
h_rollT <- bw.nrd0(expanded_data$rollT)
h_rollP <- 2
h_rollT <- 2

# KDE
kde_data <- expanded_data |> 
    filter(!year == "2016") |> # too few data in 2016
    group_by(year) |> 
    reframe({kde <- MASS::kde2d(rollP, rollT, n = 60, h = c(h_rollP, h_rollT))
            tidyr::crossing(T = kde$y, P = kde$x) |> 
            mutate(density = as.vector(kde$z))}) |> 
    group_by(year) |> 
    mutate(scaled_density = scales::rescale(density, to = c(0, 1))) |> 
    ungroup()


study_wide <- expanded_data |> 
  reframe({kde <- MASS::kde2d(rollP, rollT, n = 300)
          tidyr::crossing(T = kde$y, P = kde$x) |> 
          mutate(density = as.vector(kde$z))}) |> 
  mutate(scaled_density = scales::rescale(density, to = c(0, 1)))

# A
p1=ggplot(study_wide, aes(x = P, y = T, fill = scaled_density)) +
  geom_raster() +
    theme_apa() +
  scale_fill_viridis_c(option="H",na.value = 0, begin = 0, end = 1, name="Density") + 
  labs(x = "Precipitation (mm)",y = "Temperature (C)", fill = "Density") +
   scale_x_continuous(limits = c(0, max(study_wide$P, na.rm = TRUE)), expand = c(0, 0)) +
  theme(panel.background = element_rect(fill="#30123b"),
         axis.text = element_text(color = "black", size = 16), 
         axis.ticks = element_line(linewidth = 1), axis.title.x.bottom = element_text(size=20, color="black"),
        axis.title.y.left = element_text(size=20, color="black"), legend.position="none");p1


# B
all_T_coords <- unique(kde_data$T)
all_P_coords <- unique(kde_data$P)

# for framing the data
bounds <- kde_data |> 
  group_by(year) |> 
  summarise(xmin = min(P, na.rm = TRUE),xmax = max(P, na.rm = TRUE),
            ymin = min(T, na.rm = TRUE),ymax = max(T, na.rm = TRUE),.groups = "drop" )


p <- ggplot(kde_data, aes(x = P, y = T, fill = scaled_density)) +
  geom_raster() +
  geom_rect(data = bounds,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,fill = NA,color = "white",  linewidth = 1 , linetype = "dashed") +
  theme_apa() +
  scale_fill_viridis_c(option = "H",na.value = "black", limits = c(0, 1), begin = 0.0, end = 0.95, name = "Density" ) +
  facet_wrap(~year) +
  xlab("Precipitation (mm)") +
  ylab("Temperature (C)") +
  scale_x_continuous(limits = c(-0.2, 21), expand = c(0, 0)) +
  theme(panel.spacing = unit(1, "cm", data = NULL),
    panel.background = element_rect(fill="#30123b"),
    axis.text = element_text(color = "black", size = 16),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x.bottom = element_text(size = 20, color = "black"),
    axis.title.y.left = element_text(size = 20, color = "black"),
    strip.text.x.top = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", size = 20, vjust=0.5,margin = margin(b = 20)),
  legend.text =element_text(color="black", size=16) );p

# legend on p1
p_for_legend <- p  

p_for_legend <- p + theme(legend.background = element_rect(fill = "white", color = "black"),legend.margin= margin(t = 8, r = 12, b = 12, l = 12),
    legend.title = element_text(margin = margin(b = 14)),
    legend.text  = element_text(margin = margin(t = 2, b = 2, l=5)),
    legend.spacing.y = unit(0.3, "lines"))

p_legend <- get_legend(p_for_legend)

p1_with_legend <- p1 + 
  inset_element(p_legend, left = 0.8, bottom = 0.08, right = 0.95, top = 0.25)

# p without legend
p2 <- p+theme(legend.position = "none")

full_fig2 =plot_grid(p1_with_legend, NULL,  p2 ,
                              rel_widths = c(0.5,0.02,0.5), 
                              ncol = 3, nrow = 1,
                              align = "h", axis="tb", 
                              labels = c("A", "","B" ), label_size = 26, label_fontface = "plain"); full_fig2

gc()
ggsave(full_fig2, filename = "Figures/Figure2.pdf",  width=18, height = 9.5)






#-------------------------------------------------------------- modelling part -----------------------------------------------------------------#
library(glmmTMB)
library(parameters)
library(ggeffects)

gc()
m_quad1b <- glmmTMB(Mnb ~ rollP + poly(rollT, 2) + (1 | year),data   = climate, family =  nbinom1())

model_parameters(m_quad1b)


# -------------------------------------Prediction -----------------------------------#
rT_seq <- seq(min(climate$rollT, na.rm = TRUE),max(climate$rollT, na.rm = TRUE),length.out = 1000)

newdat <- data.frame(rollT=rT_seq) |> mutate(rollP=mean(climate$rollP, na.rm = TRUE), year=2018) # picked a year but no importance for this prediction since we look at T and P effects

# Finding the peak/optimum
newdat$mu_hat <- predict(m_quad1b, newdata = newdat, type = "response",re.form = NA)   # exclude random effects
i_max      <- which.max(newdat$mu_hat)
peak_rollT <- newdat$rollT[i_max] #13.15 rounded as 13.2 in the ms


# plot
pred_comb <- ggpredict( m_quad1b, terms = c("rollT [all]", "rollP [0, 5, 10, 15,20]"))


gc()

pred_comb2 <- pred_comb |> as_tibble() |> 
  mutate( group_f = factor( group, levels = sort(unique(group), decreasing = TRUE)))

model=ggplot(pred_comb2, aes(x = x, y = predicted)) +
  geom_vline(xintercept =13.15, linetype="dashed" )+
  geom_ribbon(aes(x=x,ymin = predicted-std.error, ymax = predicted+ std.error, fill = group_f), 
alpha = 0.8, 
color = "transparent") +
  scale_colour_brewer(palette="Blues",direction = 1,name= "Precipitation (mm)") +
  scale_fill_brewer( palette= "Blues" ,direction = -1, name= "Precipitation (mm)") +
  ylab("Sporocarp number") +
  theme_apa() +
  xlab("Temperature (C)") +
  theme(axis.title.x.bottom = element_text(size = 16, color = "black"),
    legend.title        = element_text(size = 16, color = "black"),
    axis.title.y.left   = element_text(size = 16, color = "black"))+
  annotate(geom = "text", x=17, y=-0.5, label="Optimal temperature = 13.2 C");model

gc()

ggsave(model, filename = "Figures/Figure3.pdf", width=10, height = 8)











