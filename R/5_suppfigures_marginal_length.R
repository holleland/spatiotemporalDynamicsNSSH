# ------------------
# Script for producing supplementary figures S3 - S4 
# marginal effects of length on weight and gonad weight 
# ------------------
rm(list=ls()) # clean slate

library(tidyverse)
length.models <- readRDS("output/Results_length.rds")
weight.models <- readRDS("output/Results_weight.rds")
gonad.models  <- readRDS("output/Results_gonad_weight.rds")


df <-rbind(
  tibble(
    type    = "Weight (g)",
    season  =  str_sub(names(weight.models),  8, 10),
    age     = str_sub(names(weight.models), -1, -1),
    data    = lapply(weight.models, function(x)x$marginal_length),
    formula = lapply(weight.models, function(x)x$formula)),
  tibble(
    type    = "Gonad weight (g)",
    season  = str_sub(names( gonad.models), 14, 16),
    age     = str_sub(names( gonad.models), -1, -1),
    data    = lapply(gonad.models, function(x)x$marginal_length),
    formula = lapply(gonad.models, function(x)x$formula))
) %>%
  mutate(type   = factor(type, levels = c("Length (cm)", "Weight (g)", "Gonad weight (g)")),
         season = factor(ifelse(season == "Jul", "July",
                                ifelse(season == "Nov", "Nov-Dec",
                                       ifelse(season == "Sep", "Sep-Oct", "May"))),
                         levels = c("May", "July", "Sep-Oct", "Nov-Dec")), 
         age    = factor(ifelse(age == 0, "10+",age), levels = c(3:9,"10+")))

df <- unnest(df, data)

# Transforming log-weights for seasons not being july: 
df[df$type == "Weight (g)" & df$season != "July", c("mean", "lwr","upr")]  <- 
  exp(df[df$type ==  "Weight (g)" & df$season != "July", c("mean", "lwr","upr")])

df[(df$type == "Gonad weight (g)" & df$season == "Nov-Dec") | 
     (df$type == "Weight (g)" & df$season != "July"),"length_use"] <-
  exp(df[(df$type == "Gonad weight (g)" & df$season == "Nov-Dec") | 
           (df$type == "Weight (g)" & df$season != "July"),"length_use"])

# ___ Figure S3 ____

ggplot(filter(df, type == "Weight (g)"), 
       aes(x = length_use, y = mean, color = season,
           fill = season, ymin = lwr, ymax = upr)) + 
  geom_ribbon(alpha = .4, color = NA) + 
  geom_line() + 
  facet_wrap( ~age) +
  scale_x_continuous(name = "Length (cm)",breaks = seq(15,45,10)) +
  scale_y_continuous(name = "Weight (g)", expand= c(0,0), 
                     breaks = seq(0,500,100), limits = c(0,590))+
  scale_color_viridis_d(name = "Season")+
  scale_fill_viridis_d( name = "Season")+
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
  theme_bw() + 
  theme(strip.background = element_rect(fill = "transparent", 
                                        color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))
ggsave("output/plots/S03_Simpler_with_common_scales.tiff", 
       width = 18, height=18, unit= "cm", dpi = 400)

# ___ Figure S4 ___

ggplot(filter(df, type == "Gonad weight (g)"), 
       aes(x = length_use, y = mean, color = season,
           fill = season, ymin = lwr, ymax = upr)) + 
  
  geom_ribbon(alpha = .4, color = NA) + 
  geom_line() + 
  facet_wrap( ~age) +
  scale_x_continuous(name = "Length (cm)", breaks = seq(15,40,5)) +
  scale_y_continuous(name = "Gonad weight (g)")+
  scale_color_viridis_d(begin = 2/3, name = "Season")+
  scale_fill_viridis_d(begin = 2/3,  name = "Season")+
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
  theme_bw() + 
  theme(strip.background = element_rect(fill = "transparent",
                                        color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0)) 
ggsave("output/plots/S04_Simpler_with_common_scales.tiff",
       width = 18, height=18, unit= "cm", dpi = 400)

