# -----------------------------------
# - Script for generating Figure 4: -
# -----------------------------------
rm(list=ls()) # clean slate

library(tidyverse)

length.models <- readRDS("output/Results_length.rds")
weight.models <- readRDS("output/Results_weight.rds")
gonad.models  <- readRDS("output/Results_gonad_weight.rds")


df <-rbind(
  tibble(
    type = "Length (cm)",
    season = str_sub(names(length.models),  8, 10),
    age = str_sub(names(length.models), -1, -1),
    data = lapply(length.models, function(x)x$marginal_year),
    formula = lapply(length.models, function(x)x$formula)),
tibble(
    type = "Weight (g)",
    season =  str_sub(names(weight.models),  8, 10),
    age = str_sub(names(weight.models), -1, -1),
    data = lapply(weight.models, function(x)x$marginal_year),
    formula = lapply(weight.models, function(x)x$formula)),
tibble(
    type = "Gonad weight (g)",
    season = str_sub(names( gonad.models), 14, 16),
    age    = str_sub(names( gonad.models), -1, -1),
    data = lapply(gonad.models, function(x)x$marginal_year),
    formula = lapply(gonad.models, function(x)x$formula))
) %>%
  mutate(type = factor(type, levels = c("Length (cm)", "Weight (g)", "Gonad weight (g)")),
         season = factor(ifelse(season == "Jul", "July",
                         ifelse(season == "Nov", "Nov-Dec",
                                ifelse(season == "Sep", "Sep-Oct", "May"))),
                         levels = c("May", "July", "Sep-Oct", "Nov-Dec"))) %>%
  unnest(data)

# Transforming log-weights for seasons not being july: 
df[df$type == "Weight (g)" & df$season != "July", c("mean", "lwr","upr")]  <- 
  exp(df[df$type ==  "Weight (g)" & df$season != "July", c("mean", "lwr","upr")])

# Remove negative values in confidence intervals: 
df <- df %>% mutate(lwr = pmax(0,lwr))

# For dealing with the discontinuities: 
# -----
fn <- function(data){
  y <- data$year
  idx <- c(1,diff(y));
  n <- length(y)
  i2 <- c(1,which(idx != 1), n+1)  
  return(cbind(data,"grp"=rep(1:length(diff(i2)), diff(i2))))
}

df2 <- df %>% arrange(year) %>% group_by(age, season, type) %>% nest() %>% 
  mutate(grp = map(.x=data, .f=fn)) %>%
  select(-data) %>%
  ungroup() %>% 
  mutate(row = as.numeric(row_number())) %>%
  unnest(grp) %>%
  mutate(grp = grp + 10*(row-1)) 
df <- right_join(df2,
                 df2 %>%  group_by(grp) %>% count() %>% ungroup())
# -------

df$age[df$age == 0] <- "10+"
df$age <- factor(df$age, levels = c(3:9,"10+"))

ggplot(df, aes(x = year, y = mean, ymin = lwr,ymax = upr, color = season, group = grp)) + 
    geom_vline(xintercept = 2005, linetype = "dashed",color = "red", lwd = .2)+
  geom_line()+#color = "#31688E") + 
  geom_ribbon(aes(fill = season),alpha = .4, color = NA) +
  geom_point(data= filter(df, n ==1), aes(color = season), size = .6)+
  geom_errorbar(data = filter(df, n ==1),aes(color = season))+
  facet_grid(type ~ age, scales = "free", switch = "y") + theme_bw() + 
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  scale_x_continuous(breaks = seq(1995, 2020, 10))+
  #scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linetype = "dashed", size = .3),
        strip.background = element_rect(fill = "transparent", color = "transparent"),
        strip.placement = "outside",
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10),
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        plot.margin = margin(c(10,0,0,0)))+
  guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0, alpha = 1)),
         color = "none") 
ggsave("output/plots/4_marginal_year_effects_by_age_and_season.tiff", width = 180, height = 150, unit = "mm",
       dpi = 400)
