library(tidyverse)
library(sf)


toanalyze <- "weight"  # "weight" or "gonad_weight"
#toanalyze <- "gonad_weight"
#toanalyze <- "length"

if (toanalyze=="weight") model_save <- readRDS("output/Results_weight.rds")
if (toanalyze=="gonad_weight")  model_save <- readRDS("output/Results_gonad_weight.rds")
if (toanalyze=="length") model_save <- readRDS("output/Results_length.rds")


# plot selected spatio-temporal effect 
plot_map_selected <- function(save=FALSE, percentile = 0.99){

  Norway <- st_crop(world, c(xmin = min(data$X)-15, ymin = min(data$Y)-10,
                             xmax = max(data$X)+15, ymax = max(data$Y)+10))
  
  # Now selecting only data points that do not fall on LANDS
  proj_UTM30N <- 32630 
  proj_UTM30N_km <- "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs"
  Norway_proj <- st_transform(Norway, crs=proj_UTM30N_km)
  
  to_work <- grep("age6", names(model_save))
  
  func1 <- function(x){
    years = sort(unique(model_save[[x]]$data$year))
    Nyear = length(years)
    qqq <- model_save[[x]][["epsilon_proj"]] 
    colnames(qqq) <- years
    qqq <- as.data.frame(qqq)
    qqq$LOC <- 1:nrow(qqq)
    test2 <- gather(qqq, key=year, value=Pred, 1:Nyear)
    test2$year <- as.numeric(as.character(test2$year))
    pred_grid <- model_save[[x]]$pred_grid
    pred_grid$LOC <- rep(1:(nrow(pred_grid)/Nyear), Nyear)
    test2 <- test2 %>% left_join(pred_grid)
    
    test2 <- test2 %>% select(LOC, year, Pred, X1000, Y1000, length_use, ID)
    test2 <- test2 %>% filter(year %in% c(1998,2003,2008,2013,2018))
    if (toanalyze == "weight") test2$Season <- substr(names(model_save)[x], 8,10)
    if (toanalyze == "length") test2$Season <- substr(names(model_save)[x], 8,10)
    if (toanalyze == "k") test2$Season <- substr(names(model_save)[x], 3,5)
    if (toanalyze == "gonad_weight") test2$Season <- substr(names(model_save)[x], 14,16)
    
    return(test2)
  }
  
  qqq <- do.call("rbind", lapply(to_work, function(x) func1(x)))
  qqq$Season <- as.factor(qqq$Season)
  if (toanalyze == "weight") qqq$Season <- factor(qqq$Season, labels=c("July","May","Nov-Dec","Sep-Oct"))
  if (toanalyze == "weight") qqq$Season <- factor(qqq$Season, levels=c("May","July","Sep-Oct","Nov-Dec"))
  if (toanalyze == "length") qqq$Season <- factor(qqq$Season, labels=c("July","May","Nov-Dec","Sep-Oct"))
  if (toanalyze == "length") qqq$Season <- factor(qqq$Season, levels=c("May","July","Sep-Oct","Nov-Dec"))
  if (toanalyze == "gonad_weight") qqq$Season <- factor(qqq$Season, labels=c("Nov-Dec","Sep-Oct"))
  if (toanalyze == "gonad_weight") qqq$Season <- factor(qqq$Season, levels=c("Sep-Oct","Nov-Dec"))
  
  pall <- ggplot(data=qqq) + facet_grid(Season~year) +

    geom_raster(data = filter(qqq, Season == "May") %>% mutate(May = ifelse(Pred>quantile(Pred, 1-(1-percentile)/2), 
                                                                            quantile(Pred, 1-(1-percentile)/2), 
                                                                            ifelse(Pred<quantile(Pred, (1-percentile)/2), 
                                                                                   quantile(Pred, (1-percentile)/2), Pred)) ),
                aes(x=X1000, y=Y1000, fill = May)) +
    geom_sf(data = Norway_proj$geometry, color ="grey27", size = .2)+
    scale_fill_gradient2(low="blue", mid="white", high="red", guide = guide_colorbar(order = 1)) +
    #guides(fill = guide_colorbar(order = 1))+
    ggnewscale::new_scale("fill")+
    geom_raster(data = filter(qqq, Season == "July") %>% mutate(July = ifelse(Pred>quantile(Pred, 1-(1-percentile)/2), 
                                                                            quantile(Pred, 1-(1-percentile)/2), 
                                                                            ifelse(Pred<quantile(Pred, (1-percentile)/2), 
                                                                                   quantile(Pred, (1-percentile)/2), Pred)) ),
                aes(x=X1000, y=Y1000, fill = July)) +
    geom_sf(data = Norway_proj$geometry, color ="grey27", size = .2)+
    scale_fill_gradient2(low="blue", mid="white", high="red",
                          guide = guide_colorbar(order = 2)) +
    ggnewscale::new_scale("fill")+
   # guides(fill = guide_colorbar(order = 2))+
    geom_raster(data = filter(qqq, Season == "Sep-Oct") %>% mutate("Sep_Oct" = ifelse(Pred>quantile(Pred, 1-(1-percentile)/2), 
                                                                              quantile(Pred, 1-(1-percentile)/2), 
                                                                              ifelse(Pred<quantile(Pred, (1-percentile)/2), 
                                                                                     quantile(Pred, (1-percentile)/2), Pred)) ),
                aes_string(x="X1000", y="Y1000", fill = "Sep_Oct")) +
    geom_sf(data = Norway_proj$geometry, color ="grey27", size = .2)+
    scale_fill_gradient2(low="blue", mid="white", high="red", name = "Sep-Oct",
                         guide = guide_colorbar(order = 3)) +
   # guides(fill = guide_colorbar(order = 3))+
    ggnewscale::new_scale("fill")+
    
    geom_raster(data = filter(qqq, Season == "Nov-Dec") %>% mutate("Nov_Dec" = ifelse(Pred>quantile(Pred, 1-(1-percentile)/2), 
                                                                                      quantile(Pred, 1-(1-percentile)/2), 
                                                                                      ifelse(Pred<quantile(Pred, (1-percentile)/2), 
                                                                                             quantile(Pred, (1-percentile)/2), Pred)) ),
                aes_string(x="X1000", y="Y1000", fill = "Nov_Dec")) +
    geom_sf(data = Norway_proj$geometry, color ="grey27", size = .2)+
    scale_fill_gradient2(low="blue", mid="white", high="red", name = "Nov-Dec",guide = guide_colorbar(order = 4)) +
   # guides(fill = guide_colorbar(order = 4))+
    xlab("Longitude") + ylab("Latitude") +
    theme_bw() +
    ggtitle("Estimated spatio-temporal random effect - age 6") +

    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          plot.title = element_text(hjust=0.5, size=15),
          strip.text = element_text(size=10)) + 
    coord_sf(xlim=range(qqq$X1000), ylim=range(qqq$Y1000))
  if(toanalyze != "gonad_weight"){
  ggsave(pall,#+           guides(fill = guide_colorbar(barheight = 25/2.5)),
         file = paste0(getwd(), "/plots/epsilon_selected_",toanalyze,"_", round(100*percentile,0),  ".tiff"),
         width = 10, height = 8, dpi= 400, device = "tiff")
  }else{
    ggsave(pall,#+guides(fill = guide_colorbar(barheight = 12.5/2.5)), 
           file = paste0(getwd(), "/plots/epsilon_selected_",toanalyze, "_", round(100*percentile,0), ".tiff"),
           width = 10, height = 5, dpi= 400, device = "tiff")
  }
  
}

plot_map_selected(percentile = 1)
  