library(sf)
library(tidyverse)
library(sdmTMB)
library(sfheaders)
library(ggpubr)
library(gridExtra)
library(grid)

toanalyze <- "gonad_weight"  # "weight" or "length" or "gonad_weight"

if (toanalyze=="length") model_save <- readRDS("output/Results_length.rds")
if (toanalyze=="weight") model_save <- readRDS("output/Results_weight.rds")
if (toanalyze=="gonad_weight")  model_save <- readRDS("output/Results_gonad_weight.rds")

#### Create residual QQplots
  
  plot_resid <- function(Period, Age_use){
    model_output <- model_save[[grep(paste0(toanalyze, "_", Period, "_age", Age_use),names(model_save))]]
    
    dat_use <- model_output$data
    source("R/residual_diagnostics.R")
    d <- dat_use
    
    d$residuals <- residuals.sdmTMB(model_output)
    d <- d[order(d$residuals),]
    d$qq <- (qnorm(ppoints(nrow(d))))
    
    p1 <- ggplot(d, aes(x= qq, y=residuals, col=length)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() + 
      scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample") + ggtitle(paste0(toanalyze, "_", Period, "_age ", Age_use)) +
    theme(axis.title = element_text(size=13),
            axis.text = element_text(size=10),
            axis.text.x = element_text(angle = 45, vjust=0.5),
            plot.title = element_text(hjust=0.5, size=13),
            strip.text = element_text(size=10),
            legend.text = element_text(size=11),
            legend.title = element_text(size=12),
            plot.margin=unit(c(1,1,0,0), "cm")) + ggtitle(paste0(Period, " age ", Age_use)) + labs(x=NULL)
    
    return(p1)
    
  }
  
  
  figs <- list()
  iter = 1
  
  # creating plots
  if (toanalyze %in% c("length", "weight")){
    for (Period in c("May",  "July", "Sept_Oct", "Nov_Dec")){
      for (Age_use in c(3:10)) {
        
        figs[[iter]] <- plot_resid(Period, Age_use)
        iter = iter + 1
      }
    }
  }
  if (toanalyze == "gonad_weight"){
    for (Period in c("Sept_Oct", "Nov_Dec")){
      for (Age_use in c(4:10)) {
        
        figs[[iter]] <- plot_resid(Period, Age_use)
        iter = iter + 1
      }
    }
  }
  
  
  # saving plots 
  if (toanalyze %in% c("length", "weight")){
    pp <- ggarrange(plotlist =figs[1:8], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "May.pdf"), width=250, height=180, dpi=300, unit="mm")
    pp <- ggarrange(plotlist =figs[9:16], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "July.pdf"), width=250, height=180, dpi=300, unit="mm")
    pp <- ggarrange(plotlist =figs[17:24], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "Sep_Oct.pdf"), width=250, height=180, dpi=300, unit="mm")
    pp <- ggarrange(plotlist =figs[25:31], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "Nov_Dec.pdf"), width=250, height=180, dpi=300, unit="mm")
  }
  if (toanalyze == "gonad_weight"){
    pp <- ggarrange(plotlist =figs[1:7], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "Sep_Oct.pdf"), width=250, height=180, dpi=300, unit="mm")
    pp <- ggarrange(plotlist =figs[8:14], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/QQplot_", toanalyze, "Nov_Dec.pdf"), width=250, height=180, dpi=300, unit="mm")
  }


#### Create residual QQplots
  
  plot_gear_resid <- function(Period, Age_use, xvar = "gear"){
    model_output <- model_save[[grep(paste0(toanalyze, "_", Period, "_age", Age_use),names(model_save))]]
    
    dat_use <- model_output$data
    #mutate for the gear code
    (fogear <- data.frame(country = rep('FO',8), rkode = sort(unique(data$rkode[data$country == 'FO'])), transl = c('peltrawl', 'demtrawl','peltrawl','peltrawl', 'purseseine', 'peltrawl','peltrawl',
                                                                                                                    'peltrawl')))
    
    (isgear <- data.frame(country = rep('IS',9), rkode = sort(unique(data$rkode[data$country == 'IS'])), transl = c('purseseine', NA,'peltrawl', NA, 'peltrawl', 'peltrawl', NA, NA, 'peltrawl')))
    
    (nogear <- data.frame(country = rep('NO',53), 
                          rkode = c(sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) == 31])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) == 32])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) == 34])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) == 35])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) == 37])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) %in% c(40,41)])),
                                    sort(unique(data$rkode[data$country == 'NO' & substr(data$rkode,1,2) %in% c(51,52)]))),
                          transl = c(rep('demtrawl', 3), rep('shrimptrawl', 5), rep('trawl',4), rep('peltrawl',23), rep('purseseine',7), rep('gillnet',9), rep('line',2))
    ))
    
    gear <- rbind(fogear, isgear, nogear)
    if (xvar == "gear") gear$gear <- ifelse(gear$transl %in% c("peltrawl"), "peltrawl", ifelse(gear$transl %in% c("purseseine"), "purseseine", "other"))
    if (xvar == "transl") gear$gear <- gear$transl
    
    dat_use <- dat_use %>% left_join(gear)
    if (xvar == "gear") dat_use$gear[which(is.na(dat_use$gear))] <- "other"
    
    source("R/residual_diagnostics.R")
    d <- dat_use
    
    d$residuals <- residuals.sdmTMB(model_output)
    d <- d[order(d$residuals),]
    d$qq <- (qnorm(ppoints(nrow(d))))
    
    dtable <- data.frame(table(d$gear, useNA="always"))
    colnames(dtable) <- c("gear", "n")
    dtable$y <- max(d$residuals)*1.2
    dtable$label <- paste0("n=", dtable$n)
    
    p1 <- ggplot(d, aes(x=gear, y=residuals)) + geom_boxplot() + theme_bw() + geom_hline(yintercept=0, linetype=2) +
      geom_text(data=dtable, aes(x=gear, y=y, label=label), size=2) + coord_cartesian(ylim=c(min(d$residuals),  max(d$residuals)*1.3)) +
      theme(axis.title = element_text(size=13),
            axis.text = element_text(size=10),
            axis.text.x = element_text(angle = 45, vjust=0.5),
            plot.title = element_text(hjust=0.5, size=13),
            strip.text = element_text(size=10),
            legend.text = element_text(size=11),
            legend.title = element_text(size=12),
            plot.margin=unit(c(1,1,0,0), "cm")) + ggtitle(paste0(Period, " age ", Age_use)) + labs(x=NULL)
    
    return(p1)
    
  }
  
  
  figs <- list()
  iter = 1
  
  # creating plots
  if (toanalyze == "weight"){
    for (Period in c("May",  "July", "Sept_Oct", "Nov_Dec")){
      for (Age_use in c(3:10)) {
        
        figs[[iter]] <- plot_gear_resid(Period, Age_use, xvar = "gear")
        iter = iter + 1
      }
    }
  }
  if (toanalyze == "gonad_weight"){
    for (Period in c("Sept_Oct", "Nov_Dec")){
      for (Age_use in c(4:10)) {
        
        figs[[iter]] <- plot_gear_resid(Period, Age_use, xvar = "gear")
        iter = iter + 1
      }
    }
  }
  
  # saving plots 
  if (toanalyze == "weight"){
    pp <- ggarrange(plotlist =figs[1:8], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "May.pdf"), width=180, height=180, dpi=400, unit="mm")
    pp <- ggarrange(plotlist =figs[9:16], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "July.pdf"), width=180, height=180, dpi=400, unit="mm")
    pp <- ggarrange(plotlist =figs[17:24], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "Sep_Oct.pdf"), width=180, height=180, dpi=400, unit="mm")
    pp <- ggarrange(plotlist =figs[25:31], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "Nov_Dec.pdf"), width=180, height=180, dpi=400, unit="mm")
  }
  if (toanalyze == "gonad_weight"){
    pp <- ggarrange(plotlist =figs[1:7], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "Sep_Oct.pdf"), width=180, height=180, dpi=400, unit="mm")
    pp <- ggarrange(plotlist =figs[8:14], ncol=3, nrow=3)
    ggsave(pp, file=paste0(getwd(), "/output/plots/resids/Gear_", toanalyze, "Nov_Dec.pdf"), width=180, height=180, dpi=400, unit="mm")
  }
  

