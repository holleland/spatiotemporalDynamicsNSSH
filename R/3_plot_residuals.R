library(sf)
library(tidyverse)
library(sdmTMB)
library(sfheaders)
library(ggpubr)
library(gridExtra)
library(grid)

toanalyze <- "weight"  # "weight" or "length" or "gonad_weight"

if (toanalyze=="length") model_save <- readRDS("output/Results_length.rds")
if (toanalyze=="weight") model_save <- readRDS("output/Results_weight.rds")
if (toanalyze=="gonad_weight")  model_save <- readRDS("output/Results_gonad_weight.rds")

####  Do some summary plotting 

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



