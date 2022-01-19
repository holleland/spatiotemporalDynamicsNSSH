##########################################
###  Main code for running the model
##########################################

### source data
source("R/0_data.R")

### load libraries
library(sf)
library(tidyverse)
library(sdmTMB)
library(sfheaders)
library(ggpubr)
library(gridExtra)
library(grid)

### read in shapefiles for mapping purpose
world <- st_read("./shapefile/ne_10m_land.shp")
Atlantic <- st_crop(world, c(xmin = -65, ymin = 46, xmax = 32, ymax = 80))

### Get ready to run the model 

## Empty list for saving all model results
model_save <- list()

## Specifying some general model setting
toanalyze <- "weight"   # which analysis to perform: one of "weight" or "length" or "gonad_weight" 
doresidualplot = FALSE  # create residual plots? 


if (toanalyze %in% c("length", "weight")) {
  seasons = c("May",  "July", "Sept_Oct", "Nov_Dec")
  ages = 3:10
}
if (toanalyze %in% c("gonad_weight")) {
  seasons = c("Sept_Oct", "Nov_Dec")
  ages = 4:10
}

for (Period in seasons){
  for (Age_use in ages) {
 
	#-----------------------------------------------------------------------------------------------#
  ## Stage 1: data preparation - Selecting the data to keep in the analysis
	#-----------------------------------------------------------------------------------------------#
  # converting all ages 10>= to 10+ group
  dat_use <- data %>% mutate(age= ifelse(age >= 10, 10, age))
  # Selecting only data from Norway, Faroes, and Iceland
  dat_use <- dat_use %>% 
  filter(country %in% c("NO", "FO", "IS")) %>%
  filter(!is.na(k)) %>%  
  mutate(length_std = (length - mean(length))/sd(length)) %>%
  group_by(year, station) %>% mutate(n=n()) %>% filter(n>=1)
  nrow(dat_use)
  table(dat_use$country)
  table(dat_use$month)
  dat_use <- filter(dat_use, 
                    !is.na(k) & !is.na(month)& !is.na(year) & !is.na(sex) & !is.na(X) & age==Age_use)
  
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
  
  dat_use <- dat_use %>% left_join(gear)

                                
  ## more specific filtering depending on the type of analysis and season                            
  if (toanalyze=="gonad_weight") dat_use <- filter(dat_use, !is.na(gonadweight))
  if (Period == "May") dat_use <- dat_use %>% filter(X<22, month %in% c(5))
  if (Period == "July") dat_use <- dat_use %>% filter(X<22, month %in% c(7))
  if (Period == "Sept_Oct") dat_use <- dat_use %>% filter(X<22, month %in% c(9,10), maturation <= 5, maturation >=3, !is.na(maturation))
  if (Period == "Nov_Dec") dat_use <- dat_use %>% filter(month %in% c(11,12), maturation <= 5, maturation >=3, !is.na(maturation))
  
  if (toanalyze=="gonad_weight") 
  {
    # deleting gonadweights of herring with maturity stages 4 and 5 ## rationale: the GSI needs to be at least 2.5 % of body weight in these maturity stages (the 95% percentile of fish in maturity stage 8). This seems to delete most outlier-points of maturity stages 4 and 5 
    dat_use$gonadweight[which((dat_use$maturation %in% 4:5 & dat_use$gsi <= 0.025 | (dat_use$country == 'FO' & dat_use$station == '20145029') == TRUE))] <- NA
    
    # deleting gonadweights of herring with maturity stage 3 ## rationale: the GSI needs to be at least 1.25 % of body weight in these maturity stages (1% arbitrarily chosen by fiddling back and forth - 158 gonad weights deleted with 1.25% - 107 deletions with 1%). This seems to delete most outlier-points of maturity stage 3 for older fish - but we loose some correct gonad weights for younger and smaller fish ...
    dat_use$gonadweight[which( (dat_use$maturation %in% 3 & dat_use$gsi <= 0.01 & dat_use$age > 4)== TRUE)] <- NA

    # Now filtering out these values
    dat_use <- dat_use %>% filter(!is.na(gonadweight))
  }
  
  #- load shapefiles -
  Norway <- st_crop(world, c(xmin = min(data$X)-2, ymin = min(data$Y)-2,
                              xmax = max(data$X)+2, ymax = max(data$Y)+2))
  
	# Now changing projection to UTM
  proj_UTM30N <- 32630 
  proj_UTM30N_km <- "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs"
  Norway_proj <- st_transform(Norway, crs=proj_UTM30N_km)
  survey <- dat_use %>% dplyr::select(X, Y, k, year) %>% 
    st_as_sf(crs = 4326, coords = c("X", "Y"))  %>%
    st_transform(proj_UTM30N_km)
  surv_utm_coords <- st_coordinates(survey)
   
  # Then we will scale coordinates to km so the range parameter is on a reasonable scale for estimation:
  dat_use$X1000 <- surv_utm_coords[,1]
  dat_use$Y1000 <- surv_utm_coords[,2]
  
  # change to sf object
  dat_use_sf <- st_as_sf(dat_use, coords = c("X1000", "Y1000"), 
                         crs = proj_UTM30N_km)
  
  # I am also projecting in lonlat to get the respective columns in the pred_grid data frame
  dat_use_lonlat <- st_transform(dat_use_sf,crs = 4326)

  #revert to data frame
  dat_use <- dat_use %>% 
    add_column(LON=st_coordinates(dat_use_lonlat)[,1],
               LAT=st_coordinates(dat_use_lonlat)[,2]) 
  rm(dat_use_lonlat, dat_use_sf,surv_utm_coords)
  
  # Some variable transformation 
  dat_use <- dat_use %>% mutate(log_weight = log(weight),
                                                          log_length = log(length),
                                                          log_gonad = log(gonadweight+0.1)) 
  
  # creating a new regional variable
  dat_use <- dat_use %>% mutate(Area = ifelse((LON > 10 & LAT > 67 & LAT < 68), "Lofoten", ifelse((LON < -10 & LAT < 68), "Iceland", "Nowegian")))
  
  if (toanalyze == "length") {
    dat_use$resp = dat_use$length
    dat_use$length_use = dat_use$length
    fixed_effect_formula <- formula(resp ~  as.factor(year))
  }
  
  if (toanalyze == "weight") {
    if (Period != "July") dat_use$resp = dat_use$log_weight
    if (Period == "July") dat_use$resp = dat_use$weight
    if (Period != "July") dat_use$length_use = dat_use$log_length
    if (Period == "July") dat_use$length_use = dat_use$length
    dat_use$Y1000_scaled = as.numeric(scale(dat_use$Y1000))
    if (Period != "July") fixed_effect_formula <- formula(resp ~  as.factor(year) + s(length_use, k=3))
    if (Period == "July") fixed_effect_formula <- formula(resp ~  as.factor(year) + s(length_use, k=5))# + Y1000_scaled)
  }
  
  if (toanalyze == "gonad_weight") {
    dat_use$resp = (dat_use$gonadweight)
    if (Period == "Sept_Oct") dat_use$length_use = (dat_use$length)
    if (Period != "Sept_Oct") dat_use$length_use = log(dat_use$length)
    dat_use$Y1000_scaled = as.numeric(scale(dat_use$Y1000))
    fixed_effect_formula <- formula(resp ~  as.factor(year) + as.factor(sex) + s(length_use, k=3) + s(k, k=3) + s(yday, k=3)) # + Y1000_scaled)
  }
  
	# some plotting
  ggplot(dat_use, aes(x= length_use, y=resp)) + geom_point() + geom_smooth(method="gam", formula = y ~ s(x, k=3))
  ggplot(dat_use, aes(x= yday, y=resp)) + geom_point() + geom_smooth(method="gam", formula = y ~ s(x, k=3))
  ggplot(dat_use, aes(x= k, y=resp)) + geom_point() + geom_smooth(method="gam", formula = y ~ s(x, k=3))
 
 
	#-----------------------------------------------------------------------------------------------#
	## Stage 2: preparing the model 
	#-----------------------------------------------------------------------------------------------#
	
	# Setting up the design matrix
  mf <- model.frame(mgcv::interpret.gam(fixed_effect_formula)$fake.formula, dat_use)
  mgcv_mod <- mgcv::gam(fixed_effect_formula, data = dat_use)
  X <- model.matrix(mgcv_mod)
  yobs <- model.response(mf, "numeric")
  Nobs = length(yobs)
  
  # For the mixture proportion (if using a mixture model)
  mixture_formula <-  formula(weight ~  s(resp, k=3))# + Y1000_scaled)
  mgcv_mixture <- mgcv::gam(mixture_formula, data = dat_use)
  X_mix <- model.matrix(mgcv_mixture)
  
  # Construct the base mesh structure: this is needed because using the spde approach as used for INLA spatial models
  if (nrow(dat_use) > 650) spde <- make_mesh(dat_use, xy_cols = c("X1000", "Y1000"), n_knots = 200, type = "kmeans")
  if (nrow(dat_use) < 650) spde <- make_mesh(dat_use, xy_cols = c("X1000", "Y1000"), n_knots = 150, type = "kmeans")
  
  Nyear <- length(unique(dat_use$year))

  # Now getting ready to create the barrier mesh 
	# first defining the barrier area (land + other fish zone delimitation) 
  South_end <- data.frame(lon=c(16,45,45,16,16), lat=c(67,68,55,55,67))
  South_end$ID <- "South_end"
  South_end_sf <- st_as_sf(x=South_end, coords = c("lon", "lat")) %>% st_set_crs(4326) %>%
    summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") %>%
    st_transform(crs = proj_UTM30N_km)
  Norunion <- st_union(Norway_proj, South_end_sf)
  # ggplot(Norunion) + geom_sf() + theme_bw() #+ geom_sf(data = South_end_sf, size = 0.2, fill="lightblue1")

  # then incorporating this barrier in the mesh structure
    Barrier_range_fraction = 0.1  
    bspde <- add_barrier_mesh(
      spde, Norunion, range_fraction = Barrier_range_fraction,
      proj_scaling = 1, plot = TRUE)
    
  # preparing the data to be used in the spatial model i.e. with an appropriate reference to the above mesh structure   
  dat_use$sdm_orig_id <- seq(1, nrow(dat_use))
  dat_use$sdm_x <- spde$loc_xy[,1,drop=TRUE]
  dat_use$sdm_y <- spde$loc_xy[,2,drop=TRUE]
  fake_data <- unique(data.frame(sdm_x = dat_use$sdm_x, sdm_y = dat_use$sdm_y))
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))
  dat_use <- base::merge(dat_use, fake_data, by = c("sdm_x", "sdm_y"),
                      all.x = TRUE, all.y = FALSE)
  dat_use <- dat_use[order(dat_use$sdm_orig_id),, drop=FALSE]
  
	# creating the A matrix as in INLA to link observation points to the mesh locations
  A <- spde$A
  A_st <- INLA::inla.spde.make.A(spde$mesh,
                                 loc = as.matrix(fake_data[, c("sdm_x", "sdm_y"), drop = FALSE]))
     
  n_s <- nrow(spde$mesh$loc)
  Nmesh <- n_s
  
  # Now prepare the input data for the prediction (based on 10 x 10 km grid)  
    pred_grid <- expand.grid(X1000 = seq(min(dat_use$X1000),
                                         max(dat_use$X1000), by = 10),
                             Y1000 = seq(min(dat_use$Y1000),
                                         max(dat_use$Y1000), by = 10),
                             year = sort(unique(dat_use$year)),
                             length_use = median(dat_use$length_use),
                             k = median(dat_use$k),
                             Y1000_scaled = 0,
                             rkode = "3710",
                             yday=median(dat_use$yday),
                             Area = "Norwegian", 
                             sex=2)

  # Convert to sf object to intersect the pred_grid with the convex hull.
  # The pred_grid is in km and this needs to be indicated in the crs when 
  # creating the sf object. It corresponds to "+proj=utm +zone=30 +ellps=WGS84 +units=km +no_defs". 
  pred_grid_sf <- st_as_sf(pred_grid, coords = c("X1000", "Y1000"), crs = proj_UTM30N_km)
  # get the convex hull
  hull_points <- st_as_sf(dat_use, coords = c('X1000','Y1000'), crs = proj_UTM30N_km)
  hull <- hull_points %>%
    summarise( geometry = st_combine( geometry ) ) %>%
    st_convex_hull()
 
  # cropping the prediction area
  pred_grid_croparea <- st_as_sf(subset(pred_grid, year==2019), coords = c('X1000','Y1000'), crs = proj_UTM30N_km) 
  pred_grid_croparea_sf <- st_intersection(pred_grid_croparea,hull)
  pred_grid_croparea <- sfheaders::sf_to_df(pred_grid_croparea_sf)
  pred_grid_croparea_lonlat <- st_transform(pred_grid_croparea_sf,crs = 4326)
  pred_grid_crop <- data.frame(X1000=pred_grid_croparea$x,
                               Y1000=pred_grid_croparea$y, 
                               LON=st_coordinates(pred_grid_croparea_lonlat)[,1],
                               LAT=st_coordinates(pred_grid_croparea_lonlat)[,2])
  pred_grid_crop$XY <- apply(pred_grid_crop[,c('X1000','Y1000')], 1, function(x) paste(x, collapse="_"))
  pred_grid$XY <- apply(pred_grid[,c('X1000','Y1000')], 1, function(x) paste(x, collapse="_"))
  pred_grid_crop_all <- pred_grid %>% filter(XY %in% unique(pred_grid_crop$XY))
  pred_grid_crop_all <- pred_grid_crop_all[,-grep("XY",colnames(pred_grid_crop_all))]
  
  pred_grid <- pred_grid_crop_all
  pred_grid$ID <- 1:nrow(pred_grid)

  proj_data <- subset(pred_grid, year==2012)
  
	# creating the design matrix as well as the A matrix for the projection grid
  X_proj <- mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = pred_grid)

  A_proj<- INLA::inla.spde.make.A(bspde$mesh,
                                  loc = as.matrix(proj_data[, c("X1000", "Y1000"), drop = FALSE]))
  
  # Npred
  Npred = nrow(A_proj)
  
  ### Some functions borrowed from the R package sdmTMB (https://github.com/pbs-assess/sdmTMB)
  make_barrier_spde <- function(spde) {
    if ("spde_barrier" %in% names(spde)) {
      C0 <- spde$spde_barrier$C[[1]]
      C1 <- spde$spde_barrier$C[[2]]
      D0 <- spde$spde_barrier$D[[1]]
      D1 <- spde$spde_barrier$D[[2]]
      .I <- spde$spde_barrier$I
    } else {
      C0 <- rep(1, 2)
      C1 <- rep(1, 2)
      D0 <- Matrix::Matrix(0, 1, 1)
      D1 <- Matrix::Matrix(0, 1, 1)
      .I <- Matrix::Matrix(0, 1, 1)
    }
    list(C0 = C0, C1 = C1, D0 = D0, D1 = D1, I = .I)
  }  
    
  make_year_i <- function(x) {
    x <- as.integer(as.factor(x))
    x - min(x)
  }    
  
  set_par_value <- function(opt, par) {
    as.numeric(opt$par[par == names(opt$par)])
  }   
  
  
	# Final model specification/configurations obtained from iterative model runs
		keep_omega = FALSE    # do not model the average spatial field
    family_choice = 2                                                                  # student t in most cases (see below for meaning)
    family_link = 1                                                                    # student t in most cases
    if (Period == "May") dfs = 6                                                       # for age 5, 6
    if (Period == "July") dfs = 6                                                      # for age 5, 6
    if (Period == "Sept_Oct") dfs = 20                                                 # for age 5, 6
    if (Period == "Sept_Oct" & toanalyze == "gonad_weight") { 
      family_choice = 2
      family_link = 1
      dfs = 6
      if (Age_use %in% c(5,9)) dfs = 5
    }
    if (Period == "Nov_Dec" & toanalyze == "gonad_weight") { 
      family_choice = 2
      family_link = 1
      dfs = 6
      if (Age_use %in% c(5,6,7,9)) dfs = 5
    }
    if (Period == "July" & toanalyze == "weight") { 
      family_choice = 2
      family_link = 1
    }

    # data list needed for the TMB model 
    tmb_data <- list(
      X                = as.matrix(X),
      yobs             = as.numeric(yobs),
      Nobs             = Nobs,
      spde_barrier     = make_barrier_spde(bspde),
      spde_barrier1994 = make_barrier_spde(bspde),
      spde_barrier1995 = make_barrier_spde(bspde),
      spde_barrier1996 = make_barrier_spde(bspde),
      spde_barrier1997 = make_barrier_spde(bspde),
      spde_barrier1998 = make_barrier_spde(bspde),
      spde_barrier1999 = make_barrier_spde(bspde),
      spde_barrier2000 = make_barrier_spde(bspde),
      spde_barrier2001 = make_barrier_spde(bspde),
      spde_barrier2002 = make_barrier_spde(bspde),
      spde_barrier2003 = make_barrier_spde(bspde),
      spde_barrier2004 = make_barrier_spde(bspde),
      spde_barrier2005 = make_barrier_spde(bspde),
      spde_barrier2006 = make_barrier_spde(bspde),
      spde_barrier2007 = make_barrier_spde(bspde),
      spde_barrier2008 = make_barrier_spde(bspde),
      spde_barrier2009 = make_barrier_spde(bspde),
      spde_barrier2010 = make_barrier_spde(bspde),
      spde_barrier2011 = make_barrier_spde(bspde),
      spde_barrier2012 = make_barrier_spde(bspde),
      spde_barrier2013 = make_barrier_spde(bspde),
      spde_barrier2014 = make_barrier_spde(bspde),
      spde_barrier2015 = make_barrier_spde(bspde),
      spde_barrier2016 = make_barrier_spde(bspde),
      spde_barrier2017 = make_barrier_spde(bspde),
      spde_barrier2018 = make_barrier_spde(bspde),
      spde_barrier2019 = make_barrier_spde(bspde),
      barrier_scaling  = bspde$barrier_scaling,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = dat_use$sdm_spatial_id - 1L,  
      year_i           = make_year_i(dat_use$year),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = Npred,
      do_predict       = 1,
      calc_se					 = 0,
      X_proj					 = X_proj,
      A_proj					 = A_proj,
      family           = family_choice,         # 0: tweedie, 1: gamma, 2: student, 3:lognormal, 4:gaussian, 5: skewed-normal
      link             = family_link,						# 0: identity, 1: log, 2: logit, 3: inverse
      AR1_spatial      = 0,
      df               = dfs,          					# df for the student t distribution
      exclude_RE_pred  = 0             					# exclude RE from prediction. Useful if wanting marginal effects
    )
  
    # parameter list needed for the TMB model 
    tmb_params <- list(
      beta         = c(ifelse(family_link==1, log(mean(dat_use$resp)), mean(dat_use$resp)),rep(0, ncol(X)-1)),
      omega        = rep(0, Nmesh),
      epsilon_st   = matrix(0, nrow=Nmesh, ncol=Nyear),
      transf_rho   = 0.0,
      logKappa     = 0.0,
      logTauO      = 0.0,
      logTauE      = 0.1,
      thetaf       = 0.1, 
      ln_phi       = -1    
    )  

		# if model does not converge with initial parameters
    # tmb_params$logKappa = -3
    # tmb_params$logTauE = 3
    # tmb_params$ln_phi = 3
 
 
	#-----------------------------------------------------------------------------------------------#
	## Stage 3: Running the TMB model 
	#-----------------------------------------------------------------------------------------------#
    # loading required TMB library 
    library(TMB)
    library(TMBhelper)
    
		# specifying path to the TMB model template file, compile and load it
    version <- paste0(getwd(), "/src/herring_sdm")
    #dyn.unload(dynlib(version))
    compile(paste0(version, ".cpp"), "-O1 -g", DLLFLAGS="")
    dyn.load(dynlib(version))
          
		# Specifying which parameters to fix and what not			
    Map_phase2 <- list()
    if (tmb_data$AR1_spatial == 0) Map_phase2$transf_rho <- factor(NA)

    if (keep_omega == FALSE) {
      Map_phase2$omega <- factor(rep(NA, length(tmb_params$omega)))
      Map_phase2$logTauO <- factor(NA)
      if(tmb_data$family != 0) Map_phase2$thetaf <- factor(NA)
			
      tmb_obj_phase2 <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params, map = Map_phase2,
                                       random = c("epsilon_st"), DLL = "herring_sdm", silent = FALSE)
    }
    if (keep_omega == TRUE) {
      if(tmb_data$family != 0) Map_phase2$thetaf <- factor(NA)
      tmb_obj_phase2 <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params, map = Map_phase2,
                                       random = c("omega", "epsilon_st"), DLL = "herring_sdm", silent = FALSE)
    }
    
		# running the TMB model 
    opt_phase2 <- fit_tmb(tmb_obj_phase2, lower=-15, upper=15, loopnum = 3, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          

    opt_phase2$diagnostics[which(abs(opt_phase2$diagnostics$MLE) > 14),]


	#-----------------------------------------------------------------------------------------------#
	## Stage 4: Saving model outputs
	#-----------------------------------------------------------------------------------------------#
    
    years = unique(dat_use$year)
    
    model_output <- list()
    model_output$formula <- fixed_effect_formula
    model_output$opt <- opt_phase2
    model_output$obj <- tmb_obj_phase2
    model_output$tmb_data <- tmb_data
    model_output$Period <- Period
    model_output$Age_use <- Age_use
    model_output$data <- dat_use
    model_output$response <- yobs
    model_output$pred_grid <- pred_grid
    model_output$mu <- tmb_obj_phase2$report()$mu
    model_output$mu_proj <- tmb_obj_phase2$report()$mu_proj
    model_output$epsilon_proj <- tmb_obj_phase2$report()$epsilon_st_A_proj
    model_output$omega_proj <- tmb_obj_phase2$report()$omega_A_proj
    model_output$family <- ifelse(tmb_data$family==0, "tweedie", ifelse(tmb_data$family==1, "gamma", ifelse(tmb_data$family==2, "student", ifelse(tmb_data$family==3, "lognormal", ifelse(tmb_data$family==4, "gaussian", "skewednormal")))))
    
    # Perform some residual diagnostics
    
    source("R/residual_diagnostics.R")
    d <- dat_use
    
    d$residuals <- residuals.sdmTMB(model_output)
    d <- d[order(d$residuals),]
    d$qq <- (qnorm(ppoints(nrow(d))))
    p1 <- ggplot(d, aes(x= qq, y=residuals, col=length)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() + 
      scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample") + ggtitle(paste0(toanalyze, "_", Period, "_age ", Age_use))
    p1
    ggsave(p1, file = paste0(getwd(), "/output/plots/", toanalyze, "_period", Period, "age", Age_use, "_resid_sdmTMB_noomega.pdf"),
           width = 8, height =7, dpi= "retina", device = "pdf")
    
    p2 <- ggplot(d, aes(x=residuals)) + geom_histogram() + theme_bw()
    p2
    ggsave(p2, file = paste0(getwd(), "/output/plots/", toanalyze, "_period", Period, "age", Age_use, "_resid_hist_sdmTMB_noomega.pdf"),
           width = 8, height =7, dpi= "retina", device = "pdf")

		# Now create a new prediction data frame for evaluating the marginal effects
    # For the marginal year effect
    pred_grid_year <- expand.grid(X1000 = c(median(dat_use$X1000),mean(dat_use$X1000)),
                                  Y1000 = c(median(dat_use$Y1000),mean(dat_use$Y1000)),
                                  year = sort(unique(dat_use$year)),
                                  length_use = median(dat_use$length_use),
                                  k = median(dat_use$k),
                                  Y1000_scaled = 0,
                                  yday=median(dat_use$yday),
                                  Area = "Norwegian", 
                                  sex=2)

    X_proj_year <- mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = pred_grid_year)

    A_proj_year <- INLA::inla.spde.make.A(bspde$mesh,
                                    loc = as.matrix(pred_grid_year[which(pred_grid_year$year == 2019), c("X1000", "Y1000"), drop = FALSE]))
    tmb_data_year <- list(
      X                = as.matrix(X),
      yobs             = as.numeric(yobs),
      Nobs             = Nobs,
      spde_barrier     = make_barrier_spde(bspde),
      spde_barrier1994 = make_barrier_spde(bspde),
      spde_barrier1995 = make_barrier_spde(bspde),
      spde_barrier1996 = make_barrier_spde(bspde),
      spde_barrier1997 = make_barrier_spde(bspde),
      spde_barrier1998 = make_barrier_spde(bspde),
      spde_barrier1999 = make_barrier_spde(bspde),
      spde_barrier2000 = make_barrier_spde(bspde),
      spde_barrier2001 = make_barrier_spde(bspde),
      spde_barrier2002 = make_barrier_spde(bspde),
      spde_barrier2003 = make_barrier_spde(bspde),
      spde_barrier2004 = make_barrier_spde(bspde),
      spde_barrier2005 = make_barrier_spde(bspde),
      spde_barrier2006 = make_barrier_spde(bspde),
      spde_barrier2007 = make_barrier_spde(bspde),
      spde_barrier2008 = make_barrier_spde(bspde),
      spde_barrier2009 = make_barrier_spde(bspde),
      spde_barrier2010 = make_barrier_spde(bspde),
      spde_barrier2011 = make_barrier_spde(bspde),
      spde_barrier2012 = make_barrier_spde(bspde),
      spde_barrier2013 = make_barrier_spde(bspde),
      spde_barrier2014 = make_barrier_spde(bspde),
      spde_barrier2015 = make_barrier_spde(bspde),
      spde_barrier2016 = make_barrier_spde(bspde),
      spde_barrier2017 = make_barrier_spde(bspde),
      spde_barrier2018 = make_barrier_spde(bspde),
      spde_barrier2019 = make_barrier_spde(bspde),
      barrier_scaling  = bspde$barrier_scaling,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = dat_use$sdm_spatial_id - 1L,
      year_i           = make_year_i(dat_use$year),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = nrow(A_proj_year),
      do_predict       = 1,
      calc_se					 = 1,
      X_proj					 = X_proj_year,
      A_proj					 = A_proj_year,
      family           = family_choice,            # 0: tweedie, 1: gamma, 2: student, 3:lognormal
      link             = family_link,						 # 0: identity, 1: log, 2: logit, 3: inverse
      AR1_spatial      = 0,
      df               = dfs,          # df for the student t distribution
      exclude_RE_pred  = 1             # exclude RE from prediction. Useful if wanting marginal effects
    )

    Map_phase2 <- list()
    if (tmb_data_year$AR1_spatial == 0) Map_phase2$transf_rho <- factor(NA)

       Map_phase2$omega <- factor(rep(NA, length(tmb_params$omega)))
      Map_phase2$logTauO <- factor(NA)
      if(tmb_data_year$family != 0) Map_phase2$thetaf <- factor(NA)

      tmb_obj_phase3 <- TMB::MakeADFun(data = tmb_data_year, parameters = tmb_params, map = Map_phase2,
                                       random = c("epsilon_st"), DLL = "herring_sdm", silent = FALSE)

    par_vals <- (tmb_obj_phase2$env$last.par.best)
    par_vals[grep("epsilon_st", names(par_vals))] <- 0
    marginal_year <- tmb_obj_phase3$report(par_vals)$mu_proj

    qwe <- sdreport(tmb_obj_phase3, par.fixed = opt_phase2$par)
    qqq <- summary(qwe, "report")[(0:(length(unique(dat_use$year))-1)*4)+1,]
    colnames(qqq) <- c("mean", "se")
    qqq <- as.data.frame(qqq)
    qqq <- qqq %>% mutate(lwr=mean - 1.96*se,
                          upr=mean + 1.96*se)
    qqq$year = sort(unique(dat_use$year))
 
    if (toanalyze != "length"){
    ## For the marginal length effect
    pred_grid_length <- expand.grid(X1000 = c(median(dat_use$X1000)),
                             Y1000 = c(median(dat_use$Y1000)),
                             length_use = seq(min(dat_use$length_use), max(dat_use$length_use), by=0.1),
                             year = sort(unique(dat_use$year)),
                             k = median(dat_use$k),
                             Y1000_scaled = 0,
                             yday=median(dat_use$yday),
                             Area = "Norwegian", 
                             sex=2)
    

    X_proj_length <- mgcv::predict.gam(mgcv_mod, type = "lpmatrix", newdata = pred_grid_length)

    A_proj_length <- INLA::inla.spde.make.A(bspde$mesh,
                                    loc = as.matrix(pred_grid_length[which(pred_grid_length$year == 2019), c("X1000", "Y1000"), drop = FALSE]))
    tmb_data_length <- list(
      X                = as.matrix(X),
      yobs             = as.numeric(yobs),
      Nobs             = Nobs,
      spde_barrier     = make_barrier_spde(bspde),
      spde_barrier1994 = make_barrier_spde(bspde),
      spde_barrier1995 = make_barrier_spde(bspde),
      spde_barrier1996 = make_barrier_spde(bspde),
      spde_barrier1997 = make_barrier_spde(bspde),
      spde_barrier1998 = make_barrier_spde(bspde),
      spde_barrier1999 = make_barrier_spde(bspde),
      spde_barrier2000 = make_barrier_spde(bspde),
      spde_barrier2001 = make_barrier_spde(bspde),
      spde_barrier2002 = make_barrier_spde(bspde),
      spde_barrier2003 = make_barrier_spde(bspde),
      spde_barrier2004 = make_barrier_spde(bspde),
      spde_barrier2005 = make_barrier_spde(bspde),
      spde_barrier2006 = make_barrier_spde(bspde),
      spde_barrier2007 = make_barrier_spde(bspde),
      spde_barrier2008 = make_barrier_spde(bspde),
      spde_barrier2009 = make_barrier_spde(bspde),
      spde_barrier2010 = make_barrier_spde(bspde),
      spde_barrier2011 = make_barrier_spde(bspde),
      spde_barrier2012 = make_barrier_spde(bspde),
      spde_barrier2013 = make_barrier_spde(bspde),
      spde_barrier2014 = make_barrier_spde(bspde),
      spde_barrier2015 = make_barrier_spde(bspde),
      spde_barrier2016 = make_barrier_spde(bspde),
      spde_barrier2017 = make_barrier_spde(bspde),
      spde_barrier2018 = make_barrier_spde(bspde),
      spde_barrier2019 = make_barrier_spde(bspde),
      barrier_scaling  = bspde$barrier_scaling,
      Aobs             = A,
      Ast              = A_st,
      A_spatial_index  = dat_use$sdm_spatial_id - 1L,
      year_i           = make_year_i(dat_use$year),
      Nyear            = Nyear,
      Nmesh            = Nmesh,
      Npred            = nrow(A_proj_length),
      do_predict       = 1,
      calc_se					 = 1,
      X_proj					 = X_proj_length,
      A_proj					 = A_proj_length,
      family           = family_choice,            # 0: tweedie, 1: gamma, 2: student, 3:lognormal
      link             = family_link,						 # 0: identity, 1: log, 2: logit, 3: inverse
      AR1_spatial      = 0,
      df               = dfs,          # df for the student t distribution
      exclude_RE_pred  = 1             # exclude RE from prediction. Useful if wanting marginal effects
    )

    Map_phase2 <- list()
    if (tmb_data_length$AR1_spatial == 0) Map_phase2$transf_rho <- factor(NA)
    Map_phase2$omega <- factor(rep(NA, length(tmb_params$omega)))
    Map_phase2$logTauO <- factor(NA)
    if(tmb_data_length$family != 0) Map_phase2$thetaf <- factor(NA)

    tmb_obj_phase4 <- TMB::MakeADFun(data = tmb_data_length, parameters = tmb_params, map = Map_phase2,
                                       random = c("epsilon_st"), DLL = "herring_sdm", silent = FALSE)

    marginal_length <- tmb_obj_phase4$report(par_vals)$mu_proj

    qwe <- sdreport(tmb_obj_phase4, par.fixed = opt_phase2$par)
    qqq1 <- summary(qwe, "report")[1:length(seq(min(dat_use$length_use), max(dat_use$length_use), by=0.1)),]
    colnames(qqq1) <- c("mean", "se")
    qqq1 <- as.data.frame(qqq1)
    qqq1 <- qqq1 %>% mutate(lwr=mean - 1.96*se,
                          upr=mean + 1.96*se)
    qqq1$length_use = seq(min(dat_use$length_use), max(dat_use$length_use), by=0.1)
 
    } else {
      qqq1 <- NULL 
    }
    
  ### Saving model output
    model_output$marginal_year <- qqq
    model_output$marginal_length <- qqq1


    #### Save different models
    model_save[[paste0(toanalyze, "_", Period, "_age", Age_use)]] <- model_output

}
}


# save output 
saveRDS(model_save, file=paste0("output/Results_", toanalyze, ".rds"))

