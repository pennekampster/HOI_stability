library(tidyverse)
library(effects)

library(interactions)
library(visreg)
library(cowplot)

library(sjPlot)
theme_set(theme_sjplot())
library(marginaleffects)

library(gauseR)
library(tidyr)
library(stringr)

library(ggplot2)
library(here)

# for partial regression plots
#devtools::install_github("cardiomoon/ggiraphExtra")
library(ggiraphExtra)

# function adapted from gauseR package
percap_growth_notlogged <- function (x, laggedx, dt) 
{
  xrat <- unname(x/laggedx)
  xrat[xrat <= 0 | !is.finite(xrat)] <- NA
  return((xrat) * 1/(dt))
}

load(here("data/Abundance_replicate_mean.Rdata"))

# expand the dataset to include all dates where no density was detected and hence abundance = 0.
exp_design <- replicate_mean %>% filter(NumDays==1) %>% ungroup %>% dplyr::select(replicate, treatment, predict_spec)
date_df <- data.frame(NumDays=seq(1,53, by=2))
exp_design <- full_join(exp_design, date_df, by=character())

replicate_mean <- merge(exp_design, replicate_mean, all.x = T)
replicate_mean$mean.dens.ml <- ifelse(is.na(replicate_mean$mean.dens.ml), 0, replicate_mean$mean.dens.ml)

# define the day range over which interaction strength is measured
NumDaysMax <- 53

mod_Dexio_CD <- mod_Colp_CD <- vector("list", 4)
dat_Dexio_CD <- dat_Colp_CD <- vector("list", 4)

for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D" & NumDays < NumDaysMax & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]

  # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx)
  mod_Colp_CD[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio, data=Colp_mod_dat, na.action = na.omit), conf.int = T, conf.level = 0.95 )
  #mod_Colp_CD[[1]]
  #plot(mod_Colp_CD[[1]])
  
  # print(ggplot(data=Colp_mod_dat, aes(y=Colp_dNNdt,x=Colp)) + geom_point() + stat_smooth(method="lm"))
  # print(ggplot(data=Colp_mod_dat, aes(y=Colp_dNNdt,x=Dexio)) + geom_point() + stat_smooth(method="lm"))
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx)
  mod_Dexio_CD[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp, data=Dexio_mod_dat, na.action = na.omit), conf.int = T, conf.level = 0.95 )
  
  # keep data for mixed models
  Dexio_mod_dat$replicate <- i
  Dexio_mod_dat$com <- "CD"
  dat_Dexio_CD[[i]] <- Dexio_mod_dat
  
  Colp_mod_dat$replicate <- i
  Colp_mod_dat$com <- "CD"
  
  dat_Colp_CD[[i]] <- Colp_mod_dat
  
  
  #summary(mod_Dexio_CD2)
  #plot(mod_Dexio)
}


Dexio_CD_df <- bind_rows(mod_Dexio_CD)
Dexio_CD_df$target_species <- "Dexio"
Dexio_CD_df$replicate <- rep(1:4, each=3)
Colp_CD_df <- bind_rows(mod_Colp_CD)
Colp_CD_df$target_species <- "Colp"
Colp_CD_df$replicate <- rep(1:4, each=3)


Dexio_CD_dat_df <- bind_rows(dat_Dexio_CD)
Colp_CD_dat_df <- bind_rows(dat_Colp_CD)


mod_Dexio_CDP <- mod_Colp_CDP <- mod_Para_CDP <-  vector("list", 4)
dat_Dexio_CDP <- dat_Colp_CDP <- dat_Para_CDP <-  vector("list", 4)


for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D+P" & NumDays < NumDaysMax & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]
  time<-ddt$NumDays

  # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  Para_lagged<-get_lag(x = ddt$Para, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                           Para=Para_lagged$laggedx)
  mod_Colp_CDP[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Para, data=Colp_mod_dat, na.action = na.omit), conf.int = T)
  #summary(mod_Colp_CDP2)
  
  #plot(mod_Colp)
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                            Para=Para_lagged$laggedx)
  mod_Dexio_CDP[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Para, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Dexio_CDP)
  #plot(mod_Dexio)
  
  
  Para_dNNdt<-percap_growth_notlogged(x = Para_lagged$x,
                            laggedx = Para_lagged$laggedx, dt = Para_lagged$dt)
  Para_mod_dat<-data.frame(Para_dNNdt=Para_dNNdt,
                           Para=Para_lagged$laggedx, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx)
  mod_Para_CDP[[i]] <- broom::tidy(lm(Para_dNNdt~Para+Colp+Dexio, data=Para_mod_dat), conf.int = T)
  #summary(mod_Para_CDP)
  
  
  # keep data for mixed models
  Dexio_mod_dat$replicate <- i
  Dexio_mod_dat$com <- "CDP"
  dat_Dexio_CDP[[i]] <- Dexio_mod_dat
  
  Colp_mod_dat$replicate <- i
  Colp_mod_dat$com <- "CDP"
  dat_Colp_CDP[[i]] <- Colp_mod_dat
  
  
}

Dexio_CDP_df <- bind_rows(mod_Dexio_CDP)
Dexio_CDP_df$target_species <- "Dexio"
Dexio_CDP_df$replicate <- rep(1:4, each=4)
Colp_CDP_df <- bind_rows(mod_Colp_CDP)
Colp_CDP_df$target_species <- "Colp"
Colp_CDP_df$replicate <- rep(1:4, each=4)
Para_CDP_df <- bind_rows(mod_Para_CDP)
Para_CDP_df$target_species <- "Para"
Para_CDP_df$replicate <- rep(1:4, each=4)

Dexio_CDP_dat_df <- bind_rows(dat_Dexio_CDP)
Colp_CDP_dat_df <- bind_rows(dat_Colp_CDP)



mod_Dexio_CDS <- mod_Colp_CDS <- mod_Spiro_CDS <- vector("list", 4)
dat_Dexio_CDS <- dat_Colp_CDS <- dat_Spiro_CDS <- vector("list", 4)

for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D+S" & NumDays < NumDaysMax & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]
  time<-ddt$NumDays

    # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  Spiro_lagged<-get_lag(x = ddt$Spiro, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                           Spiro=Spiro_lagged$laggedx)
  mod_Colp_CDS[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Spiro, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Colp_CDS)
  #plot(mod_Colp)
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                            Spiro=Spiro_lagged$laggedx)
  mod_Dexio_CDS[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Spiro, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Dexio_CDS)
  #plot(mod_Dexio)
  
  
  Spiro_dNNdt<-percap_growth_notlogged(x = Spiro_lagged$x,
                             laggedx = Spiro_lagged$laggedx, dt = Spiro_lagged$dt)
  Spiro_mod_dat<-data.frame(Spiro_dNNdt=Spiro_dNNdt,
                            Spiro=Spiro_lagged$laggedx, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx)
  mod_Spiro_CDS[[i]] <- broom::tidy(lm(Spiro_dNNdt~Spiro+Colp+Dexio, data=Spiro_mod_dat), conf.int=T)
  #summary(mod_Spiro_CDS)
  
  
  # keep data for mixed models
  Dexio_mod_dat$replicate <- i
  Dexio_mod_dat$com <- "CDS"
  dat_Dexio_CDS[[i]] <- Dexio_mod_dat
  Colp_mod_dat$replicate <- i
  Colp_mod_dat$com <- "CDS"
  dat_Colp_CDS[[i]] <- Colp_mod_dat
}

Dexio_CDS_df <- bind_rows(mod_Dexio_CDS)
Dexio_CDS_df$target_species <- "Dexio"
Dexio_CDS_df$replicate <- rep(1:4, each=4)
Colp_CDS_df <- bind_rows(mod_Colp_CDS)
Colp_CDS_df$target_species <- "Colp"
Colp_CDS_df$replicate <- rep(1:4, each=4)
Spiro_CDS_df <- bind_rows(mod_Spiro_CDS)
Spiro_CDS_df$target_species <- "Spiro"
Spiro_CDS_df$replicate <- rep(1:4, each=4)


Dexio_CDS_dat_df <- bind_rows(dat_Dexio_CDS)
Colp_CDS_dat_df <- bind_rows(dat_Colp_CDS)

Dexio_com_dat_df <- bind_rows(Dexio_CD_dat_df, Dexio_CDP_dat_df, Dexio_CDS_dat_df)
Colp_com_dat_df <- bind_rows(Colp_CD_dat_df, Colp_CDP_dat_df, Colp_CDS_dat_df)

Dexio_com_dat_df$Para <- ifelse(is.na(Dexio_com_dat_df$Para), 0, Dexio_com_dat_df$Para)
Dexio_com_dat_df$Spiro <- ifelse(is.na(Dexio_com_dat_df$Spiro), 0, Dexio_com_dat_df$Spiro)

Dexio_com_dat_df$com <- as.factor(Dexio_com_dat_df$com)



mod_Dexio_CDPd <- mod_Colp_CDPd <- mod_predator_CDPd <- vector("list", 4)

for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D+Pd" & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]

  # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                           predator=predator_lagged$laggedx)
  mod_Colp_CDPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Colp_CDPd)
  #plot(mod_Colp)
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                            predator=predator_lagged$laggedx)
  mod_Dexio_CDPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Dexio_CDPd)
  #plot(mod_Dexio)
  
  predator_dNNdt<-percap_growth_notlogged(x = predator_lagged$x,
                                laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)
  predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
                               predator=predator_lagged$laggedx, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx)
  mod_predator_CDPd[[i]] <- broom::tidy(lm(predator_dNNdt~predator+Colp+Dexio, data=predator_mod_dat), conf.int=T)
  #summary(mod_predator_CDPd)
  #plot(mod_predator)
}


Dexio_CDPd_df <- bind_rows(mod_Dexio_CDPd)
Dexio_CDPd_df$target_species <- "Dexio"
Dexio_CDPd_df$replicate <- rep(1:4, each=4)
Colp_CDPd_df <- bind_rows(mod_Colp_CDPd)
Colp_CDPd_df$target_species <- "Colp"
Colp_CDPd_df$replicate <- rep(1:4, each=4)




mod_Dexio_CDPPd <- mod_Colp_CDPPd <- mod_Para_CDPPd <- mod_predator_CDPPd <- vector("list", 4)

for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D+P+Pd" & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]

  # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  Para_lagged<-get_lag(x = ddt$Para, time = ddt$NumDays)
  predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Para=Para_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Colp_CDPPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Para+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Colp_CDPPd)
  #plot(mod_Colp)
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,Para=Para_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Dexio_CDPPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Para+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Dexio_CDPPd)
  #plot(mod_Dexio)
  
  
  Para_dNNdt<-percap_growth_notlogged(x = Para_lagged$x, laggedx = Para_lagged$laggedx, dt = Para_lagged$dt)
  Para_mod_dat<-data.frame(Para_dNNdt=Para_dNNdt, Para=Para_lagged$laggedx, Colp=Colp_lagged$laggedx,Dexio=Dexio_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Para_CDPPd[[i]] <- broom::tidy(lm(Para_dNNdt~Para+Colp+Dexio+predator, data=Para_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Para_CDPPd)
  #plot(mod_Para)
  
  predator_dNNdt<-percap_growth_notlogged(x = predator_lagged$x,
                                laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)
  predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
                               predator=predator_lagged$laggedx, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Para=Para_lagged$laggedx)
  mod_predator_CDPPd[[i]] <- broom::tidy(lm(predator_dNNdt~predator+Colp+Dexio+Para, data=predator_mod_dat), conf.int=T)
  #summary(mod_predator_CDPPd)
  #plot(mod_predator)
}

Dexio_CDPPd_df <- bind_rows(mod_Dexio_CDPPd)
Dexio_CDPPd_df$target_species <- "Dexio"
Dexio_CDPPd_df$replicate <- rep(1:4, each=5)
Colp_CDPPd_df <- bind_rows(mod_Colp_CDPPd)
Colp_CDPPd_df$target_species <- "Colp"
Colp_CDPPd_df$replicate <- rep(1:4, each=5)

mod_Dexio_CDSPd <- mod_Colp_CDSPd <- mod_Spiro_CDSPd <- mod_predator_CDSPd <- vector("list", 4)

for (i in 1:4){
  
  # manual fitting of dynamics
  dd <- subset(replicate_mean, treatment == "C+D+S+Pd" & replicate == i)
  ddt <- pivot_wider(dd, 
                     names_from = "predict_spec",
                     values_from = "mean.dens.ml")
  ddt$Noise <- NULL
  ddt <- ddt[complete.cases(ddt),]

  # get time-lagged observations for each species
  Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
  Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
  Spiro_lagged<-get_lag(x = ddt$Spiro, time = ddt$NumDays)
  predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)
  
  # calculate per-capita growth rates
  Colp_dNNdt<-percap_growth_notlogged(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
  Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Spiro=Spiro_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Colp_CDSPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Spiro+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Colp_CDSPd)
  #plot(mod_Colp)
  
  Dexio_dNNdt<-percap_growth_notlogged(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
  Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,Spiro=Spiro_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Dexio_CDSPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Spiro+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Dexio_CDSPd)
  #plot(mod_Dexio)
  
  Spiro_dNNdt<-percap_growth_notlogged(x = Spiro_lagged$x, laggedx = Spiro_lagged$laggedx, dt = Spiro_lagged$dt)
  Spiro_mod_dat<-data.frame(Spiro_dNNdt=Spiro_dNNdt, Spiro=Spiro_lagged$laggedx, Colp=Colp_lagged$laggedx,Dexio=Dexio_lagged$laggedx, predator=predator_lagged$laggedx)
  mod_Spiro_CDSPd[[i]] <- broom::tidy(lm(Spiro_dNNdt~Spiro+Colp+Dexio+predator, data=Spiro_mod_dat, na.action = na.omit), conf.int=T)
  #summary(mod_Spiro_CDSPd)
  #plot(mod_Spiro)
  
  predator_dNNdt<-percap_growth_notlogged(x = predator_lagged$x,
                                laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)
  predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
                               predator=predator_lagged$laggedx, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Spiro=Spiro_lagged$laggedx)
  mod_predator_CDSPd[[i]] <- broom::tidy(lm(predator_dNNdt~predator+Colp+Dexio+Spiro, data=predator_mod_dat), conf.int=T)
  #summary(mod_predator_CDSPd)
  #plot(mod_predator)
}

Dexio_CDSPd_df <- bind_rows(mod_Dexio_CDSPd)
Dexio_CDSPd_df$target_species <- "Dexio"
Dexio_CDSPd_df$replicate <- rep(1:4, each=5)
Colp_CDSPd_df <- bind_rows(mod_Colp_CDSPd)
Colp_CDSPd_df$target_species <- "Colp"
Colp_CDSPd_df$replicate <- rep(1:4, each=5)

# Lotka-Volterra
m1 <- glm(Dexio_dNNdt~Dexio + Colp, data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian)
m2 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2), data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian)
m3 <- glm(Dexio_dNNdt~Dexio + Colp + Dexio:Colp, data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian)
m4 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2) + Dexio:Colp, data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian)
# Ricker
m5 <- glm(Dexio_dNNdt~Dexio + Colp, data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Dexio_dNNdt~Dexio + Colp + Dexio:Colp, data=Dexio_CD_dat_df[complete.cases(Dexio_CD_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)

modelsummary::modelsummary(models, output = here("other_output/model_output_Dexio_GR_CD.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Dexio_GR_CD.png"), bg = "white", width=7, height=5)

AIC_Dexio_CD <- AICcmodavg::aictab(models, second.ord = F)
plot_model(m4, type = "pred", terms = c("Dexio[0:400]", "Colp"),  show.data = F)
ggsave(here("other_output/Dexio_GR_CD.png"), bg = "white", height=5, width=7)
best_Dexio_CD <- m4

# Colp in CD

m1 <- glm(Colp_dNNdt~Colp + Dexio, data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian)
m2 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2), data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian)
m3 <- glm(Colp_dNNdt~Colp + Dexio + Colp:Dexio, data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian)
m4 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2) + Colp:Dexio, data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian)
m5 <- glm(Colp_dNNdt~Colp + Dexio, data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Colp_dNNdt~Colp + Dexio + Colp:Dexio, data=Colp_CD_dat_df[complete.cases(Colp_CD_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)

modelsummary::modelsummary(models, output = here("other_output/model_output_Colp_GR_CD.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Colp_GR_CD.png"), bg = "white", width=7, height=5)

AIC_Colp_CD <- AICcmodavg::aictab(models, second.ord = F)
plot_model(m4, type = "pred", terms = c("Colp", "Dexio"),  show.data = F)
ggsave(here("other_output/Colp_GR_CD.png"), width=7,height=5, bg = "white")
best_Colp_CD <- m4

# Dexio in CDP

m1 <- glm(Dexio_dNNdt~Dexio + Colp + Para, data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian)
m2 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2) + Para + I(Para^2), data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian)
m3 <- glm(Dexio_dNNdt~Dexio + Colp + Para + Dexio:Colp + Dexio:Para + Colp:Para, data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian)
m4 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2) + Para + I(Para^2) + Dexio:Colp + Dexio:Para + Colp:Para, data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian)
m5 <- glm(Dexio_dNNdt~Dexio + Colp + Para, data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Dexio_dNNdt~Dexio + Colp + Para + Dexio:Colp + Dexio:Para + Colp:Para, data=Dexio_CDP_dat_df[complete.cases(Dexio_CDP_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)

modelsummary::modelsummary(models, output = here("other_output/model_output_Dexio_GR_CDP.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Dexio_GR_CDP.png"), bg = "white", width=10, height=7)

AIC_Dexio_CDP <- AICcmodavg::aictab(models, second.ord = F)
plot_model(m5, type = "pred", terms = c("Dexio", "Colp", "Para"))
ggsave(here("other_output/Dexio_GR_CDP.png"), width=7,height=5, bg = "white")
best_Dexio_CDP <- m5

# Colp in CDP

m1 <- glm(Colp_dNNdt~Colp + Dexio + Para, data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian)
m2 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2)+ Para + I(Para^2), data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian)
m3 <- glm(Colp_dNNdt~Colp + Dexio + Para + Colp:Dexio + Colp:Para + Dexio:Para, data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian)
m4 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2)+ Para + I(Para^2) + Colp:Dexio + Colp:Para + Dexio:Para, data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian)
m5 <- glm(Colp_dNNdt~Colp + Dexio + Para, data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Colp_dNNdt~Colp + Dexio + Para + Colp:Dexio + Colp:Para + Dexio:Para, data=Colp_CDP_dat_df[complete.cases(Colp_CDP_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)

modelsummary::modelsummary(models, output = here("other_output/model_output_Colp_GR_CDP.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Colp_GR_CDP.png"), bg = "white", width=10, height=7)

AIC_Colp_CDP <- AICcmodavg::aictab(models, second.ord = F)
plot_model(m5, type = "pred", terms = c("Colp", "Dexio", "Para"), bg = "white")
ggsave(here("other_output/Colp_GR_CDP.png"), width=7,height=5, bg = "white")
best_Colp_CDP <- m5

# Dexio in CDS

m1 <- glm(Dexio_dNNdt~Dexio + Colp + Spiro, data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian)
m2 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2) + Spiro + I(Spiro^2), data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian)
m3 <- glm(Dexio_dNNdt~Dexio + Colp + Spiro + Dexio:Colp + Dexio:Spiro + Colp:Spiro, data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian)
m4 <- glm(Dexio_dNNdt~Dexio + I(Dexio^2) + Colp + I(Colp^2) + Spiro + I(Spiro^2) + Dexio:Colp + Dexio:Spiro + Colp:Spiro, data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian)
m5 <- glm(Dexio_dNNdt~Dexio + Colp + Spiro, data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Dexio_dNNdt~Dexio + Colp + Spiro + Dexio:Colp + Dexio:Spiro + Colp:Spiro, data=Dexio_CDS_dat_df[complete.cases(Dexio_CDS_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)

modelsummary::modelsummary(models, output = here("other_output/model_output_Dexio_GR_CDS.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Dexio_GR_CDS.png"), bg = "white", width=10, height=7)

AIC_Dexio_CDS <- AICcmodavg::aictab(models, second.ord = F)
plot_model(m5, type = "pred", terms = c("Dexio", "Colp", "Spiro"))
ggsave(here("other_output/Dexio_GR_CDS.png"), width=7,height=5, bg = "white")
best_Dexio_CDS <- m5

# Colp  in CDS

m1 <- glm(Colp_dNNdt~Colp + Dexio + Spiro, data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian)
m2 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2) + Spiro + I(Spiro^2), data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian)
m3 <- glm(Colp_dNNdt~Colp + Dexio + Spiro + Colp:Dexio + Colp:Spiro + Dexio:Spiro, data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian)
m4 <- glm(Colp_dNNdt~Colp + I(Colp^2) + Dexio + I(Dexio^2) + Spiro + I(Spiro^2) + Colp:Dexio + Colp:Spiro + Dexio:Spiro, data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian)
m5 <- glm(Colp_dNNdt~Colp + Dexio + Spiro, data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian(link = "log"))
m6 <- glm(Colp_dNNdt~Colp + Dexio + Spiro + Colp:Dexio + Colp:Spiro + Dexio:Spiro, data=Colp_CDS_dat_df[complete.cases(Colp_CDS_dat_df), ], family = gaussian(link = "log"))

models <- list(
  "additive LV"    = m1,
  "interactive LV (intra HOI)" = m2,
  "interactive LV (inter HOI)" = m3,
  "interactive LV (full HOI)" = m4,
  "additive Ricker"     = m5,
  "interactive Ricker" = m6)


modelsummary::modelsummary(models, output = here("other_output/model_output_Colp_GR_CDS.docx"), fmt = 5, statistic = 'conf.int', gof_map = c("nobs", "aic", "r.squared"))
modelsummary::modelplot(models, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("other_output/model_output_Colp_GR_CDS.png"), bg = "white", width=10, height=7)

AIC_Colp_CDS <- AICcmodavg::aictab(models, second.ord = F)
best_Colp_CDS <- m6
plot_model(best_Colp_CDS, type = "pred", terms = c("Colp", "Dexio", "Spiro"))
ggsave(here("other_output/Colp_GR_CDS.png"), width=7,height=5, bg = "white")


AIC_table <- rbind(AIC_Colp_CD, AIC_Colp_CDP, AIC_Colp_CDS, AIC_Dexio_CD, AIC_Dexio_CDP, AIC_Dexio_CDS)
AIC_table$growth <- rep(c("Colp", "Dexio"), each=18)
AIC_table$community <- rep(c("CD", "CDP", "CDS"), each=6)
write_excel_csv(AIC_table, here("AIC_table_HOIs.csv"))

best_models_Dexio <- list(
  "Dexio CD"    = best_Dexio_CD,
  "Dexio CDP"     = best_Dexio_CDP,
  "Dexio CDS" = best_Dexio_CDS)

modelsummary::modelplot(best_models_Dexio, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("MS_figures/coef_plot_Dexio.png"), width=10,height=6, bg = "white")

best_models_Colp <- list(
  "Colp CD" = best_Colp_CD,
  "Colp CDP" = best_Colp_CDP,
  "Colp CDS" = best_Colp_CDS)

modelsummary::modelplot(best_models_Colp, coef_omit = 'Interc') + facet_wrap(~term, scales = "free") + geom_vline(xintercept=0, linetype="dashed") + theme(axis.text.y = element_blank(), legend.position = "bottom")
ggsave(here("MS_figures/coef_plot_Colp.png"), width=10,height=6, bg = "white")




