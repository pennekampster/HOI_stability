library(gauseR)
library(tidyr)
library(stringr)

library(ggplot2)
library(here)

# for partial regression plots
#devtools::install_github("cardiomoon/ggiraphExtra")
library(ggiraphExtra)

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
time<-ddt$NumDays
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio
species$Colp <- species$Colp

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx)
mod_Colp_CD[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio, data=Colp_mod_dat, na.action = na.omit), conf.int = T, conf.level = 0.95 )
#mod_Colp_CD[[1]]
#plot(mod_Colp_CD[[1]])

# print(ggplot(data=Colp_mod_dat, aes(y=Colp_dNNdt,x=Colp)) + geom_point() + stat_smooth(method="lm"))
# print(ggplot(data=Colp_mod_dat, aes(y=Colp_dNNdt,x=Dexio)) + geom_point() + stat_smooth(method="lm"))

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
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
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio/1000
species$Colp <- species$Colp/1000
species$Para <- species$Para/1000

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
Para_lagged<-get_lag(x = ddt$Para, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                         Para=Para_lagged$laggedx)
mod_Colp_CDP[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Para, data=Colp_mod_dat, na.action = na.omit), conf.int = T)
#summary(mod_Colp_CDP2)

#plot(mod_Colp)

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                          Para=Para_lagged$laggedx)
mod_Dexio_CDP[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Para, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Dexio_CDP)
#plot(mod_Dexio)


Para_dNNdt<-percap_growth(x = Para_lagged$x,
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
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio/1000
species$Colp <- species$Colp/1000
species$Spiro <- species$Spiro/1000

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
Spiro_lagged<-get_lag(x = ddt$Spiro, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                         Spiro=Spiro_lagged$laggedx)
mod_Colp_CDS[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Spiro, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Colp_CDS)
#plot(mod_Colp)

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                          Spiro=Spiro_lagged$laggedx)
mod_Dexio_CDS[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Spiro, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Dexio_CDS)
#plot(mod_Dexio)


Spiro_dNNdt<-percap_growth(x = Spiro_lagged$x,
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
time<-ddt$NumDays
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio/1000
species$Colp <- species$Colp/1000
species$Spath <- species$Spath/1000

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx,
                         predator=predator_lagged$laggedx)
mod_Colp_CDPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Colp_CDPd)
#plot(mod_Colp)

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,
                          predator=predator_lagged$laggedx)
mod_Dexio_CDPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Dexio_CDPd)
#plot(mod_Dexio)

predator_dNNdt<-percap_growth(x = predator_lagged$x,
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
time<-ddt$NumDays
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio/1000
species$Colp <- species$Colp/1000
species$Para <- species$Para/1000
species$Spath <- species$Spath/1000

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
Para_lagged<-get_lag(x = ddt$Para, time = ddt$NumDays)
predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Para=Para_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Colp_CDPPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Para+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Colp_CDPPd)
#plot(mod_Colp)

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,Para=Para_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Dexio_CDPPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Para+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Dexio_CDPPd)
#plot(mod_Dexio)


Para_dNNdt<-percap_growth(x = Para_lagged$x, laggedx = Para_lagged$laggedx, dt = Para_lagged$dt)
Para_mod_dat<-data.frame(Para_dNNdt=Para_dNNdt, Para=Para_lagged$laggedx, Colp=Colp_lagged$laggedx,Dexio=Dexio_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Para_CDPPd[[i]] <- broom::tidy(lm(Para_dNNdt~Para+Colp+Dexio+predator, data=Para_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Para_CDPPd)
#plot(mod_Para)

predator_dNNdt<-percap_growth(x = predator_lagged$x,
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
time<-ddt$NumDays
species <- as.data.frame(ddt[,3:ncol(ddt)])
species$Dexio <- species$Dexio/1000
species$Colp <- species$Colp/1000
species$Spiro <- species$Spiro/1000
species$Spath <- species$Spath/1000

# get time-lagged observations for each species
Colp_lagged<-get_lag(x = ddt$Colp, time = ddt$NumDays)
Dexio_lagged<-get_lag(x = ddt$Dexio, time = ddt$NumDays)
Spiro_lagged<-get_lag(x = ddt$Spiro, time = ddt$NumDays)
predator_lagged<-get_lag(x = ddt$Spath, time = ddt$NumDays)

# calculate per-capita growth rates
Colp_dNNdt<-percap_growth(x = Colp_lagged$x, laggedx = Colp_lagged$laggedx, dt = Colp_lagged$dt)
Colp_mod_dat<-data.frame(Colp_dNNdt=Colp_dNNdt, Colp=Colp_lagged$laggedx, Dexio=Dexio_lagged$laggedx, Spiro=Spiro_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Colp_CDSPd[[i]] <- broom::tidy(lm(Colp_dNNdt~Colp+Dexio+Spiro+predator, data=Colp_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Colp_CDSPd)
#plot(mod_Colp)

Dexio_dNNdt<-percap_growth(x = Dexio_lagged$x, laggedx = Dexio_lagged$laggedx, dt = Dexio_lagged$dt)
Dexio_mod_dat<-data.frame(Dexio_dNNdt=Dexio_dNNdt, Dexio=Dexio_lagged$laggedx, Colp=Colp_lagged$laggedx,Spiro=Spiro_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Dexio_CDSPd[[i]] <- broom::tidy(lm(Dexio_dNNdt~Dexio+Colp+Spiro+predator, data=Dexio_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Dexio_CDSPd)
#plot(mod_Dexio)

Spiro_dNNdt<-percap_growth(x = Spiro_lagged$x, laggedx = Spiro_lagged$laggedx, dt = Spiro_lagged$dt)
Spiro_mod_dat<-data.frame(Spiro_dNNdt=Spiro_dNNdt, Spiro=Spiro_lagged$laggedx, Colp=Colp_lagged$laggedx,Dexio=Dexio_lagged$laggedx, predator=predator_lagged$laggedx)
mod_Spiro_CDSPd[[i]] <- broom::tidy(lm(Spiro_dNNdt~Spiro+Colp+Dexio+predator, data=Spiro_mod_dat, na.action = na.omit), conf.int=T)
#summary(mod_Spiro_CDSPd)
#plot(mod_Spiro)

predator_dNNdt<-percap_growth(x = predator_lagged$x,
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

# find all Dexio models
mod_Dexio_names <- ls(pattern = "^mod_Dexio_CD*")
mod_Dexio_list <- lapply(mod_Dexio_names, get)
names(mod_Dexio_list) <- mod_Dexio_names
mod_Dexio_df <- bind_rows(mod_Dexio_list)
mod_Dexio_df$combination = c(rep("CD", 12), rep("CDP", 16), rep("CDPd", 16), rep("CDPPd", 20), rep("CDS", 16), rep("CDSPd", 20))
mod_Dexio_df$predation = ifelse(grepl("Pd", mod_Dexio_df$combination), "yes", "no")
mod_Dexio_df$replicate = c(rep(1:4, each=3), rep(1:4, each=4), rep(1:4, each=4), rep(1:4, each=5), rep(1:4, each=4), rep(1:4, each=5))

# find all Colp models
mod_Colp_names <- ls(pattern = "^mod_Colp_CD*")
mod_Colp_list <- lapply(mod_Colp_names, get)
names(mod_Colp_list) <- mod_Colp_names
mod_Colp_df <- bind_rows(mod_Colp_list)
mod_Colp_df$combination = c(rep("CD", 12), rep("CDP", 16), rep("CDPd", 16), rep("CDPPd", 20), rep("CDS", 16), rep("CDSPd", 20))
mod_Colp_df$predation = ifelse(grepl("Pd", mod_Colp_df$combination), "yes", "no")
mod_Colp_df$replicate = c(rep(1:4, each=3), rep(1:4, each=4), rep(1:4, each=4), rep(1:4, each=5), rep(1:4, each=4), rep(1:4, each=5))
