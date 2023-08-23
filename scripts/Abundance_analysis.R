library(tidyverse)
library(lubridate)
library(here)

library(gtsummary)
library(officer);
library(flextable);

replicate_mean <- read_csv(here("data/Abundance_replicate_mean.csv"))
replicate_mean <- replicate_mean %>% mutate(species = case_when(predict_spec == "Colp" ~ "Colpidium",
                           predict_spec == "Dexio" ~ "Dexiostoma",
                           predict_spec == "Para" ~ "Paramecium",
                           predict_spec == "Spiro" ~ "Spirostomum",
                           predict_spec == "Spath" ~ "Spathidium"))

replicate_mean$treatment <- factor(replicate_mean$treatment, levels=c("C+D", "C+D+P", "C+D+S", "C+D+Pd", "C+D+P+Pd", "C+D+S+Pd"))

# expand the dataset to include all dates where no density was detected and hence abundance = 0.
exp_design <- replicate_mean %>% filter(NumDays==1) %>% ungroup %>% select(replicate, treatment, predict_spec, species)
date_df <- data.frame(NumDays=seq(1,53, by=2))
exp_design <- full_join(exp_design, date_df, by=character())

replicate_mean <- merge(exp_design, replicate_mean, all.x = T)
replicate_mean$mean.dens.ml <- ifelse(is.na(replicate_mean$mean.dens.ml), 0, replicate_mean$mean.dens.ml)
replicate_mean$species <- factor(replicate_mean$species, levels=c("Colpidium", "Dexiostoma", "Paramecium", "Spirostomum", "Spathidium"))

ggplot(data=subset(replicate_mean, predict_spec != "Noise"), aes(x=NumDays, y=mean.dens.ml+1, group=species, colour=species)) + geom_point() +
  geom_path(aes(group=species)) + facet_grid(treatment~replicate) + scale_y_continuous(trans = "sqrt", breaks=c(0,100,200,400, 600))  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  ylab("Mean density per mL") + xlab("Time (in days)") + guides(colour = guide_legend("Species")) + scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50))
ggsave(here("MS_figures/figure2.png"), width=7, height=10)
ggsave(here("MS_figures/figure2.pdf"), dpi=600, width=3800, height=5400, units = "px")


replicate_mean <- replicate_mean %>% filter(predict_spec != "Noise") %>% 
  group_by(treatment, replicate, predict_spec) %>% 
  summarize(mean_dens = mean(mean.dens.ml, na.rm=T))

#replicate_mean <- bind_rows(mono, replicate_mean)

replicate_mean <- replicate_mean %>% mutate(predation =
  case_when(grepl("Pd", treatment) ~ "yes",
            TRUE ~ "no")) %>%
            mutate(competition =
           case_when(!grepl("Pd", treatment) ~ "yes",
                     TRUE ~ "no")) %>%
           mutate(competitor =
           case_when(treatment %in% c("C+D", "C+D+Pd") ~ "none",
                     treatment %in% c("C+D+S", "C+D+S+Pd") ~ "Spiro",
                     treatment %in% c("C+D+P", "C+D+P+Pd") ~ "Para"))
replicate_mean$competitor <- factor(replicate_mean$competitor, level=c("none", "Colp/Dexio", "Para", "Spiro"))
replicate_mean$interaction <- ifelse(replicate_mean$predation == "yes", "With predator and competition", "only competition")

replicate_mean <- replicate_mean %>% mutate(species = case_when(predict_spec == "Colp" ~ "Colpidium",
                                                                            predict_spec == "Dexio" ~ "Dexiostoma",
                                                                            predict_spec == "Para" ~ "Paramecium",
                                                                            predict_spec == "Spiro" ~ "Spirostomum",
                                                                            predict_spec == "Spath" ~ "Spathidium"),
                                            treatment2 =
                                              case_when(treatment == "C+D" ~ "CD",
                                                        treatment == "C+D+P" ~ "CDP",
                                                        treatment == "C+D+S" ~ "CDS",
                                                        treatment == "C+D+Pd" ~ "CD",
                                                        treatment == "C+D+P+Pd" ~ "CDP",
                                                        treatment == "C+D+S+Pd" ~ "CDS"))
replicate_mean$treatment2 <- factor(replicate_mean$treatment2, levels=c("CD", "CDP", "CDS"))

# split competition and predation treatments
replicate_mean_competition <- subset(replicate_mean, predation == "no")
replicate_mean_competition <- replicate_mean_competition %>% mutate(species =
                    case_when(predict_spec == "Colp" ~ "Colpidium",
                              predict_spec == "Dexio" ~ "Dexiostoma",
                              predict_spec == "Para" ~ "Paramecium",
                              predict_spec == "Spiro" ~ "Spirostomum"),
                    treatment2 =
                      case_when(treatment == "C+D" ~ "CD",
                                treatment == "C+D+P" ~ "CDP",
                                treatment == "C+D+S" ~ "CDS"))
         

replicate_mean_predation <- subset(replicate_mean, predation == "yes")
replicate_mean_predation$treatment <- factor(replicate_mean_predation$treatment, levels=c("C+D+Pd", "C+D+P+Pd", "C+D+S+Pd"))
replicate_mean_predation <- replicate_mean_predation %>% mutate(species =
                                                                      case_when(predict_spec == "Colp" ~ "Colpidium",
                                                                                predict_spec == "Dexio" ~ "Dexiostoma",
                                                                                predict_spec == "Para" ~ "Paramecium",
                                                                                predict_spec == "Spiro" ~ "Spirostomum",
                                                                                predict_spec == "Spath" ~ "Spathidium"),
                                                                    treatment2 =
                                                                      case_when(treatment == "C+D+Pd" ~ "focal pair and predator",
                                                                                treatment == "C+D+P+Pd" ~ "focal pair, Paramecium and predator",
                                                                                treatment == "C+D+S+Pd" ~ "focal pair, Spirostomum and predator"))


         
 
plot_df <- subset(replicate_mean)
plot_df$species <- factor(plot_df$species, levels=c("Colpidium", "Dexiostoma", "Paramecium", "Spirostomum", "Spathidium"))

ggplot(data=plot_df, aes(x=predation, y=log(mean_dens), group=treatment, colour=treatment2, shape = interaction)) +  
  geom_point(position = position_dodge(width = 0.4), alpha=.5) +
  stat_summary(
    geom = "pointrange",
    fun.data = mean_se,
    position = position_dodge(width = 0.4)
  ) + scale_y_log10()  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c("bottom")) +
  ylab("log density per mL") + xlab("Predation") +  
  guides(
    colour = guide_legend("Food web configuration"),
    shape = guide_none("Interaction")
  ) + facet_wrap(.~species, ncol=5)
ggsave(here("MS_figures/figure3.png"), width=8, height=4)
ggsave(here("MS_figures/figure3.pdf"), dpi=600, width=4000, height=2000, units = "px")


mod_dexio <- lm(log(mean_dens)~predation*competitor, data=subset(replicate_mean, predict_spec == "Dexio"))
mod_colp <- lm(log(mean_dens)~predation*competitor, data=subset(replicate_mean, predict_spec == "Colp"))

mod_spath <- (lm(log(mean_dens)~competitor, data=subset(replicate_mean, predict_spec == "Spath")))
mod_para <-(lm(log(mean_dens)~predation, data=subset(replicate_mean, predict_spec == "Para")))
mod_spiro <- (lm(log(mean_dens)~predation, data=subset(replicate_mean, predict_spec == "Spiro")))

# model diagnostics
par(mfrow=c(2,2))
plot(lm(log(mean_dens)~predation*competitor, data=subset(replicate_mean, predict_spec == "Dexio")))
plot(lm(log(mean_dens)~predation*competitor, data=subset(replicate_mean, predict_spec == "Colp")))
plot(lm(log(mean_dens)~competitor, data=subset(replicate_mean, predict_spec == "Spath")))
plot(lm(log(mean_dens)~treatment, data=subset(replicate_mean, predict_spec == "Para")))
plot(lm(log(mean_dens)~treatment, data=subset(replicate_mean, predict_spec == "Spiro")))
par(mfrow=c(1,1))


ft1 <- tbl_regression(mod_dexio, intercept=T) %>% as_flex_table()
flextable::save_as_docx(ft1, path = here("MS_tables/dexio_table.docx"))

ft2 <- tbl_regression(mod_colp, intercept=T)%>% as_flex_table()
flextable::save_as_docx(ft2, path = here("MS_tables/colp_table.docx"))

ft3 <- tbl_regression(mod_spath, intercept=T)%>% as_flex_table()
flextable::save_as_docx(ft3, path = here("MS_tables/spath_table.docx"))

ft4 <- tbl_regression(mod_para, intercept=T)%>% as_flex_table()
flextable::save_as_docx(ft4, path = here("MS_tables/para_table.docx"))

ft5 <- tbl_regression(mod_spiro, intercept=T)%>% as_flex_table()
flextable::save_as_docx(ft5, path = here("MS_tables/spiro_table.docx"))

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", width = 15, height = 8.3),
  type = "continuous",
  page_margins = page_mar()
)


tbl_merge_ex1 <-
  tbl_merge(
    tbls = list(tbl_regression(mod_colp, intercept=T), tbl_regression(mod_dexio, intercept=T), 
                tbl_regression(mod_para, intercept=T), tbl_regression(mod_spiro, intercept=T), 
                tbl_regression(mod_spath, intercept=T)), 
    tab_spanner = c("**Colpidium**", "**Dexiostoma**", "**Paramecium**", "**Spirostomum**", "**Spathidium**")
  )
tbl_merge_ex1
flextable::save_as_docx(tbl_merge_ex1 %>% as_flex_table(), path = here("MS_tables/Abundance_ANOVA_table.docx"), pr_section = sect_properties)
flextable::save_as_html(tbl_merge_ex1 %>% as_flex_table(), path = here("MS_tables/Abundance_ANOVA_table.html"))


