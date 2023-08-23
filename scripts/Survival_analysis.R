library("survival")
library("survminer")
library(tidyverse)
library(here)

Extinction <- read_csv(here("data/Extinctions.csv"))
Extinction <- Extinction %>% group_by(treatment, predict_spec, replicate) %>% summarize(extincation_days = mean(extincation_days, na.rm = T))

Extinction$status <- ifelse(Extinction$extincation_days == 53, 0, 1)
Extinction$Predation <- ifelse(grepl("Pd", Extinction$treatment), "yes", "no")
Extinction$Competition <- ifelse(!grepl("Pd", Extinction$treatment), "yes", "no")

Extinction$Competitor <- ifelse(Extinction$treatment %in% c("C+D", "C+D+Pd"), "none", ifelse(Extinction$treatment %in% c("C+D+P", "C+D+P+Pd"), "P", "S"))

fit <- survfit(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Colp"& Predation == "no"))
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table


survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Colp"& Predation == "no"))
pairwise_survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Colp"& Predation == "no"), p.adjust.method = "BH")

p1 <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           legend.title = "treatment",
           legend.labs = c("CD", "CDP", "CDS"),
           pval.method=T,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_survminer(), 
           conf.int.alpha =.1,
           title = "Colpidium (competition)"  # Change ggplot2 theme
           )

fit <- survfit(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Dexio" & Predation == "no"))
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Dexio" & Predation == "no"))
pairwise_survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Dexio" & Predation == "no"), p.adjust.method = "BH")

p2 <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           legend.title = "treatment",
           legend.labs = c("CD", "CDP", "CDS"),
           pval.method=T,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_survminer(), 
           conf.int.alpha =.1,
           title = "Dexiostoma (competition)"  # Change ggplot2 theme
)


fit <- survfit(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = subset(Extinction, predict_spec == "Colp" & Predation == "yes"))
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

surv_data <- subset(Extinction, predict_spec == "Colp" & Predation == "yes")
surv_data <- as.data.frame(surv_data)
survdiff(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = surv_data)
pairwise_survdiff(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = surv_data, p.adjust.method = "BH")

p3 <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
         #  linetype = "strata", # Change line type by groups
           legend.title = "treatment",
           legend.labs = c("CDPd", "CDPPd", "CDSPd"),
           pval.method=T,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_survminer(), 
           conf.int.alpha =.1,
           title = "Colpidium (competition + predation)"  # Change ggplot2 theme
)


fit <- survfit(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = subset(Extinction, predict_spec == "Dexio" & Predation == "yes"))
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

surv_data <- subset(Extinction, predict_spec == "Dexio" & Predation == "yes")
surv_data <- as.data.frame(surv_data)
survdiff(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = surv_data)
pairwise_survdiff(Surv(extincation_days, status, type = "right") ~ Competitor + Predation, data = surv_data, p.adjust.method = "BH")

p4 <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
         #  linetype = "strata", # Change line type by groups
           legend.title = "treatment",
           legend.labs = c("CDPd", "CDPPd", "CDSPd"),
           pval.method=T,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_survminer(), 
           conf.int.alpha =.1,
           title = "Dexiostoma (competition + predation)"  # Change ggplot2 theme
)





fit <- survfit(Surv(extincation_days, status, type = "right") ~ Competitor, data = subset(Extinction, predict_spec == "Spath" & Predation == "yes"))
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

surv_data <- subset(Extinction, predict_spec == "Spath" & Predation == "yes")
surv_data <- as.data.frame(surv_data)
survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = surv_data)
pairwise_survdiff(Surv(extincation_days, status, type = "right") ~ Competitor, data = surv_data, p.adjust.method = "BH")


p5 <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
          # linetype = "strata", # Change line type by groups
           legend.title = "treatment",
           legend.labs = c("CDPd", "CDPPd", "CDSPd"),
           pval.method=T,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_survminer(), 
           conf.int.alpha =.1,
           title = "Spathidium (competition + predation)"  # Change ggplot2 theme
)








splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
splots[[3]] <- p5
splots[[4]] <- p3
splots[[5]] <- p4



# access the plot objects and add a tag with labs()
splots[[1]]$plot<-splots[[1]]$plot + labs(tag="A")
splots[[2]]$plot<-splots[[2]]$plot + labs(tag="C")
splots[[3]]$plot<-splots[[3]]$plot + labs(tag="E")
splots[[4]]$plot<-splots[[4]]$plot + labs(tag="B")
splots[[5]]$plot<-splots[[5]]$plot + labs(tag="D")

# Arrange and save into pdf file
surv_plot <- arrange_ggsurvplots(splots, print = F, ncol = 2, nrow = 3)
ggsave(here("MS_figures/figure4.png"), surv_plot,  width=10, height=15)
ggsave(here("MS_figures/figure4.pdf"), surv_plot, dpi=600, width=6400, height=9200, units = "px")

