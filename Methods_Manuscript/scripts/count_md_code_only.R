## Packages

library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggnewscale)
library(data.table)
library(igraph)
library(scales)
library(vegan)
library(RColorBrewer)
library(lme4)
library(DHARMa)
library(performance)
library(mvabund)
library(emmeans)
library(glmmTMB)
library(dlookr)
library(multcomp)
library(multcompView)
library(stringr)
library(car)
library(superb)


## Setup

### Count data download


all_pitfall_names <- fread("../raw_data/out.csv", header = FALSE)

setnames(all_pitfall_names, "V1", "FullName")

all_pitfall_names <- all_pitfall_names[grep("debris", FullName, invert = TRUE)]
all_pitfall_names <- all_pitfall_names[grep("multiple_insects", FullName, invert = TRUE)]
all_pitfall_names$FullName <- gsub(" ", "_", all_pitfall_names$FullName)
all_pitfall_names[, c("Root","Filename") := tstrsplit(FullName, "\\/+(?!\\S*\\/)", perl = TRUE)]
all_pitfall_names[, c("Root","Phylum","Class","Order","Family","Subfamily","Genus","Species") := tstrsplit(Root, "/")]
all_pitfall_names[, c("StartDate", "Site", "PrescribedTreatment", "Distance") := tstrsplit(Filename, "_", fixed = TRUE, type.convert = TRUE)[1:4]]
all_pitfall_names[, Distance := strtrim(Distance, 2)]
all_pitfall_names[, StartDate := as.Date(StartDate)]
all_pitfall_names[, id := .I]
all_pitfall_names <- all_pitfall_names[StartDate > as.Date("2023-06-21")]

metatable <- fread("../cleaned_data/metadata.csv")
metatable <- metatable[!is.na(Site)]

pitfall_names <- all_pitfall_names[metatable, Treatment := i.Treatment, on=.(Site, StartDate)]
pitfall_names <- pitfall_names[Treatment %in% c("A","C","W")]
pitfall_names$Treatment <- factor(pitfall_names$Treatment, levels=c("C","A","W"), ordered = FALSE)


## Functions for formatting model outputs

summaries_file <- "../cleaned_data/summaries.csv"
pairs_file <- "../cleaned_data/pairs.csv"

if (file.exists(summaries_file)) {
  file.remove(summaries_file)
}

if (file.exists(pairs_file)) {
  file.remove(pairs_file)
}

write_summary <- function(summary_obj) {
  summary_table <- data.table(names(summary_obj[,1]),
                              sapply(as.data.table(summary_obj), function(N) prettyNum(signif(N, digits = 3))))
  setnames(summary_table, names(summary_table), append(deparse(substitute(summary_obj)),colnames(summary_obj)))
  fwrite(summary_table, file = summaries_file, append = TRUE,  col.names = TRUE)
}

write_anova <- function(anova_obj) {
  anova_table <- data.table(rownames(anova_obj),
                            sapply(as.data.table(anova_obj), function(N) prettyNum(signif(N, digits = 3))))
  setnames(anova_table, names(anova_table), append(deparse(substitute(anova_obj)),names(anova_obj)))
  fwrite(anova_table, file = summaries_file, append = TRUE,  col.names = TRUE)
}

write_pairs <- function(pairs_obj) {
  pairs_table <- as.data.table(pairs_obj)
  numcolnames = pairs_table[, names(.SD), .SDcols = is.numeric]
  pairs_table[, (numcolnames) := lapply(.SD, function(N) prettyNum(signif(N, digits = 3))), .SDcols = is.numeric]
  fwrite(pairs_table, file = pairs_file, append = TRUE,  col.names = TRUE)
}



### Biomass data download

biomass_c <- fread("../cleaned_data/biomass.csv")
biomass_c <- biomass_c[, StartDate := as.Date(StartDate)]
setnames(biomass_c, "Near", "0m")
setnames(biomass_c, "Far", "5m")
biomass_c <- melt(setDT(biomass_c[, !c("Total")]), id.vars = c("StartDate","Site", "Treatment"), variable.name = "Distance")
setnames(biomass_c, "value", "Biomass")

first_julian_biomass <- julian(min(biomass_c$StartDate))
biomass_c[, StartDateJulian := as.numeric(julian(StartDate)) - as.numeric(first_julian_biomass)]

# unique(biomass_c[, (as.numeric(julian(StartDate)) - min(as.numeric(julian(StartDate)))) / (max(as.numeric(julian(StartDate))) - min(as.numeric(julian(StartDate))))])
# 
# unique(wide_counts[, (as.numeric(julian(StartDate)) - min(as.numeric(julian(StartDate)))) / (max(as.numeric(julian(StartDate))) - min(as.numeric(julian(StartDate))))])

#biomass_c[,StartDateJulian := scale(StartDateJulian)]
biomass_c

biomass_zeroes <- biomass_c[Biomass == 0]
biomass_zeroes[, Phylum := "Arthropod"]
setnames(biomass_zeroes, "Biomass", "Count")

merge(pitfall_names, biomass_zeroes, all.x = TRUE, all.y = TRUE)

biomass_c$Treatment <- factor(biomass_c$Treatment, levels=c("C","W","A"), ordered = FALSE)

taxa_cols <- c("Class","Order","Family","Subfamily","Genus","Species")

my_comparisons <- list(c("C","W"), c("C","A"), c("W","A"))


### Add rows for each nested level

wide_counts <- pitfall_names[, !c("FullName", "Root", "Filename", "PrescribedTreatment", "id")][, .N, by=.(Phylum, Class, Order, Family, Subfamily, Genus, Species, Site, Treatment, Distance, StartDate)]
cols <- c("Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species")
f <- function(dt) {
  i <- sum(!is.na(dt[,..cols]))
  if (i > 1) {
    dt <- dt[rep(1L, i)]
    for (j in 2:i) set(dt, i = 1:(j - 1), j = j, NA)
  }
  dt
}
wide_counts <- wide_counts[,f(.SD), 1:nrow(wide_counts)][,nrow := NULL][
  ,.(N = sum(N)), Phylum:StartDate
]

setnames(wide_counts, old = "N", new = "Count")




## Biomass boxplot


pdf("../plots/biomass_box_all.pdf", w = 15, h = 9)
ggplot(bsc[Measure == "Biomass"], aes(Treatment, value)) +
  geom_boxplot() +
  #geom_point(alpha = 0.3, position = position_jitter(seed = 1, width = 0.2)) +
  theme(panel.grid = element_blank()) +
  labs(y = "Biomass (grams)", title = "Biomass captured across treatments") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=22, face = "bold"),
        strip.text = element_text(size = 13)) +
  ylim(0, 15) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("No light", "White", "Amber"))
dev.off()



## Simple biomass model

lme_biomass_simple <- glmmTMB(sqrt(Biomass) ~ Treatment * Distance  + (1|Site) + (1|Site:Distance),
                              data = capture_size)
plot(simulateResiduals(lme_biomass_simple))
summary(lme_biomass_simple)
emmgrid <- regrid(emmeans(lme_biomass_simple, specs=pairwise~Treatment, type="response"))
#emmgrid$contrasts

write_summary(coef(summary(lme_biomass_simple))$cond)

lme_biomass_simple_c <- glmmTMB(sqrt(Biomass) ~ Treatment * Distance  + (1|Site) + (1|Site:Distance),
                                data = capture_size,
                                contrasts = list(Treatment = contr.sum))
Anova(lme_biomass_simple_c, type ="III")
write_anova(Anova(lme_biomass_simple_c, type = "III"))


simple_biomass_emmeans <- as.data.table(regrid(emmeans(lme_biomass_simple, list(pairwise ~ Treatment * Distance), adjust = "tukey", type = "response")))

pdf("../plots/simple_regrid_biomass_emmeans.pdf", w = 15, h = 9)
ggplot(simple_biomass_emmeans, aes(Distance, response, group = Treatment, fill = Treatment)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = response-SE, ymax = response+SE), position = "dodge") +
  scale_fill_manual(values = c("C" = "black",
                               "A" = "orange",
                               "W" = "grey")) +
  theme(panel.grid = element_blank()) +
  labs(y = "Biomass (grams)", title = "Predicted biomass captured across treatments and distances") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=18),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=22, face = "bold"),
        legend.text = element_text(size = 20)) +
  ylim(0, 5)
dev.off()

