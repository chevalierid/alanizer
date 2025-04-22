library(ggplot2)
library(ggpubr)
library(ggsignif)
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

setwd("Documents/alan-setup/prog/refactored")

all_pitfall_names <- fread("./raw_data/pitfall_names_241023.csv", header = FALSE)
#all_pitfall_names <- fread("./raw_data/pitfall_classification/out.csv", header = FALSE)

setnames(all_pitfall_names, "V1", "FullName")
#2023-08-05_1_A_5m002_0040_5763_8740_0004.png
#StartDate_Site_Treatment_Distance

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

metatable <- fread("./cleaned_data/metadata.csv")
metatable <- metatable[!is.na(Site)]

pitfall_names <- all_pitfall_names[metatable, Treatment := i.Treatment, on=.(Site, StartDate)]
pitfall_names <- pitfall_names[Treatment %in% c("A","C","W")]
pitfall_names$Treatment <- factor(pitfall_names$Treatment, levels=c("C","A","W"), ordered = FALSE)

pitfall_names <- pitfall_names[Class != "unclassified_arthropod"]



biomass_c <- fread("./cleaned_data/biomass.csv")
biomass_c <- biomass_c[, StartDate := as.Date(StartDate)]
setnames(biomass_c, "Near", "0m")
setnames(biomass_c, "Far", "5m")
biomass_c <- melt(setDT(biomass_c[, !c("Total")]), id.vars = c("StartDate","Site", "Treatment"), variable.name = "Distance")
setnames(biomass_c, "value", "Biomass")

biomass_zeroes <- biomass_c[Biomass == 0]
biomass_zeroes[, Phylum := "Arthropod"]
setnames(biomass_zeroes, "Biomass", "Count")

merge(pitfall_names, biomass_zeroes, all.x = TRUE, all.y = TRUE)

taxa_cols <- c("Class","Order","Family","Subfamily","Genus","Species")

my_comparisons <- list(c("C","W"), c("C","A"), c("W","A"))

########################
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

first_julian <- julian(min(wide_counts$StartDate))
wide_counts[, StartDateJulian := as.numeric(julian(StartDate)) - as.numeric(first_julian)]
wide_counts[,StartDateJulian := scale(StartDateJulian)]
wide_counts

ggplot(wide_counts[Phylum == "Arthropoda" & is.na(Class)], aes(x = StartDate, y = Count)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar") + 
  labs(y="Mean Total Count") +
  facet_grid(rows = vars(Treatment))

wide_counts$Treatment <- factor(wide_counts$Treatment, levels=c("C","W","A"), ordered = FALSE)

########### BOXPLOTS ###########
capture_size <- wide_counts[Phylum == "Arthropoda" & is.na(Class)]
capture_size <- capture_size[biomass_c, on = .(StartDate, Site, Treatment, Distance)]
capture_size <- capture_size[Phylum == "Arthropoda"]

capture_size[, Size := Biomass/Count * 1000]

capture_size$Treatment <- factor(capture_size$Treatment, levels=c("C","W","A"), ordered = FALSE)

pdf("./plots/size_box_all.pdf", w = 15, h = 9)
ggplot(capture_size, aes(Treatment,Size)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=24, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Larger arthropods (on average) captured under lit conditions") +
  ylab("Size (mg/arthropod)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)
dev.off()

pdf("./plots/size_box_near.pdf", w = 15, h = 9)
ggplot(capture_size[Distance == "0m"], aes(Treatment,Size)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=22, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Size of captured arthropods for closest traps only (biomass divided by count)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter() +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)
dev.off()

########### COUNT BOXPLOTS ###########

pdf("./plots/count_box_all.pdf", w = 15, h = 9)
ggplot(capture_size, aes(Treatment,Count)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),
        title = element_text(size=25, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "More arthropods captured under white light than in darkness") +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)
dev.off()

# pdf("./plots/count_box_near.pdf", w = 15, h = 9)
# ggplot(capture_size[Distance == "0m"], aes(Treatment,Count)) +
#   geom_boxplot(outlier.shape = NA) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=22,face="bold"),
#         title = element_text(size=24, face = "bold")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   labs(title = "Arthropod count per pitfall trap for closest traps only") +
#   geom_jitter() +
#   stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)
# dev.off()

########### FACETED BOXPLOTS ########### 

bsc <- melt(capture_size[, c("Site", "Treatment", "Distance", "StartDateJulian", "Biomass", "Count", "Size")], 
     id.vars = c("Site", "Treatment", "Distance", "StartDateJulian"), variable.name = "Measure")

box_labels <- c('Biomass' = 'Arthropod biomass (g)', 'Count' = 'Total arthropod count', 'Size' = 'Arthropod biomass over count (mg)')

ggplot(bsc, aes(Treatment, value)) +
  geom_boxplot() +
  geom_point(alpha = 0.3, position = position_jitter(seed = 1, width = 0.2)) +
  theme(panel.grid = element_blank()) +
  labs(y="") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=22, face = "bold"),
        strip.text = element_text(size = 13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans='log10') +
  #stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
  facet_grid(rows = vars(Measure), scales = "free_y", labeller = labeller(Measure = box_labels))


pdf("./plots/countsize_box_all.pdf", w = 15, h = 9)
ggplot(bsc[Measure == "Count" | Measure == "Size"], aes(Treatment, value)) +
  geom_boxplot() +
  geom_point(alpha = 0.3, position = position_jitter(seed = 1, width = 0.2)) +
  theme(panel.grid = element_blank()) +
  labs(y="") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=22, face = "bold"),
        strip.text = element_text(size = 13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(labels = c("No light", "White", "Amber")) +
  #stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
  facet_grid(rows = vars(Measure), scales = "free_y", labeller = labeller(Measure = box_labels))
dev.off()

pdf("./plots/biomass_box_all.pdf", w = 15, h = 9)
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

pdf("./plots/count_box_all.pdf", w = 15, h = 9)
ggplot(capture_size, aes(Treatment,Count)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),
        title = element_text(size=25, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "More arthropods captured under white light than in darkness") +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)
dev.off()

########### POOLED COUNT MODEL ###########
lme_count <- glmmTMB(Count ~ Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                      data = wide_counts[Phylum == "Arthropoda" & is.na(Class)],
                     family = "nbinom2")
plot(simulateResiduals(lme_count))
#r2_nakagawa(lme_count)
summary(lme_count)
emmeans(lme_count, list(pairwise ~ Treatment), adjust = "tukey")

########### POOLED SIZE MODEL ###########
lme_size <- lmer(sqrt(Size) ~ Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance), data = capture_size)
summary(lme_size)
plot(simulateResiduals(lme_size))
emmeans(lme_size, list(pairwise ~ Treatment), adjust = "tukey")
anova(lme_size)
r2_nakagawa(lme_size)


########### ANTS ###########
pdf("./plots/ants.pdf", w = 15, h = 9)
ggplot(wide_counts[Family == "Formicidae"], aes(y = Count, x = Treatment)) +
  geom_col() +
#  scale_fill_brewer(palette = "YlGnBu") +
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        panel.grid = element_blank()) +
  ylab("Number of ants caught") +
  scale_x_discrete(labels = c("No light", "White", "Amber"))
  # stat_compare_means(comparisons = list(c("C","W"), c("C","A"), c("W","A"))) +
  # stat_compare_means(method = "wilcox.test", label.y = 14) +

ggplot(wide_counts[Family == "Carabidae"], aes(y = Count, x = Treatment)) +
  geom_col() +
  #  scale_fill_brewer(palette = "YlGnBu") +
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        panel.grid = element_blank()) +
  ylab("Number of ants caught") +
  scale_x_discrete(labels = c("No light", "White", "Amber"))

########### MULTITAXON MODEL ###########
# for all combinations of Phylum, Class, Order, Family, Subfamily, Genus, Species, Site, Distance, StartDate, and StartDateJulian
# , create a row with a zero count
# wide_counts[,.SD, .SDcols = !c("StartDate")]
# wide_counts[, 
#   CJ(,SD$Phylum,
#      .SD$Class,
#      .SD$Order,
#      .SD$Family,
# #     Subfamily = Subfamily,
# #     Genus = Genus,
# #     Species = Species,
#      .SD$Site),
# #     Distance = Distance,
# #     StartDateJulian = StartDateJulian,
# #     unique = TRUE), 
#   .(StartDateJulian)
# #  on = .(Phylum, Class, Order, Family, Subfamily, Genus, Species, Site, Distance, StartDate, StartDateJulian)
#   ]
# 
# 
# # df <- data.frame(
# #   a = as.factor(c("D", "E", "D", "F", "D", "D")),
# #   x = as.factor(c("A", "A", "A", "B", "B", "C")),
# #   y = as.factor(c("i", "ii", "ii", "i", "i", "i")),
# #   z = 1:6
# # )
# lme_class_1 <- glmmTMB(Count ~ Class + Treatment : Class + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
#                        data = classes,
#                        ziformula= ~1,
#                        family=poisson)
# df <- data.frame(
#   Phylum = as.factor(c("Arthropod", "Arthropod", "Arthropod", "Arthropod", "Arthropod", "Arthropod")),
#   Class = as.factor(c("Insecta", "Insecta", "Insecta", "Arachnida", "Arachnida", NA)),
#   Order = as.factor(c("Coleoptera", "Coleoptera", "Hymenoptera", "Aranaeae", NA, NA)),
#   Site = c(1, 1, 2, 5, 6, 1)
# )
# 
# #nullgrid_expand_dt <- function(df){
#   tf_list <- rep(list(c(TRUE, FALSE)), ncol(df))
#   tmp     <- as.matrix(do.call(CJ, c(tf_list, unique=TRUE))[-.N])
#   OUT     <- lapply(
#     seq_len(nrow(df)),
#     \(i) {
#       x      <- df[rep(i, nrow(tmp)), ]
#       x[tmp] <- NA_character_
#       x
#     }
#   ) |> rbindlist() |> unique()
#   setorderv(OUT, names(OUT))
#   OUT
# #}
# 
# #nullgrid_expand_dt(in_df)
# 
# X = data.table(X.raw)[
#   CJ(y = y,
#      x = x,
#      unique = TRUE), 
#   on = .(x, y)
# ][ , .(z = sum(z)), .(x, y) ][ order(x, y) ]
# X
# 



## CLASSES
zero_classes <- wide_counts[!is.na(Class) & is.na(Order), .SD, .SDcols = c("Class", "Site", "Treatment", "Distance", "StartDateJulian", "Count")]

# add zero rows
zero_classes <- zero_classes[
  CJ(Class = Class,
     Site = Site,
     Treatment = Treatment,
     Distance = Distance,
     StartDateJulian = StartDateJulian,
     unique = TRUE), 
  on = .(Class, Site, Treatment, Distance, StartDateJulian)
]
setnafill(zero_classes, fill = 0, cols=c("Count"))
zero_classes


classes <- zero_classes[, sum(Count), by=Class][zero_classes, on=c("Class")][V1>0.03*sum(unique(V1))][,.SD, .SDcols = !c('V1')]
factor_cols <- c("Class", "Distance", "Site")
classes[, (factor_cols) := lapply(.SD, factor, ordered = FALSE), .SDcols = factor_cols]

lme_class_1 <- glmmTMB(Count ~ Class + Treatment : Class + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = classes,
                     ziformula= ~1,
                     family="nbinom2")
plot(simulateResiduals(lme_class_1))
testZeroInflation(simulateResiduals(lme_class_1))
testDispersion(simulateResiduals(lme_class_1), type = "DHARMa")
testOutliers(simulateResiduals(lme_class_1))

emmeans(lme_class_1, ~Treatment * Class, type="response")
plot(emmeans(lme_class_1, ~Treatment | Class, type="response"), comparisons = TRUE)
pairs(emmeans(lme_class_1, ~Class * Treatment, type="response"), simple = "each", adjust = "holm")






lme_class <- glmmTMB(Count ~ Class * Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = classes,
                     ziformula=~1,
                     family="nbinom2")
plot(simulateResiduals(lme_class))
summary(lme_class)
plot(emmeans(lme_class_1, ~Treatment | Class, type="response"), comparisons = TRUE)

pairs(emmeans(lme_class, ~Class * Treatment, type="response"), simple = "each", adjust = "holm")
emms <-  emmeans(lme_class, ~ Treatment * Distance | Class)
contrast(emms, interaction = "pairwise")



## FAMILIES
families <- wide_counts[!is.na(Family) & is.na(Subfamily), .SD, .SDcols = c("Family", "Site", "Treatment", "Distance", "StartDateJulian", "Count")]
families <- families[
  CJ(Family = Family,
     Site = Site,
     Treatment = Treatment,
     Distance = Distance,
     StartDateJulian = StartDateJulian,
     unique = TRUE), 
  on = .(Family, Site, Treatment, Distance, StartDateJulian)
]
setnafill(families, fill = 0, cols=c("Count"))
families


families_freq <- families[families[, sum(Count), by=Family][V1>0.01*sum(unique(V1))], on=c("Family")]

lme_fam_1 <- glmmTMB(Count ~ Family + Treatment : Family + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = families_freq,
                     ziformula=~1,
                     family=nbinom2)

plot(simulateResiduals(lme_fam_1))
summary(lme_fam_1)

plot(emmeans(lme_fam_1, ~Treatment | Family, type="response"), comparisons = TRUE)
pairs(emmeans(lme_fam_1, ~Family * Treatment, type="response"), simple = "each", adjust = "holm")




lme_fam_dist <- glmmTMB(Count ~ Family * Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                         data = families_freq,
                         ziformula=~1,
                         family=nbinom2)
summary(lme_fam_dist)
plot(emmeans(lme_fam_dist, ~Treatment | Family, type="response"), comparisons = TRUE)

pairs(emmeans(lme_fam_dist, ~Family * Treatment, type="response"), simple = "each", adjust = "holm")
fam_dist_emms <-  emmeans(lme_fam_dist, ~ Treatment * Distance | Family)
contrast(fam_dist_emms, interaction = "pairwise")


## SUBFAMILIES
subfamilies <- wide_counts[!is.na(Subfamily) & is.na(Genus), .SD, .SDcols = c("Subfamily", "Site", "Treatment", "Distance", "StartDateJulian", "Count")]
subfamilies <- subfamilies[
  CJ(Subfamily = Subfamily,
     Site = Site,
     Treatment = Treatment,
     Distance = Distance,
     StartDateJulian = StartDateJulian,
     unique = TRUE), 
  on = .(Subfamily, Site, Treatment, Distance, StartDateJulian)
]
setnafill(subfamilies, fill = 0, cols=c("Count"))
subfamilies


subfamilies_freq <- subfamilies[subfamilies[, sum(Count), by=Subfamily][V1>0.01*sum(unique(V1))], on=c("Subfamily")]

lme_sub_1 <- glmmTMB(Count ~ Subfamily + Treatment : Subfamily + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = subfamilies_freq,
                     ziformula=~1,
                     family=nbinom2)
diagnose(lme_sub_1)
plot(simulateResiduals(lme_sub_1))
testZeroInflation(simulateResiduals(lme_sub_1))
testDispersion(simulateResiduals(lme_sub_1), type = "DHARMa")
testOutliers(simulateResiduals(lme_sub_1))
emmeans(lme_sub_1, ~Treatment * Subfamily, type="response")
pairs(emmeans(lme_sub_1, ~Subfamily * Treatment, type="response"), simple = "each", adjust = "holm")

plot(emmeans(lme_sub_1, ~Treatment | Subfamily, type="response"), comparisons = TRUE)



## GENERA
genera <- wide_counts[!is.na(Genus) & is.na(Species), .SD, .SDcols = c("Genus", "Site", "Treatment", "Distance", "StartDateJulian", "Count")]
genera <- genera[
  CJ(Genus = Genus,
     Site = Site,
     Treatment = Treatment,
     Distance = Distance,
     StartDateJulian = StartDateJulian,
     unique = TRUE), 
  on = .(Genus, Site, Treatment, Distance, StartDateJulian)
]
setnafill(genera, fill = 0, cols=c("Count"))
genera


genera_freq <- genera[genera[, sum(Count), by=Genus][V1>0.03*sum(unique(V1))], on=c("Genus")]

lme_gen_1 <- glmmTMB(Count ~ Genus + Treatment : Genus + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = genera_freq,
                     ziformula=~1,
                     family=nbinom2)
plot(simulateResiduals(lme_gen_1))
testZeroInflation(simulateResiduals(lme_gen_1))
testDispersion(simulateResiduals(lme_gen_1), type = "DHARMa")
testOutliers(simulateResiduals(lme_gen_1))
emmeans(lme_gen_1, ~Treatment * Genus, type="response")
plot(emmeans(lme_gen_1, ~Treatment | Genus, type="response"), comparisons = TRUE)
pairs(emmeans(lme_gen_1, ~Genus * Treatment, type="response"), simple = "each", adjust = "holm")










## SPECIES
genera <- wide_counts[!is.na(Species), .SD, .SDcols = c("Species", "Site", "Treatment", "Distance", "StartDateJulian", "Count")]
genera <- genera[
  CJ(Species = Species,
     Site = Site,
     Treatment = Treatment,
     Distance = Distance,
     StartDateJulian = StartDateJulian,
     unique = TRUE), 
  on = .(Species, Site, Treatment, Distance, StartDateJulian)
]
setnafill(genera, fill = 0, cols=c("Count"))
genera


genera_freq <- genera[genera[, sum(Count), by=Species][V1>0.03*sum(unique(V1))], on=c("Species")]

lme_gen_1 <- glmmTMB(Count ~ Species + Treatment : Species + Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance),
                     data = genera_freq,
                     ziformula=~1,
                     family=nbinom2)
plot(simulateResiduals(lme_gen_1))
testZeroInflation(simulateResiduals(lme_gen_1))
testDispersion(simulateResiduals(lme_gen_1), type = "DHARMa")
testOutliers(simulateResiduals(lme_gen_1))
emmeans(lme_gen_1, ~Treatment * Species, type="response")
plot(emmeans(lme_gen_1, ~Treatment | Species, type="response"), comparisons = TRUE)
pairs(emmeans(lme_gen_1, ~Species * Treatment, type="response"), simple = "each", adjust = "holm")







## 







# this is a long-to-wide restructure
families <- dcast(families, Site + Treatment + Distance + StartDateJulian ~ Family, value.var = "Count")
families[is.na(families)] <- 0
families.fams <- families[, !c("Site", "Treatment", "Distance", "StartDateJulian")]
families.factors <- families[, c("Site", "Treatment", "Distance", "StartDateJulian")]
#families.factors[, Distance := as.numeric(gsub("m","",Distance))]

#families.factors$Distance <- as.factor(families.factors$Distance)
famSite <- as.factor(families.factors$Site)
famDistance <- as.factor(families.factors$Distance)
famDistanceNum <- as.numeric(gsub("m", "", families.factors$Distance))
famTreatment <- families.factors$Treatment
famStartDateJulian <- families.factors$StartDateJulian

fam.abund <- mvabund(families.fams)
fam.abund.nb <- manyglm(fam.abund ~ famTreatment * famDistance +
                          famStartDateJulian,
                        family = "negative.binomial")
summary(fam.abund.nb, p.uni = "adjusted")




fam.abund.nb <- manyglm(fam.abund ~ famTreatment * famDistance +
                          famStartDateJulian + (1|famSite)
                        + (1|famSite:famDistance),
                        family = "negative.binomial")
fam.abund.nb <- manyglm(fam.abund ~ families.factors$Treatment * families.factors$Distance +
                                     families.factors$StartDateJulian + (1|families.factors$Site)
                                     + (1|families.factors$Site:families.factors$Distance),
                        family = "negative.binomial")

families.abund.nb <- glmer.nb(Aphididae ~ Treatment * Distance +
                          StartDateJulian + (1|Site)
                          + (1|Site:Distance),
                          data = families,
                          family = "negative.binomial")

carabidae.nb <- glmer.nb(Carabidae ~ Treatment * Distance +
                                StartDateJulian + (1|Site)
                              + (1|Site:Distance),
                              data = families,
                              family = "negative.binomial")

data.frame(Site = families$Site, Treatment = families$Treatment, 
           Distance = families$Distance, StartDateJulian = families$StartDateJulian)

try_data <- copy(families[, c("Site", "Treatment", "Distance", "StartDateJulian", "Carabidae")])

try_data[, Predictions := fitted(carabidae.nb)]










# fam.abund.nb <- manyglm(fam.abund ~ families.factors$Treatment * families.factors$Distance + 
#                                     families.factors$StartDateJulian + (1|families.factors$Site)
#                                     + (1|families.factors$Site:families.factors$Distance),
#                        family = "negative.binomial")
#, p.uni = "adjusted")





ggplot(biomass, aes(Treatment, Near)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),
        title = element_text(size=32, face = "bold")) +
  labs(y = "Biomass (grams)", title = "Biomass captured at pitfall traps 0m from lights") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter() +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 8)



ggboxplot(capture_size, x = "Treatment", y = "Size") +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(y = "Biomass captured divided by number of arthropods") +
  geom_jitter()

lmec_size <- glmer.nb(Size ~ Treatment * Distance + StartDateJulian + (1|Site) + (1|Site:Distance), data = capture_size)
summary(lmec_size)






# all_pitfall_names[, id := .I]
# all_pitfall_names <- all_pitfall_names[StartDate > as.Date("2023-06-21")]

######################## STACKED BAR PLOT ########################
stacked_bar_prep <- function(pitfall_taxa) {
  #pitfall_taxa[, .N, by=.(StartDate, Site, Distance, Treatment, Family)]
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDate", "Distance", "Treatment")])
  
  pitfall_taxa <- reshape(pitfall_taxa[, .N, by=names(pitfall_taxa)], 
                          idvar = c("Site", "StartDate", "Distance", "Treatment"), 
                          timevar = taxon, direction = "wide")
  
  setnafill(pitfall_taxa, cols = names(pitfall_taxa)[grepl("N.", names(pitfall_taxa))], fill = 0)
  
  pitfall_taxa_n <- pitfall_taxa[, lapply(.SD, sum), by=.(Treatment), .SDcols = !c("StartDate", "Distance", "Site", "Treatment")]
  pitfall_taxa_n <- transpose(pitfall_taxa_n, keep.names = "Treatment", make.names = "Treatment")
  
  pitfall_taxa_n[, (names(pitfall_taxa_n[, .SD, .SDcols = !c("Treatment")])) := lapply(.SD, function(x) x/sum(x)), .SDcols = !c("Treatment")]
  setnames(pitfall_taxa_n, "Treatment", taxon)
  pitfall_taxa_n <- melt(pitfall_taxa_n, id.vars = taxon, variable.name = "Treatment")
  setnames(pitfall_taxa_n, "value", "Percent")
  pitfall_taxa_n$Treatment <- factor(pitfall_taxa_n$Treatment, levels=c("C","W","A"), ordered = FALSE)
  return(
    ggplot(pitfall_taxa_n, aes(fill = get(taxon), y = Percent, x = Treatment)) +
      geom_bar(position = "stack", stat = "identity") +
      geom_text(aes(label=scales::percent(Percent)), position = position_stack(vjust = .5)) +
      scale_y_continuous(labels = scales::percent) +
      labs(fill=taxon) +
      scale_fill_brewer(palette = "YlGnBu"))
}

wide_counts[, .SD, .SDcols = ]
wide_counts[!is.na(Class) & is.na(Order), .(Class, StartDate, Site, Distance, Treatment, Count)]


stacked_bar_prep <- function(pitfall_taxa) {
  #pitfall_taxa[, .N, by=.(StartDate, Site, Distance, Treatment, Family)]
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDate", "Distance", "Treatment")])
  
  pitfall_taxa <- reshape(pitfall_taxa[, .N, by=names(pitfall_taxa)], 
                          idvar = c("Site", "StartDate", "Distance", "Treatment"), 
                          timevar = taxon, direction = "wide")
  
  setnafill(pitfall_taxa, cols = names(pitfall_taxa)[grepl("N.", names(pitfall_taxa))], fill = 0)
  
  pitfall_taxa_n <- pitfall_taxa[, lapply(.SD, sum), by=.(Treatment), .SDcols = !c("StartDate", "Distance", "Site", "Treatment")]
  pitfall_taxa_n <- transpose(pitfall_taxa_n, keep.names = "Treatment", make.names = "Treatment")
  
  pitfall_taxa_n[, (names(pitfall_taxa_n[, .SD, .SDcols = !c("Treatment")])) := lapply(.SD, function(x) x/sum(x)), .SDcols = !c("Treatment")]
  setnames(pitfall_taxa_n, "Treatment", taxon)
  pitfall_taxa_n <- melt(pitfall_taxa_n, id.vars = taxon, variable.name = "Treatment")
  setnames(pitfall_taxa_n, "value", "Percent")
  pitfall_taxa_n$Treatment <- factor(pitfall_taxa_n$Treatment, levels=c("C","W","A"), ordered = FALSE)
  return(
    ggplot(pitfall_taxa_n, aes(fill = get(taxon), y = Percent, x = Treatment)) +
      geom_bar(position = "stack", stat = "identity") +
      geom_text(aes(label=scales::percent(Percent)), position = position_stack(vjust = .5)) +
      scale_y_continuous(labels = scales::percent) +
      labs(fill=taxon) +
      scale_fill_brewer(palette = "YlGnBu"))
}

### below is old version (percent of catch per treatment)
pitfall_class <- copy(wide_counts[!is.na(Class) & is.na(Order), .(Class, StartDate, Site, Distance, Treatment)])
stacked_bar_prep(pitfall_class) + labs(title = "Proportion of observed taxonomic classes by treatment")


### below divides by total number of traps *where creature was present*. need to divide by total number of traps
pitfall_classes <- copy(wide_counts[!is.na(Class) & is.na(Order) & Distance == "0m", .(Class, StartDate, Site, Treatment, Count)])
taxon <- names(pitfall_classes[, .SD, .SDcols = !c("Site", "StartDate", "Treatment", "Count")])
pitfall_classes <- pitfall_classes[, .(Traps = .N, Number = sum(Count), PerTrap = sum(Count)/.N), by=.(get(taxon), Treatment)]
setnames(pitfall_classes,"get", taxon)


### Treatment NTraps      NCaptured     NCapturedPerTrap
### C         212         150           
### W         110         80
### A         107         70


pitfall_classes <- copy(wide_counts[!is.na(Class) & is.na(Order) & Distance == "0m", .(Class, StartDate, Site, Treatment, Count)])

stacked_bar_per_trap <- function(pitfall_taxa) {
  pitfall_taxa[pitfall_taxa[, .N, by=Treatment], TotalTraps := i.N, on = c("Treatment")]
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDate", "Treatment", "Count", "TotalTraps")])
  pitfall_taxa <- unique(pitfall_taxa[, .(Number = sum(Count), PerTrap = sum(Count)/TotalTraps), by=.(get(taxon), Treatment)])
  setnames(pitfall_taxa,"get", taxon)
  return(ggplot(pitfall_taxa[PerTrap > 0.1], aes(fill = get(taxon), y = PerTrap, x = Treatment)) +
           geom_bar(position = "stack", stat = "identity") +
           geom_text(aes(label=round(PerTrap, digits = 2)), position = position_stack(vjust = .5)) +
           labs(fill=taxon) +
           scale_fill_brewer(palette = "YlGnBu")
  )
}

stacked_bar_per_trap(wide_counts[!is.na(Class) & is.na(Order) & Distance == "0m", .(Class, StartDate, Site, Treatment, Count)])
stacked_bar_per_trap(wide_counts[!is.na(Order) & is.na(Family) & Distance == "0m", .(Order, StartDate, Site, Treatment, Count)])
stacked_bar_per_trap(wide_counts[!is.na(Family) & is.na(Subfamily) & Distance == "0m", .(Family, StartDate, Site, Treatment, Count)])
stacked_bar_per_trap(wide_counts[!is.na(Subfamily) & is.na(Genus) & Distance == "0m", .(Subfamily, StartDate, Site, Treatment, Count)])
stacked_bar_per_trap(wide_counts[!is.na(Genus) & is.na(Species) & Distance == "0m", .(Genus, StartDate, Site, Treatment, Count)])



faceted_bar_per_trap <- function(pitfall_taxa) {
  pitfall_taxa[pitfall_taxa[, .N, by=Treatment], TotalTraps := i.N, on = c("Treatment")]
  
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDate", "Treatment", "Count", "TotalTraps")])
  pitfall_taxa <- unique(pitfall_taxa[, .(Number = sum(Count), PerTrap = sum(Count)/TotalTraps), by=.(get(taxon), Treatment)])
  setnames(pitfall_taxa,"get", taxon)
  return(ggplot(pitfall_taxa[PerTrap > 0.1], aes(fill = get(taxon), y = PerTrap, x = Treatment)) +
           geom_bar(stat = "identity") +
           geom_text(aes(label=round(PerTrap, digits = 2))) +
           labs(fill=taxon) +
           scale_fill_brewer(palette = "YlGnBu") +
           theme(title = element_text(size = 20),
                 axis.text = element_text(size = 12),
                 legend.text = element_text(size = 14),
                 strip.text.y = element_text(size = 14)) +
           # stat_compare_means(comparisons = list(c("C","W"), c("C","A"), c("W","A"))) +
           # stat_compare_means(method = "wilcox.test", label.y = 14) +
           facet_grid(rows = vars(get(taxon)))
  )
}

faceted_bar_per_trap(wide_counts[!is.na(Class) & is.na(Order) & Distance == "0m", .(Class, StartDate, Site, Treatment, Count)])
faceted_bar_per_trap(wide_counts[!is.na(Order) & is.na(Family) & Distance == "0m", .(Order, StartDate, Site, Treatment, Count)])
faceted_bar_per_trap(wide_counts[!is.na(Family) & is.na(Subfamily) & Distance == "0m", .(Family, StartDate, Site, Treatment, Count)])
faceted_bar_per_trap(wide_counts[!is.na(Subfamily) & is.na(Genus) & Distance == "0m", .(Subfamily, StartDate, Site, Treatment, Count)])
# staphylininae have a BUNCH of diversity in how they live...carabids on the other hand are mostly visual predators
faceted_bar_per_trap(wide_counts[!is.na(Genus) & is.na(Species) & Distance == "0m", .(Genus, StartDate, Site, Treatment, Count)])







faceted_bar_common <- function(pitfall_taxa) {
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDate", "Treatment", "Count")])
  pitfall_taxa <- pitfall_taxa[, sum(Count), by=taxon][pitfall_taxa, on=c(taxon)][V1>0.01*sum(unique(V1))][,.SD, .SDcols = !c('V1')]
  pitfall_taxa[pitfall_taxa[, .N, by=Treatment], TotalTraps := i.N, on = c("Treatment")]
  
  pitfall_taxa <- unique(pitfall_taxa[, .(Number = sum(Count), PerTrap = sum(Count)/TotalTraps), by=.(get(taxon), Treatment)])
  setnames(pitfall_taxa,"get", taxon)
  return(
  ggplot(pitfall_taxa[PerTrap > 0.1], aes(fill = get(taxon), y = PerTrap, x = Treatment)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=round(PerTrap, digits = 2))) +
    labs(fill=taxon) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme(title = element_text(size = 20),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          panel.grid = element_blank()) +
    # stat_compare_means(comparisons = list(c("C","W"), c("C","A"), c("W","A"))) +
    # stat_compare_means(method = "wilcox.test", label.y = 14) +
    facet_grid(rows = vars(get(taxon)), scales = "free_y")
  )
}

faceted_bar_common(wide_counts[!is.na(Class) & is.na(Order) & Distance == "0m", .(Class, StartDate, Site, Treatment, Count)])
faceted_bar_common(wide_counts[!is.na(Order) & is.na(Family) & Distance == "0m", .(Order, StartDate, Site, Treatment, Count)])
faceted_bar_common(wide_counts[!is.na(Family) & is.na(Subfamily) & Distance == "0m", .(Family, StartDate, Site, Treatment, Count)])
faceted_bar_common(wide_counts[Order == "Coleoptera" & !is.na(Family) & is.na(Subfamily), .(Family, StartDate, Site, Treatment, Count)])

#summary(emmeans(lme_fam_1), type="response")
lme_fam_1_summary <- summary(emmeans(lme_fam_1, ~Treatment : Family, type="response"), 
                             comparisons = TRUE)

lme_fam_1_means <- emmeans(object = lme_fam_1,
                             specs = ~ Treatment : Family,
                           type = "response")

lme_fam_1_cld <- cld(object = lme_fam_1_means,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
lme_fam_1_cld
library(stringr)

ggplot(data = lme_fam_1_cld,
       aes(x = Treatment, y = response, fill = Family)) +
  facet_grid(rows = vars(Family), scales = "free_y") +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE), width = 0.2) +
  scale_y_continuous(
    name = "EMMEANS",
    expand = expansion(mult = c(0,0.1)),
                       labels = scales::number_format(accuracy = 0.1)
    ) +
  scale_fill_brewer(palette = "YlGnBu") +
  xlab(NULL) +
  geom_text(aes(label = str_trim(.group), y = response + SE), vjust = -0.5) +
  labs(title = "Families") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(),
    axis.title.x = element_text(vjust = 0),
    axis.title.y = element_text(size = 8)
  ) +
  theme(legend.title = element_blank(),
        legend.position = "none")






faceted_bar_per_trap_carabidae <- function(pitfall_taxa) {
  pitfall_taxa[pitfall_taxa[, .N, by=Treatment], TotalTraps := i.N, on = c("Treatment")]
  taxon <- names(pitfall_taxa[, .SD, .SDcols = !c("Site", "StartDateJulian", "Treatment", "Predictions", "TotalTraps")])
  pitfall_taxa <- unique(pitfall_taxa[, .(Number = sum(Predictions), PerTrap = sum(Predictions)/TotalTraps), by=.(get(taxon), Treatment)])
  setnames(pitfall_taxa,"get", taxon)
  return(ggplot(pitfall_taxa[PerTrap > 0.1], aes(fill = get(taxon), y = PerTrap, x = Treatment)) +
           geom_bar(stat = "identity") +
           geom_text(aes(label=round(PerTrap, digits = 2))) +
           labs(fill=taxon) +
           stat_compare_means(comparisons = list(c("C","W"), c("C","A"), c("W","A"))) +
           stat_compare_means(method = "anova", label.y = 14, )
  )
}

faceted_bar_per_trap_carabidae(unique(wide_counts[try_data, on = c("Site", "Treatment", "Distance", "StartDateJulian")][Distance == "0m" & Family == "Carabidae", .(Family, StartDateJulian, Site, Treatment, Predictions)]))

# look at seasonal component for some groups as well


###################


### get total number of traps experiencing each treatment (converted to function stacked_bar_per_trap)
# pitfall_classes[pitfall_classes[, .N, by=Treatment], TotalTraps := i.N, on = c("Treatment")]
# pitfall_classes <- unique(pitfall_classes[, .(Number = sum(Count), PerTrap = sum(Count)/TotalTraps), by=.(get(taxon), Treatment)])
# setnames(pitfall_classes,"get", taxon)
# ggplot(pitfall_classes, aes(fill = get(taxon), y = PerTrap, x = Treatment)) +
#   geom_bar(position = "stack", stat = "identity") +
#   geom_text(aes(label=round(PerTrap, digits = 2)), position = position_stack(vjust = .5)) +
#   labs(fill=taxon) +
#   scale_fill_brewer(palette = "YlGnBu")




pitfall_fams <- copy(wide_counts[!is.na(Family) & is.na(Subfamily), .(Family, StartDate, Site, Distance, Treatment)])




pitfall_ords <- copy(wide_counts[Class == "Insecta" & !is.na(Order) & is.na(Family), .(Order, StartDate, Site, Distance, Treatment)])
stacked_bar_prep(pitfall_ords) + labs(title = paste0("Proportion of observed insect orders by treatment (n = ", dim(pitfall_fams)[1], ")"))

stacked_bar_prep(pitfall_fams) + labs(title = paste0("Proportion of observed taxonomic families by treatment (n = ",dim(pitfall_fams)[1],")"))

pitfall_genus <- copy(wide_counts[Family == "Carabidae" & !is.na(Genus) & is.na(Species), .(Genus, StartDate, Site, Distance, Treatment)])
stacked_bar_prep(pitfall_genus) + labs(title = paste0("Proportion of observed beetle genera by treatment (n = ",dim(pitfall_genus)[1],")"))


########### HISTOGRAM ############

filtered_pitfall <- metatable[pitfall_names, on = .(StartDate, Site), nomatch = 0]
filtered_pitfall[metatable, Treatment := i.Treatment, on=.(Site, StartDate)]

# divide the number of observations by the number of sites for each StartDate

obs_col_plot <- function(data) {
  count_bar <- as.data.table(data[, .N, by=.(Site, StartDate)][, .N, by=StartDate][order(StartDate)])
  setnames(count_bar, "N", "NumberOfSites")
  count_bar
  count_bar[, TotalObservations := data[, .N, by=StartDate][order(StartDate)]$N]
  count_bar[, Normalised := TotalObservations/NumberOfSites]
  return(ggplot(count_bar, aes(x=StartDate, y=Normalised)) +
    geom_line())
}

library(gridExtra)

all_col <- obs_col_plot(filtered_pitfall)
all_col + labs(y = "Total count divided by number of sites")

c_col <- obs_col_plot(filtered_pitfall[Treatment == "C"])
c_col + labs(y = "Total count under control (darkness) divided by number of sites")

a_col <- obs_col_plot(filtered_pitfall[Treatment == "A"])
a_col + labs(y = "Total count under amber light divided by number of sites")

w_col <- obs_col_plot(filtered_pitfall[Treatment == "W"])
w_col + labs(y = "Total count under white light divided by number of sites")

col_list <- c(all_col, c_col, a_col, w_col)

count_bar <- as.data.table(filtered_pitfall[, .N, by=.(Site, StartDate)][, .N, by=StartDate][order(StartDate)])
setnames(count_bar, "N", "NumberOfSites")
ggplot(count_bar, aes(x=StartDate, y=NumberOfSites)) +
         geom_col()

mcount_bar <- as.data.table(metatable[, .N, by=.(Site, StartDate)][, .N, by=StartDate][order(StartDate)])
setnames(mcount_bar, "N", "NumberOfSites")
ggplot(mcount_bar, aes(x=StartDate, y=NumberOfSites)) +
  geom_point()



filtered_pitfall[, .N, by=.(StartDate, Site, Treatment)][, .N, by=.(Treatment)]
filtered_pitfall[Treatment == "", .N, by=.(StartDate, Site)][order(StartDate)][order(StartDate, Site)]





########### INDEX ############

counts[, Index := Near/(Near+Far), by = .(StartDate,Site)]

ggplot(counts, aes(Treatment, Index)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Count statistic") +
  geom_jitter() +
  geom_hline(aes(yintercept = 0.5),colour="red",linetype="dotted")

wilcoxon_count <- data.frame(W = wilcox.test(counts[Treatment == ("W"),Index], counts[Treatment == ("C"),Index], 
            alternative = "two.sided", exact = FALSE)$p.value,
            A = wilcox.test(counts[Treatment == ("A"),Index], counts[Treatment == ("C"),Index], 
            alternative = "two.sided", exact = FALSE)$p.value)
wilcoxon_count



############# N_ISOPODA ######################

counts <- dcast(raw_counts, StartDate+Site~Distance, value.var=c("N_Isopoda"))
setnames(counts, c("0m","5m"), c("Near","Far"))
counts[metatable, Treatment := i.Treatment, on=.(Site, StartDate)]

counts[, Index := Near/(Near+Far), by = .(StartDate,Site)]

counts <- counts[Treatment %in% c("C","A","W")]

counts$Treatment <- factor(counts$Treatment, levels=c("C","W","A"), ordered = TRUE)

ggplot(counts, aes(Treatment, Index)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Count statistic") +
  geom_jitter() +
  geom_hline(aes(yintercept = 0.5),colour="red",linetype="dotted")

wilcoxon_count <- data.frame(W = wilcox.test(counts[Treatment == ("W"),Index], counts[Treatment == ("C"),Index], 
                                             alternative = "two.sided", exact = FALSE)$p.value,
                             A = wilcox.test(counts[Treatment == ("A"),Index], counts[Treatment == ("C"),Index], 
                                             alternative = "two.sided", exact = FALSE)$p.value)
wilcoxon_count




############# N_ISOPODA ######################

counts <- dcast(raw_counts, StartDate+Site~Distance, value.var=c("N_Isopoda"))
setnames(counts, c("0m","5m"), c("Near","Far"))
counts[metatable, Treatment := i.Treatment, on=.(Site, StartDate)]

counts[, Index := Near/(Near+Far), by = .(StartDate,Site)]

counts <- counts[Treatment %in% c("C","A","W")]

counts$Treatment <- factor(counts$Treatment, levels=c("C","W","A"), ordered = TRUE)

ggplot(counts, aes(Treatment, Index)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Count statistic") +
  geom_jitter() +
  geom_hline(aes(yintercept = 0.5),colour="red",linetype="dotted")

wilcoxon_count <- data.frame(W = wilcox.test(counts[Treatment == ("W"),Index], counts[Treatment == ("C"),Index], 
                                             alternative = "two.sided", exact = FALSE)$p.value,
                             A = wilcox.test(counts[Treatment == ("A"),Index], counts[Treatment == ("C"),Index], 
                                             alternative = "two.sided", exact = FALSE)$p.value)
wilcoxon_count


######################## HIERARCHICAL VISUALISATION ######################## 

library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)
library(gridExtra)

prepare_graph <- function(pitfall_subset) {
  
  count_pn <- pitfall_subset[, .N, by=.(Treatment, Root,Phylum,Class,Order,Family,Subfamily,Genus,Species)]
  
  
  #count_pn <- pitfall_names[Treatment == "C"][, .N, by=.(Treatment, Root,Phylum,Class,Order,Family,Subfamily,Genus,Species)]
  
  count_edges <- na.omit(rbind(
    data.frame(from=count_pn$Root, to=count_pn$Phylum),
    data.frame(from=count_pn$Phylum, to=count_pn$Class),
    data.frame(from=count_pn$Class, to=count_pn$Order),
    data.frame(from=count_pn$Order, to=count_pn$Family),
    data.frame(from=count_pn$Family, to=count_pn$Subfamily),
    data.frame(from=count_pn$Subfamily, to=count_pn$Genus),
    data.frame(from=count_pn$Genus, to=count_pn$Species)))
  count_edges <- unique(count_edges)
  
  count_name <- unique(c(as.character(count_edges$from),as.character(count_edges$to)))
  
  count_pn[, Depth := Reduce(`+`, lapply(.SD,function(x) !is.na(x))),.SDcols=taxa_cols]
  
  # get MSTL
  count_pn[, Name := NA_character_]
  for (v in rev(names(count_pn[,..taxa_cols]))) 
    count_pn[is.na(Name), Name := get(v)]
  
  setdiff(count_pn$Name,unique(c(count_edges$from,count_edges$to)))
  names_to_add <- setdiff(unique(c(count_edges$from,count_edges$to)),count_pn$Name)
  
  n_to_add <- rep(0,length(names_to_add))
  
  # add new rows in count_pn for the elements that aren't currently in there
  # get column # in count_pn[,..taxa_cols] for each name
  names_to_add
  # we have to lapply the column-finding across the whole array of names
  
  get_depth <- function(name, data) {
    return(data[, which(.SD==name, arr.ind = TRUE)][,2][[1]])
  }
  depths_to_add <- unlist(lapply(names_to_add,get_depth,data=count_pn[,..taxa_cols]))
  depths_to_add
  
  to_add <- data.frame(name = names_to_add,
                       value = n_to_add,
                       depth = depths_to_add)
  
  c_vertices <- data.frame(
    name=count_pn$Name,
    value=count_pn$N,
    depth=count_pn$Depth
  )
  
  c_vertices <- rbindlist(list(c_vertices,to_add))
  
  c_vertices$id=NA
  c_leaves=which(is.na(match(c_vertices$name, count_edges$from)))
  c_nleaves = length(c_leaves)
  c_vertices$id[c_leaves] = seq(1:c_nleaves)
  c_vertices$angle= 90 - 360 * c_vertices$id / c_nleaves
  c_vertices$hjust<-ifelse( c_vertices$angle < -90, 1, 0)
  c_vertices$angle<-ifelse(c_vertices$angle < -90, c_vertices$angle+180, c_vertices$angle)
  
  setdiff(c_vertices$name,unique(c(count_edges$from,count_edges$to)))
  setdiff(unique(c(count_edges$from,count_edges$to)),c_vertices$name)
  
  my_graph <- graph_from_data_frame(count_edges,vertices = c_vertices)
  
  return(my_graph)
}

set.seed(1)
control_graph <- prepare_graph(pitfall_names[Treatment == "C"])
set.seed(1)
white_graph <- prepare_graph(pitfall_names[Treatment == "W"])
set.seed(1)
amber_graph <- prepare_graph(pitfall_names[Treatment == "A"])



plot_list <- list(ggraph(control_graph, layout = 'circlepack', weight=value) + 
                    geom_node_circle(aes(fill=depth)) +
                    geom_node_label(aes(label=name,size=value),position=position_jitter(width=.2,height=.2),repel=TRUE,max.overlaps=10) +
                    theme_void() +
                    scale_fill_viridis(direction=1),
                  ggraph(white_graph, layout = 'circlepack', weight=value) + 
                    geom_node_circle(aes(fill=depth)) +
                    geom_node_label(aes(label=name,size=value),position=position_jitter(width=.2,height=.2),repel=TRUE,max.overlaps=10) +
                    theme_void() +
                    scale_fill_viridis(direction=1),
                  ggraph(amber_graph, layout = 'circlepack', weight=value) + 
                    geom_node_circle(aes(fill=depth)) +
                    geom_node_label(aes(label=name,size=value),position=position_jitter(width=.2,height=.2),repel=TRUE,max.overlaps=10) +
                    theme_void() +
                    scale_fill_viridis(direction=1))

# pdf("plots.pdf", onefile = TRUE)
# for (i in seq(length(p))) {
#   do.call("grid.arrange", p[[i]])  
# }
# dev.off()

pdf("./plots/hierarch_vis_by_treatments.pdf", onefile = TRUE, w=12, h=12)
par(mfrow=c(1,3))
#layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
ggraph(control_graph, layout = 'circlepack', weight=value) + 
  geom_node_circle(aes(fill=depth)) +
  geom_node_label(aes(label=name,size=value),repel=TRUE,max.overlaps=10) +
  theme_void() +
  scale_fill_viridis(direction=1)
ggraph(white_graph, layout = 'circlepack', weight=value) + 
  geom_node_circle(aes(fill=depth)) +
  geom_node_label(aes(label=name,size=value),repel=TRUE,max.overlaps=10) +
  theme_void() +
  scale_fill_viridis(direction=1)
ggraph(amber_graph, layout = 'circlepack', weight=value) + 
  geom_node_circle(aes(fill=depth)) +
  geom_node_label(aes(label=name,size=value),repel=TRUE,max.overlaps=10) +
  theme_void() +
  scale_fill_viridis(direction=1)
dev.off()
p 

###################### TAXONOMIC DISTANCE###################### 
taxon_table <- unique(copy(pitfall_names[, ..taxa_cols]))

taxdis <- taxa2dist(taxon_table, varstep = FALSE, check = TRUE)

plot(hclust(taxdis), hang = -1)





test_data_table <- data.table(StartDate = rep(seq(1:5),each = 2),        # Create first data.table
                   Site = rep(seq(1:5),2),
                   FullName = letters[1:10],
                   Family = c("Apple", "Apple", "Banana", "Apple", "Banana", "Apple", "Banana", "Banana", "Banana", "Banana"))
reshape(test_data_table[, .N, by=.(StartDate, Site, Family)], idvar = c("StartDate", "Site"), timevar = "Family", direction = "wide")


# setnames(name_cols, "V1", "FullName")
# name_cols[, Parent := sub("\\/[A-Za-z0-9?._-]+$", "", FullName)]
# name_cols <- subset(name_cols, Parent != FullName)
# name_graph <- graph.data.frame(as.data.frame(name_cols))
# 
# node <- "./Arthropoda"
# setdiff(names(subcomponent(name_graph, node, mode = "out")), node)
# 
# setdiff(names(subcomponent(name_graph, "./Arthropoda", mode = "out")), "./Arthropoda")

#, c_vertices = as.data.frame(pitfall_names)





########### BOXPLOT 0m ############

ggboxplot(wide_counts[is.na(Class) & Distance == "0m"], x = "Treatment", y = "Count") +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(y = "Number of arthropods captured") +
  geom_jitter() +
  geom_hline(aes(yintercept = 0.5),colour="red",linetype="dotted") +
  stat_compare_means(method = "anova", label.y = 14) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")


ggboxplot(counts, x = "Treatment", y = "Near") +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(y = "Number of trapped arthropods") +
  geom_jitter() +
  geom_hline(aes(yintercept = 0.5),colour="red",linetype="dotted") +
  stat_compare_means(method = "anova", label.y = 14) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")

wilcoxon_count <- data.frame(W = wilcox.test(counts[Treatment == ("W"),Near], counts[Treatment == ("C"),Near], 
                                             alternative = "two.sided", exact = FALSE)$p.value,
                             A = wilcox.test(counts[Treatment == ("A"),Near], counts[Treatment == ("C"),Near], 
                                             alternative = "two.sided", exact = FALSE)$p.value)
wilcoxon_count

########### BOXPLOT 5m ############


ggplot(counts, aes(Treatment, Far)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Count statistic")

wilcoxon_count <- data.frame(W = wilcox.test(counts[Treatment == ("W"),Far], counts[Treatment == ("C"),Far], 
                                             alternative = "two.sided", exact = FALSE)$p.value,
                             A = wilcox.test(counts[Treatment == ("A"),Far], counts[Treatment == ("C"),Far], 
                                             alternative = "two.sided", exact = FALSE)$p.value)
wilcoxon_count




