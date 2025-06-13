## This R script describes the code used for Graphical and Statistical analysis associated with the manuscript: ##
## Cerebellum mitochondrial DNA copy number is increased in Parkinson’s Disease ##
## published in Brain Communications - 2025 ##
## author: Brooke Hjelm, PhD (bhjelm@usc.edu) ##

## Note: Data provided in Extended Data Tables 1 and 2 with manuscript do not contain exact age for subjects > 90yo ##
## Age was not significant but is tested in code and included as co-variate in rfit/rlm models ##
## These commands will not result in exact same p-values in manuscript due to missing age values ##
## Code where age is removed as co-variate from rfit models is provided as well ##

library(ggplot2)
library(Rfit)
library(FSA)
library(dplyr)
library(ggcorrplot)
library(sjPlot)
library(sfsmisc)
library(MASS)

##############################################################################################
## Loading Data for Figure 1 ##
PDvCTRL <- read.table('~/Downloads/Extended_Data_Table1_v3.txt', header=TRUE, sep="\t")

## Log10 transformation of mtDNA copy number ##
PDvCTRL$mt_copy_number_avg_Log10 <- log10(PDvCTRL$mt_copy_number_avg)

## Graph and Stats - Figure 1a - Age ##
## This won't run correctly if using Extended Data Table 1 since some age values are binned as characters ##
AGE_PLOT <- ggplot(PDvCTRL, aes(x = Age, y = mt_copy_number_avg_Log10)) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  geom_smooth(method = 'lm', se = T, colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

AGE_PLOT + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Age", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

cor.test(PDvCTRL$Age, PDvCTRL$mt_copy_number_avg_Log10, alternative = c("two.sided"), method=c("spearman"))

## Graph and Stats - Figure 1b - Sex ##
SEX_PLOT <- ggplot(PDvCTRL, aes(x = SEX, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SEX_PLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Sex", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.sex = rfit(mt_copy_number_avg_Log10 ~ SEX + Age, data = PDvCTRL)
summary(model.sex)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.sex.noAGE = rfit(mt_copy_number_avg_Log10 ~ SEX, data = PDvCTRL)
summary(model.sex.noAGE)

## Graph and Stats - Figure 1c - Autosomal Coverage ##
AUTO_PLOT <- ggplot(PDvCTRL, aes(x = autosomal_coverage, y = mt_copy_number_avg_Log10)) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  geom_smooth(method = 'lm', se = T, colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

AUTO_PLOT + theme(panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                     panel.grid.major = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Autosomal Coverage", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

cor.test(PDvCTRL$autosomal_coverage, PDvCTRL$mt_copy_number_avg_Log10, alternative = c("two.sided"), method=c("spearman"))

## Graph and Stats - Figure 1d - Diagnosis (All Brain Banks) ##
DIAGNOSIS_PLOT <- ggplot(PDvCTRL, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

DIAGNOSIS_PLOT + theme(panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                         panel.grid.major = element_blank(), 
                         axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.dx = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX + Source, data = PDvCTRL)
summary(model.dx)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.dx.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX + Source, data = PDvCTRL)
summary(model.dx.noAGE)

## Separating PD and Control Samples ##
PD <- subset(PDvCTRL, Status == "PD")
CONTROLS <- subset(PDvCTRL, Status == "Control")

## Separating PD Samples from Each Brain Bank & Combining with Controls ##
BANNER <- subset(PD, Source == "BSHRI")
HOPKINS <- subset(PD, Source == "JHU")
SEPULV <- subset(PD, Source == "Sepulveda")
MARYL <- subset(PD, Source == "UofMaryland")
HARVARD <- subset(PD, Source == "Harvard")

BANNER2 <- rbind(BANNER,CONTROLS)
HOPKINS2 <- rbind(HOPKINS,CONTROLS)
SEPULV2 <- rbind(SEPULV,CONTROLS)
MARYL2 <- rbind(MARYL,CONTROLS)
HARVARD2 <- rbind(HARVARD,CONTROLS)

## Graph and Stats - Figure 1e - Diagnosis (Banner) ###
BANNER_PLOT <- ggplot(BANNER2, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

BANNER_PLOT + theme(panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                      panel.grid.major = element_blank(), 
                      axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.bann = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = BANNER2)
summary(model.bann)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.bann.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = BANNER2)
summary(model.bann.noAGE)

## Graph and Stats - Figure 1f - Diagnosis (JHU) ###
JHU_PLOT <- ggplot(HOPKINS2, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

JHU_PLOT + theme(panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                   panel.grid.major = element_blank(), 
                   axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.jhu = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = HOPKINS2)
summary(model.jhu)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.jhu.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = HOPKINS2)
summary(model.jhu.noAGE)

## Graph and Stats - Figure 1g - Diagnosis (Sepulveda) ###
SEPULV_PLOT <- ggplot(SEPULV2, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SEPULV_PLOT + theme(panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                      panel.grid.major = element_blank(), 
                      axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.sepul = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = SEPULV2)
summary(model.sepul)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.sepul.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = SEPULV2)
summary(model.sepul.noAGE)

## Graph and Stats - Figure 1h - Diagnosis (U Maryland) ###
MARY_PLOT <- ggplot(MARYL2, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

MARY_PLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.mary = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = MARYL2)
summary(model.mary)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.mary.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = MARYL2)
summary(model.mary.noAGE)

## Graph and Stats - Figure 1i - Diagnosis (Harvard) ###
HARV_PLOT <- ggplot(HARVARD2, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 2) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

HARV_PLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.harv = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = HARVARD2)
summary(model.harv)
## model without age as co-variate - run if using Extended Data Table 1 as input ##
model.harv.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = HARVARD2)
summary(model.harv.noAGE)

##############################################################################################
## Loading Data for Figure 2 ##
PATH <- read.table('~/Downloads/Extended_Data_Table2_v3-2.txt', header=TRUE, sep="\t")

## Log10 transformation of mtDNA copy number ##
PATH$mt_copy_number_avg_Log10 <- log10(PATH$mt_copy_number_avg)

## Graph and Stats - Figure 2a - Diagnosis (Banner) ###
DIAGNOSIS_BANNER_PATH <- ggplot(PATH, aes(x = Status, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

DIAGNOSIS_BANNER_PATH + theme(panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                         panel.grid.major = element_blank(), 
                         axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Diagnosis", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

model.bannerpath = rfit(mt_copy_number_avg_Log10 ~ Status + Age + SEX, data = PATH)
summary(model.bannerpath)
## model without age as co-variate - run if using Extended Data Table 2 as input ##
model.bannerpath.noAGE = rfit(mt_copy_number_avg_Log10 ~ Status + SEX, data = PATH)
summary(model.bannerpath.noAGE)

## Graph and Stats - Figure 2b - USSLB Stage ##
UNIFIED_LBPLOT <- ggplot(PATH, aes(x = Unified_LB_Stage, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

UNIFIED_LBPLOT + theme(panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                       panel.grid.major = element_blank(), 
                       axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Unified LB Density Stage", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ Unified_LB_Stage, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ Unified_LB_Stage,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2c - UPDRS (off meds) ##
UDPRS_LBPLOT <- ggplot(PATH, aes(x = motor_updrs_off, y = mt_copy_number_avg_Log10)) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  geom_smooth(method = 'lm', se = T, colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

UDPRS_LBPLOT + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                     axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="UDPRS (off)", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

cor.test(PATH$mt_copy_number_avg_Log10, PATH$motor_updrs_off, alternative = c("two.sided"), method = c("spearman"))

## Graph and Stats - Figure 2d - LB Scores - Cranial Nerves ix x ##
PATH$brain_stem_ix_x_character <- as.character(PATH$brain_stem_ix_x)

CRANIAL_LBPLOT <- ggplot(PATH, aes(x = brain_stem_ix_x_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

CRANIAL_LBPLOT + theme(panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                       panel.grid.major = element_blank(), 
                       axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Cranial Nerves ix x", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ brain_stem_ix_x_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ brain_stem_ix_x_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2e - LB Scores - Substantia Nigra ##
PATH$brain_stem_sn_character<- as.character(PATH$brain_stem_sn)

SN_LBPLOT <- ggplot(PATH, aes(x = brain_stem_sn_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SN_LBPLOT + theme(panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                  panel.grid.major = element_blank(), 
                  axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Substantia Nigra", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ brain_stem_sn_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ brain_stem_sn_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2f - LB Scores - Locus Coeruleus ##
PATH$brain_stem_lc_character<- as.character(PATH$brain_stem_lc)

LC_LBPLOT <- ggplot(PATH, aes(x = brain_stem_lc_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

LC_LBPLOT + theme(panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                  panel.grid.major = element_blank(), 
                  axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Locus Coeureleus", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ brain_stem_lc_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ brain_stem_lc_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Amygdala ##
PATH$bf_amygdala_character<- as.character(PATH$bf_amygdala)

AMY_LBPLOT <- ggplot(PATH, aes(x = bf_amygdala_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

AMY_LBPLOT + theme(panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                   panel.grid.major = element_blank(), 
                   axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Amygdala", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ bf_amygdala_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ bf_amygdala_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Transentorhinal Cortex ##
PATH$bf_trans_character<- as.character(PATH$bf_trans)

TRANS_LBPLOT <- ggplot(PATH, aes(x = bf_trans_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

TRANS_LBPLOT + theme(panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                     panel.grid.major = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Transent Cort", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ bf_trans_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ bf_trans_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Cingulate Gyrus ##
PATH$bf_cing_character<- as.character(PATH$bf_cing)

CING_LBPLOT <- ggplot(PATH, aes(x = bf_cing_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

CING_LBPLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Cingulatet", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ bf_cing_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ bf_cing_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Temporal Lobe ##
PATH$nctx_temporal_character<- as.character(PATH$nctx_temporal)

TEMP_LBPLOT <- ggplot(PATH, aes(x = nctx_temporal_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

TEMP_LBPLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Temporal Lobe", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ nctx_temporal_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ nctx_temporal_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Frontal Lobe ##
PATH$nctx_frontal_character<- as.character(PATH$nctx_frontal)

FRONT_LBPLOT <- ggplot(PATH, aes(x = nctx_frontal_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

FRONT_LBPLOT + theme(panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                     panel.grid.major = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Frontal Lobe", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ nctx_frontal_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ nctx_frontal_character,
         data=PATH,
         method="bonferroni")

## Graph and Stats - Figure 2g - LB Scores - Parietal Lobe ##
PATH$nctx_parieta_character<- as.character(PATH$nctx_parietal)

PARI_LBPLOT <- ggplot(PATH, aes(x = nctx_parieta_character, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

PARI_LBPLOT + theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                    panel.grid.major = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Parietal Lobe", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ nctx_parieta_character, data = PATH)

dunnTest(mt_copy_number_avg_Log10 ~ nctx_parieta_character,
         data=PATH,
         method="bonferroni")

##############################################################################################
## Graphs and Stats - Figure 3 ##
## Loading Data for Figure 3 ##
PATH <- read.table('~/Downloads/Extended_Data_Table2_v3.txt', header=TRUE, sep="\t")

## Log10 transformation of mtDNA copy number & selection of PD subjects##
PATH$mt_copy_number_avg_Log10 <- log10(PATH$mt_copy_number_avg)
PATH$SEX <- as.factor(PATH$SEX)
PATH.PD <- subset(PATH, Stage_early_mid_late!="NA")

## robust linear models (rlm) of LB scores and predicitions of CB mtDNA CN ##
rlm_obt = rlm(log(mt_copy_number_avg) ~ obt + SEX + Age, PATH.PD)
rlm_bs_ixx = rlm(log(mt_copy_number_avg) ~ brain_stem_ix_x + SEX + Age, PATH.PD)
rlm_bs_lc = rlm(log(mt_copy_number_avg) ~ brain_stem_lc + SEX + Age, PATH.PD)
rlm_bs_sn = rlm(log(mt_copy_number_avg) ~ brain_stem_sn + SEX + Age, PATH.PD)
rlm_bf_amy = rlm(log(mt_copy_number_avg) ~ bf_amygdala + SEX + Age, PATH.PD)
rlm_bf_trans = rlm(log(mt_copy_number_avg) ~ bf_trans + SEX + Age, PATH.PD)
rlm_bf_cing = rlm(log(mt_copy_number_avg) ~ bf_cing + SEX + Age, PATH.PD)
rlm_nctx_temp = rlm(log(mt_copy_number_avg) ~ nctx_temporal + SEX + Age, PATH.PD)
rlm_nctx_front = rlm(log(mt_copy_number_avg) ~ nctx_frontal + SEX + Age, PATH.PD)
rlm_nctx_pari = rlm(log(mt_copy_number_avg) ~ nctx_parietal + SEX + Age, PATH.PD)

## summary of rlm results ##
summary(rlm_obt)
summary(rlm_bs_ixx)
summary(rlm_bs_lc)
summary(rlm_bs_sn)
summary(rlm_bf_amy)
summary(rlm_bf_cing)
summary(rlm_bf_trans)
summary(rlm_nctx_front)
summary(rlm_nctx_pari)
summary(rlm_nctx_temp)

## Wald tests of LB score on rlms ##
f.robftest(rlm_obt, var = "obt")
f.robftest(rlm_bs_ixx, var = "brain_stem_ix_x")
f.robftest(rlm_bs_lc, var = "brain_stem_lc")
f.robftest(rlm_bs_sn, var = "brain_stem_sn")
f.robftest(rlm_bf_amy, var = "bf_amygdala")
f.robftest(rlm_bf_cing, var = "bf_cing")
f.robftest(rlm_bf_trans, var = "bf_trans")
f.robftest(rlm_nctx_front, var = "nctx_frontal")
f.robftest(rlm_nctx_pari, var = "nctx_parietal")
f.robftest(rlm_nctx_temp, var = "nctx_temporal")

## plotting models in Figure 3 ##

##bs_ixx##
p2 <- plot_model(rlm_bs_ixx, type = "pred", terms = ("brain_stem_ix_x"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="bs_ixx", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##bs_lc##
p3 <- plot_model(rlm_bs_lc, type = "pred", terms = ("brain_stem_lc"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="bs_lc", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##bs_sn##
p4 <- plot_model(rlm_bs_sn, type = "pred", terms = ("brain_stem_sn"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="brain_stem_sn", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500) 

##bf_amy##
p5 <- plot_model(rlm_bf_amy, type = "pred", terms = ("bf_amygdala"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p5 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="bf_amygdala", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##bf_cing##
p6 <- plot_model(rlm_bf_cing, type = "pred", terms = ("bf_cing"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p6 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="bf_cing", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

## bf_trans ##

p7 <- plot_model(rlm_bf_trans, type = "pred", terms = ("bf_trans"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p7 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="bf_trans", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##nctx_front##
p8 <- plot_model(rlm_nctx_front, type = "pred", terms = ("nctx_frontal"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p8 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="nctx_frontal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##nctx_temp##
p9 <- plot_model(rlm_nctx_temp, type = "pred", terms = ("nctx_temporal"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p9 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
           axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="nctx_temporal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##nctx_pari##
p10 <- plot_model(rlm_nctx_pari, type = "pred", terms = ("nctx_parietal"), colors="black", line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p10 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="nctx_parietal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2500)

##############################################################################################
## Graphs - Supplementary Figure 3 ##

p_a <- plot_model(rlm_obt, type = "pred", terms = c("obt", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="obt", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_b <- plot_model(rlm_bs_ixx, type = "pred", terms = c("brain_stem_ix_x", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_b + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="bs_ixx", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_c <- plot_model(rlm_bs_lc, type = "pred", terms = c("brain_stem_lc", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_c + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="brain_stem_lc", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_d <- plot_model(rlm_bs_sn, type = "pred", terms = c("brain_stem_sn", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="brain_stem_sn", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_e <- plot_model(rlm_bf_amy, type = "pred", terms = c("bf_amygdala", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_e + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="bf_amygdala", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_f <- plot_model(rlm_bf_cing, type = "pred", terms = c("bf_cing", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="bf_cing", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_g <- plot_model(rlm_bf_trans, type = "pred", terms = c("bf_trans", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="bf_trans", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_h <- plot_model(rlm_nctx_front, type = "pred", terms = c("nctx_frontal", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_h + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="nctx_frontal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_i <- plot_model(rlm_nctx_temp, type = "pred", terms = c("nctx_temporal", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_i + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="nctx_temporal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)

p_j <- plot_model(rlm_nctx_pari, type = "pred", terms = c("nctx_parietal", "Age", "SEX"), colors = c("#287D8EFF","#B8DE29FF", "#DCE319FF"), line.size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_j + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
            axis.line = element_line(colour = "black"), axis.text = element_text(size= 20)) +
  labs(x ="nctx_parietal", y = "Cerebellum Mitochondrial Copy Number") +
  ylim(500, 2800)


##############################################################################################
## Graphs and Stats - Supplementary Figure 1 ##

## Subsetting only Samples with Clinical Stage (Years of Disease before death) ##
PATH_clin <- subset(PATH, Stage_early_mid_late!="NA")
PATH_clin$Stage_early_mid_late <- as.factor(PATH_clin$Stage_early_mid_late)
levels(PATH_clin$Stage_early_mid_late)
PATH_clin$Stage_early_mid_late <- factor(PATH_clin$Stage_early_mid_late, levels=c('early', 'mid', 'late'))

## Graphs and Stats - Supplementary Figure 1a - Clinical Stage (early, mid, late) ##
CLINSTAGE_PLOT <- ggplot(PATH_clin, aes(x = Stage_early_mid_late, y = mt_copy_number_avg_Log10)) +
  geom_boxplot(width = 0.3, lwd=0.8,outlier.shape = 4, outlier.size = 4) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=6, color="black", fill="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

CLINSTAGE_PLOT + theme(panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                       panel.grid.major = element_blank(), 
                       axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="Stage based on Years since Dx", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

kruskal.test(mt_copy_number_avg_Log10 ~ Stage_early_mid_late, data = PATH_clin)

## Graph and Stats - Figure 1b - UPDRS (off meds) - PD only ##
UDPRS_PDonly <- ggplot(PATH_clin, aes(x = motor_updrs_off, y = mt_copy_number_avg_Log10)) +
  geom_jitter(aes(color= mt_copy_number_avg_Log10), width=0.1, show.legend = FALSE) +
  geom_smooth(method = 'lm', se = T, colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

UDPRS_PDonly + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "azure1", colour = "black",size = 1.5, linetype = "solid"),
                     axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  labs(x ="UDPRS (off)", y = "Cerebellum Mitochondrial Copy Number") +
  scale_color_viridis_c(option="viridis",limits= c(2.0,4.0)) +
  ylim(2.53, 3.88)

cor.test(PATH_clin$mt_copy_number_avg_Log10, PATH_clin$motor_updrs_off, alternative = c("two.sided"), method = c("spearman"))

##############################################################################################
## Stats for Table 2 and Table 3 ##

## Calculating Sum LB Scores ##
PATH_clin$sumbrainstem <- PATH_clin$brain_stem_ix_x + PATH_clin$brain_stem_lc + PATH_clin$brain_stem_sn
PATH_clin$sumlimbic <- PATH_clin$bf_amygdala + PATH_clin$bf_trans + PATH_clin$bf_cing
PATH_clin$sumneoc <- PATH_clin$nctx_temporal + PATH_clin$nctx_frontal + PATH_clin$nctx_parietal

## Binning PD samples for CB mtDNA CN < or > 2000 ##
less2000 <- subset(PATH_clin, mt_copy_number_avg < 2000)
over2000 <- subset(PATH_clin, mt_copy_number_avg > 2000)

## Binning PD samples for CB mtDNA CN < or > 4000 ##
less4000 <- subset(PATH_clin, mt_copy_number_avg < 4000)
over4000 <- subset(PATH_clin, mt_copy_number_avg > 4000)

## Stats Table 2 - LB Sum scores by bin (< or > 2000) ##
wilcox.test(less2000$Sum_LB_score, over2000$Sum_LB_score, alternative = c("two.sided"))
wilcox.test(less2000$sumbrainstem, over2000$sumbrainstem, alternative = c("two.sided"))
wilcox.test(less2000$sumlimbic, over2000$sumlimbic, alternative = c("two.sided"))
wilcox.test(less2000$sumneoc, over2000$sumneoc, alternative = c("two.sided"))

## Stats Table 2 - Individual Brain Regions - LB scores by bin (< or > 2000) ##
wilcox.test(less2000$obt, over2000$obt, alternative = c("two.sided"))
wilcox.test(less2000$brain_stem_ix_x, over2000$brain_stem_ix_x, alternative = c("two.sided"))
wilcox.test(less2000$brain_stem_sn, over2000$brain_stem_sn, alternative = c("two.sided"))
wilcox.test(less2000$brain_stem_lc, over2000$brain_stem_lc, alternative = c("two.sided"))
wilcox.test(less2000$bf_amygdala, over2000$bf_amygdala, alternative = c("two.sided"))
wilcox.test(less2000$bf_trans, over2000$bf_trans, alternative = c("two.sided"))
wilcox.test(less2000$bf_cing, over2000$bf_cing, alternative = c("two.sided"))
wilcox.test(less2000$nctx_temporal, over2000$nctx_temporal, alternative = c("two.sided"))
wilcox.test(less2000$nctx_frontal, over2000$nctx_frontal, alternative = c("two.sided"))
wilcox.test(less2000$nctx_parietal, over2000$nctx_parietal, alternative = c("two.sided"))

## Stats Table 2 - UPDRS (off meds) - LB scores by bin (< or > 2000) ##
wilcox.test(less2000$motor_updrs_off, over2000$motor_updrs_off, alternative = c("two.sided"))
wilcox.test(less2000$motor_updrs_score_months_prior_to_death, over2000$motor_updrs_score_months_prior_to_death, alternative = c("two.sided"))

## Stats Supplementary Table 2 - LB Sum scores by bin (< or > 4000) ##
wilcox.test(less4000$Sum_LB_score, over4000$Sum_LB_score, alternative = c("two.sided"))
wilcox.test(less4000$sumbrainstem, over4000$sumbrainstem, alternative = c("two.sided"))
wilcox.test(less4000$sumlimbic, over4000$sumlimbic, alternative = c("two.sided"))
wilcox.test(less4000$sumneoc, over4000$sumneoc, alternative = c("two.sided"))

## Stats Supplementary Table 2 - Individual Brain Regions - LB scores by bin (< or > 4000) ##
wilcox.test(less4000$obt, over4000$obt, alternative = c("two.sided"))
wilcox.test(less4000$brain_stem_ix_x, over4000$brain_stem_ix_x, alternative = c("two.sided"))
wilcox.test(less4000$brain_stem_sn, over4000$brain_stem_sn, alternative = c("two.sided"))
wilcox.test(less4000$brain_stem_lc, over4000$brain_stem_lc, alternative = c("two.sided"))
wilcox.test(less4000$bf_amygdala, over4000$bf_amygdala, alternative = c("two.sided"))
wilcox.test(less4000$bf_trans, over4000$bf_trans, alternative = c("two.sided"))
wilcox.test(less4000$bf_cing, over4000$bf_cing, alternative = c("two.sided"))
wilcox.test(less4000$nctx_temporal, over4000$nctx_temporal, alternative = c("two.sided"))
wilcox.test(less4000$nctx_frontal, over4000$nctx_frontal, alternative = c("two.sided"))
wilcox.test(less4000$nctx_parietal, over4000$nctx_parietal, alternative = c("two.sided"))

## Stats Supplementary Table 2 - UPDRS (off meds) - LB scores by bin (< or > 4000) ##
wilcox.test(less4000$motor_updrs_off, over4000$motor_updrs_off, alternative = c("two.sided"))
wilcox.test(less4000$motor_updrs_score_months_prior_to_death, over4000$motor_updrs_score_months_prior_to_death, alternative = c("two.sided"))

##############################################################################################
## Plot for Supplementary Figure 2 ##

## seleciton of LB scores from PD subjects ##
LB_scores_only <- select(PATH.PD, obt, brain_stem_ix_x, brain_stem_lc, brain_stem_sn, bf_amygdala, bf_trans, bf_cing, nctx_temporal, nctx_frontal, nctx_parietal)  

## Spearman's correlation plot of LB scores ##
LB_cor <- cor(LB_scores_only, method="spearman")
ggcorrplot(LB_cor,hc.order = TRUE, type = "lower",
           outline.col = "white",  lab=TRUE)

