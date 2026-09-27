# Drafted by Nichole Ginnan (nginn001@ucr.edu)
# August 2025
# Analyses for CellTiter Blue quantitative fluorence growth assay data for Luteibacter strains collected across Kansas, USA

# load libraries ####
library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
set.seed(4444)

# load data ####
NaCl_data <- read_excel("Supp.Table.S3.CellTiterBlue.assay.measurements.that.passed.qualifty.filtering.xlsx", 
                        sheet = "NaCl")
PEG_data <- read_excel("Supp.Table.S3.CellTiterBlue.assay.measurements.that.passed.qualifty.filtering.xlsx", 
                       sheet = "PEG")

##############################################################
# Prepare CellTiter Blue assay data ##########################
##############################################################
# Replace negative values in the "RFU" column with zero, biologically a negative value = no growth
NaCl_data$RFU[NaCl_data$RFU < 0] <- 0
PEG_data$RFU[PEG_data$RFU < 0] <- 0
# Remove blanks/controls from in vivo growth assay data
NaCl_data <- NaCl_data %>% filter(Strain != "BLANK")
PEG_data <- PEG_data %>% filter(Strain != "BLANK")
##### Subset controls (0.0%) ######
control <- NaCl_data %>% filter(NaCl.concentration == 0.0) # all strains
control.seq <- control %>% filter(seq == "sequenced") # only strains that we sequenced whole genomes for

##############################################################
# Fig. 1b. Luteibacter growth under controlled conditions ####
##############################################################
#### NaCl control plot grouped by Site - ALL ####
tehc.rep.count<-control %>% group_by(Strain) %>% summarize(count = n()) # sample size, biological reps; 3-6 technical reps per biological rep
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_tech.rep.count.xlsx")
control %>% group_by(Site.2) %>% summarize(count = n()) # sample size, includes all technical reps separatly
#KNZ      332
#TLI      342
#WK       345

# Average the technical reps for each biological rep
avg.tech.reps.df <- control %>% group_by(Strain, Site.2) %>% summarise(avg_RFU = mean(RFU, na.rm = TRUE),
    across(-RFU, ~ first(.)),   # Keeps the first value of all other columns except RFU
    .groups = "drop")
avg.tech.reps.df %>% group_by(Site.2) %>% summarize(count = n()) # sample size, biological reps per site
#KNZ       60
#TLI       63
#WK        63
#### model the biological reps ##
model <- lm(sqrt(avg_RFU)~Site.2, data = avg.tech.reps.df)
plot(resid(model)~fitted(model))
anova(model)
#           Df     Sum Sq    Mean Sq  F value    Pr(>F)    
#Site.2     2    244152501 122076250  49.117 < 2.2e-16 ***
#Residuals 183 454826570   2485391 
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit (RFU)")+
  theme_classic()
ggsave(filename= "control.RFU.site.BioReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=70, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

#######################################################################
# Fig. 2a. Luteibacter growth under increasing NaCl concentrations ####
#######################################################################
# Calculate the RFU ratio for each strain: growth under stress/ growth under control conditions
baseline_rfu <- NaCl_data %>% # average the technical replicates for each strain's growth under standard media ("growth baseline")
  filter(NaCl.concentration == 0.0) %>%
  group_by(Strain) %>%
  summarize(RFU_baseline = mean(RFU, na.rm = TRUE))
NaCl_data_with_baseline <- NaCl_data %>% # Merge the baseline data with the original data to have the baseline RFU for each strain
  left_join(baseline_rfu, by = "Strain")
NaCl_data_relative.growth <- NaCl_data_with_baseline %>% # For NaCl concentrations of 0.3, 1.8, and 3.0, divide the RFU by the baseline RFU for each strain
  mutate(
    RFU_ratio = ifelse(NaCl.concentration %in% c(0.3, 1.8, 3), RFU / RFU_baseline, NA)
  )
NaCl_data_relative.growth
# subset by NaCl concentration
NaCl.0.0.ratio <- NaCl_data_relative.growth %>% filter(NaCl.concentration == 0.0)
NaCl.0.3.ratio <- NaCl_data_relative.growth %>% filter(NaCl.concentration == 0.3)
NaCl.1.8.ratio <- NaCl_data_relative.growth %>% filter(NaCl.concentration == 1.8)
NaCl.3.0.ratio<- NaCl_data_relative.growth %>% filter(NaCl.concentration == 3)
# Count technical rep per biological rep in each concentration
tehc.rep.count<-NaCl.0.3.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_0.3_tech.rep.count.xlsx")
tehc.rep.count<-NaCl.1.8.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_1.8_tech.rep.count.xlsx")
tehc.rep.count<-NaCl.3.0.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_3.0_tech.rep.count.xlsx")
# Average the technical reps for each biological rep for each concentration
avg.tech.reps.df.NaCl.0.0.ratio <- NaCl.0.0.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_baseline = mean(RFU_baseline, na.rm = TRUE),across(-RFU_baseline, ~ first(.)),.groups = "drop")
avg.tech.reps.df.NaCl.0.3.ratio <- NaCl.0.3.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
avg.tech.reps.df.NaCl.1.8.ratio <- NaCl.1.8.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
write.xlsx(avg.tech.reps.df.NaCl.1.8.ratio, "/Users/nicholeginnan/Documents/KU/Luteibacter/mGWAS_nichole/avg.tech.reps.df.NaCl.1.8.ratio.xlsx")
avg.tech.reps.df.NaCl.3.0.ratio <- NaCl.3.0.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
# Count sample size, biological reps per site
avg.tech.reps.df.NaCl.0.3.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=60; TLI=63; WK=63
avg.tech.reps.df.NaCl.1.8.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=60; TLI=63; WK=63
avg.tech.reps.df.NaCl.3.0.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=60; TLI=63; WK=63
##############################################################################
### NaCl 3.0% model ####
model <- lm(log10(avg_RFU_ratio)~Site.2, data = avg.tech.reps.df.NaCl.3.0.ratio)
plot(resid(model)~fitted(model))
anova(model)
#           Df Sum Sq Mean Sq F value Pr(>F)
#Site.2      2  1.340 0.67025  1.9806 0.1409
#Residuals 183 61.928 0.33840               
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "NaCl.3.RFU.site.biologicalReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

### NaCl 1.8% model ####
outliers <- c("SVR6-11D")  # outlier (all had very small avg_RFU_ration (0.002158))
avg.tech.reps.df.NaCl.1.8.ratio.2 <- avg.tech.reps.df.NaCl.1.8.ratio %>% filter(!Strain %in% outliers) #remove outlier
model <- lm(log10(avg_RFU_ratio)~Site.2, data = avg.tech.reps.df.NaCl.1.8.ratio.2)
plot(resid(model)~fitted(model))
anova(model)
#           Df Sum Sq Mean Sq F value  Pr(>F)  
#Site.2      2  2.998 1.49878  3.3471 0.03737 *
#Residuals 182 81.498 0.44779                 
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.08, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "NaCl.1.8.RFU.site.bioReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

### NaCl 0.3% model ####
outliers <- c("SVR6-11D","SVR5-9B","SVR6-11A","SVR3-9F","SVR5-1B","SVR3-5C","SVR3-7D","SVR3-12F")  # outliers (very large an very small RFU ratios)
avg.tech.reps.df.NaCl.0.3.ratio.2 <- avg.tech.reps.df.NaCl.0.3.ratio %>% filter(!Strain %in% outliers) #remove outlier
model <- lm(log10(avg_RFU_ratio)~Site.2, data = avg.tech.reps.df.NaCl.0.3.ratio.2)
plot(resid(model)~fitted(model))
anova(model)
#           Df  Sum Sq  Mean Sq F value Pr(>F)
#Site.2      2  0.0153 0.007665  0.0726   0.93
#Residuals 175 18.4872 0.105641               
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "NaCl.0.3.RFU.site.bioReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

#######################################################################
# Fig. 2c. Growth trade offs NaCl 1.8% and control conditions #########
#######################################################################
##### Subset NaCL (1.8%) ######
NaCl.1.8.raw <- NaCl_data %>% filter(NaCl.concentration == 1.8) # all strains
tehc.rep.count.NaCl.1.8.raw<-NaCl.1.8.raw %>% group_by(Strain) %>% summarize(count = n()) # sample size, biological reps; 3-6 technical reps per biological rep
write.xlsx(tehc.rep.count.NaCl.1.8.raw, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_tech.rep.count.NaCl.1.8.raw.xlsx")
NaCl.1.8.raw %>% group_by(Site.2) %>% summarize(count = n()) # sample size, includes all technical reps separatly
# Average the technical reps for each biological rep
avg.tech.reps.df.NaCl.1.8.raw <- NaCl.1.8.raw %>% group_by(Strain, Site.2) %>% summarise(avg_RFU = mean(RFU, na.rm = TRUE),across(-RFU, ~ first(.)),.groups = "drop")
avg.tech.reps.df.NaCl.1.8.raw %>% group_by(Site.2) %>% summarize(count = n()) # sample size, biological reps per site
#KNZ       60
#TLI       63
#WK        63
# check dfs that will be included in this analysis
avg.tech.reps.df
avg.tech.reps.df$avg_RFU # average tech rep per biological reps for control conditions growth

avg.tech.reps.df.NaCl.1.8.raw
avg.tech.reps.df.NaCl.1.8.raw$avg_RFU #average tech rep per biological rep growth in NaCl 1.8% 

#Add suffixes before joining dfs
names(avg.tech.reps.df)[-1] <- paste0(names(avg.tech.reps.df)[-1], "_control")
names(avg.tech.reps.df.NaCl.1.8.raw)[-1] <- paste0(names(avg.tech.reps.df.NaCl.1.8.raw)[-1], "_NaCl")
# Then join the df
merged_df <- left_join(avg.tech.reps.df, avg.tech.reps.df.NaCl.1.8.raw, by = "Strain")
merged_df.filtered <- merged_df %>%filter(Strain != "SVR5-12G") #outlier
#### Regression of control and NaCl 1.8% (averaged technical reps)  ####
model<-lm(log10(avg_RFU_NaCl)~sqrt(avg_RFU_control), data=merged_df.filtered)
plot(resid(model)~fitted(model))
summary(model)
summary(model)$coefficients
#                         Estimate   Std. Error    t value     Pr(>|t|)
#(Intercept)            6.579762e+00 1.134915e-01 57.975832 1.034893e-119
#sqrt(avg_RFU_control) -1.908187e-05 2.034427e-05 -0.937948  3.495074e-01
summary(model)$r.squared
# 0.004784358
confint(model)
#                           2.5 %        97.5 %
#(Intercept)            6.355842e+00 6.803682e+00
#sqrt(avg_RFU_control) -5.922136e-05 2.105762e-05
anova(model)
#                       Df  Sum Sq  Mean Sq  F value  Pr(>F)    
#sqrt(avg_RFU_control)   1  0.245 0.24520    0.8797 0.3495
#Residuals             183 51.005 0.27871               

# Create a new column with the transformed variables for plotting
merged_df.filtered <- merged_df.filtered %>%
  mutate(
    sqrt_control = sqrt(avg_RFU_control),
    log_NaCl = log10(avg_RFU_NaCl) )
# Plot
ggplot(merged_df.filtered, aes(x = sqrt_control, y = log_NaCl, color=Site.2_control)) +
  geom_point(size = 3, alpha = 0.9) +  # data points
  geom_smooth(method = "lm", color = "#2e294e", se = TRUE) +  # regression line + CI
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  #scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  labs(x = "Control conditions growth (sqrt)",y = "NaCl-stressed growth (log10)") +
  theme_classic()+theme(legend.position = "none")
ggsave(filename= "Control.NaCl1.8.tradeoff.regression.nolegend.raw_NaCl_values.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

#######################################################################
# Fig. 2b. Luteibacter growth under increasing PEG concentrations #####
#######################################################################
# Calculate the RFU ratio for each strain: growth under stress/ growth under control conditions
baseline_rfu <- PEG_data %>% # average the technical replicates for each strain's growth under standard media ("growth baseline")
  filter(PEG.concentration == 0) %>%
  group_by(Strain) %>%
  summarize(RFU_baseline = mean(RFU, na.rm = TRUE))
PEG_data_with_baseline <- PEG_data %>% # Merge the baseline data with the original data to have the baseline RFU for each strain
  left_join(baseline_rfu, by = "Strain")
PEG_data_relative.growth <- PEG_data_with_baseline %>% # For PEG concentrations of 10, 30, and 40, divide the RFU by the baseline RFU for each strain
  mutate(
    RFU_ratio = ifelse(PEG.concentration %in% c(10, 30, 40), RFU / RFU_baseline, NA)
  )
PEG_data_relative.growth
# subset by PEG concentration
PEG.0.baseline <- PEG_data_relative.growth %>% filter(PEG.concentration == 0)
PEG.10.ratio <- PEG_data_relative.growth %>% filter(PEG.concentration == 10)
PEG.30.ratio <- PEG_data_relative.growth %>% filter(PEG.concentration == 30)
PEG.40.ratio<- PEG_data_relative.growth %>% filter(PEG.concentration == 40)
# Count technical rep per biological rep in each concentration
tehc.rep.count<-PEG.10.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/PEG_10_tech.rep.count.xlsx")
tehc.rep.count<-PEG.30.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/PEG_30_tech.rep.count.xlsx")
tehc.rep.count<-PEG.40.ratio %>% group_by(Strain) %>% summarize(count = n()) 
write.xlsx(tehc.rep.count, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/PEG_40_tech.rep.count.xlsx")
# Average the technical reps for each biological rep for each concentration
avg.tech.reps.df.PEG.0.baseline <- PEG.0.baseline %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_baseline = mean(RFU_baseline, na.rm = TRUE),across(-RFU_baseline, ~ first(.)),.groups = "drop")
avg.tech.reps.df.PEG.10.ratio <- PEG.10.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
avg.tech.reps.df.PEG.30.ratio <- PEG.30.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
write.xlsx(avg.tech.reps.df.PEG.30.ratio, "/Users/nicholeginnan/Documents/KU/Luteibacter/mGWAS_nichole/avg.tech.reps.df.PEG.30.ratio.xlsx")
avg.tech.reps.df.PEG.40.ratio <- PEG.40.ratio %>% group_by(Strain, Site.2) %>% summarise(avg_RFU_ratio = mean(RFU_ratio, na.rm = TRUE),across(-RFU_ratio, ~ first(.)),.groups = "drop")
# Count sample size, biological reps per site
avg.tech.reps.df.PEG.10.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=50; TLI=51; WK=54
avg.tech.reps.df.PEG.30.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=50; TLI=51; WK=54
avg.tech.reps.df.PEG.40.ratio %>% group_by(Site.2) %>% summarize(count = n()) # KNZ=50; TLI=51; WK=54
##############################################################################
### PEG 40% model ####
model <- lm(log10(avg_RFU_ratio+0.00001)~Site.2, data = avg.tech.reps.df.PEG.40.ratio)
plot(resid(model)~fitted(model))
anova(model)
#           Df  Sum Sq  Mean Sq  F value  Pr(>F)
#Site.2      2   5.645  2.8227  2.1423 0.1209
#Residuals 152  200.279  1.3176              
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "PEG.40.RFU.site.biologicalReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

### PEG 30% model ####
outliers <- c("KNZ2-4C","SVR3-9D","KNZ1-7G","SVR5-5F")  # outliers
avg.tech.reps.df.PEG.30.ratio.2 <- avg.tech.reps.df.PEG.30.ratio %>% filter(!Strain %in% outliers) #remove outlier

model <- lm(log10(avg_RFU_ratio)~Site.2, data = avg.tech.reps.df.PEG.30.ratio.2)
plot(resid(model)~fitted(model))
anova(model)
#           Df Sum Sq Mean Sq F value  Pr(>F)  
#Site.2     2  2.954  1.4770  3.5642 0.03078 *
#Residuals 148 61.332  0.4144                
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.08, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "PEG.30.RFU.site.bioReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")

### PEG 10% model ####
model <- lm(log10(avg_RFU_ratio)~Site.2, data = avg.tech.reps.df.PEG.10.ratio)
plot(resid(model)~fitted(model))
anova(model)
#           Df  Sum Sq  Mean Sq F value Pr(>F)
#Site.2      2  0.803 0.40169  1.2623  0.286
#Residuals 152 48.369 0.31822               
emm_g1<-emmeans(model, specs=~Site.2, type="response")
cld<-cld(emm_g1, method='tukey', alpha=0.05, adjust='fdr', Letters=LETTERS, sort=TRUE, reversed=TRUE) # fdr or sidak
cld %>% as.data.frame %>%
  ggplot(., aes(x=Site.2, y=response, color=Site.2, fill=Site.2))+
  geom_pointrange(aes(ymin=response-SE, ymax=response+SE),color="black", size=1, position=position_dodge2(0.65)) +
  geom_point(aes(group=Site.2),size=3.9, position=position_dodge2(width=0.65))+
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  geom_text(aes(label=.group, y = response+SE, group=Site.2),color="black",size=4, vjust = -0.5,hjust=0.5, position=position_dodge(0.7))+
  scale_x_discrete(limit=c('WK', 'TLI', 'KNZ'))+
  xlab("Origin Site") +
  ylab("Relative fluorescence unit in the context of control")+
  theme_classic()
ggsave(filename= "PEG.10.RFU.site.bioReps.pdf", device="pdf", units="mm", dpi=300, width=70, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")
#######################################################################
# Fig. 2d. Growth trade offs PEG 30% and control conditions #########
#######################################################################
##### Subset PEG (30%) ######
PEG.30.raw <- PEG_data %>% filter(PEG.concentration == 30) # all strains
tehc.rep.count.PEG.30.raw<-PEG.30.raw %>% group_by(Strain) %>% summarize(count = n()) # sample size, biological reps; 3-6 technical reps per biological rep
write.xlsx(tehc.rep.count.PEG.30.raw, "/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/NaCl_tech.rep.count.PEG.30.raw.xlsx")
PEG.30.raw %>% group_by(Site.2) %>% summarize(count = n()) # sample size, includes all technical reps separatly
# Average the technical reps for each biological rep
avg.tech.reps.df.PEG.30.raw <- PEG.30.raw %>% group_by(Strain, Site.2) %>% summarise(avg_RFU = mean(RFU, na.rm = TRUE),across(-RFU, ~ first(.)),.groups = "drop")
avg.tech.reps.df.PEG.30.raw %>% group_by(Site.2) %>% summarize(count = n()) # sample size, biological reps per site
#KNZ       50
#TLI       51
#WK        54
# check dfs that will be included in this analysis
avg.tech.reps.df
avg.tech.reps.df$avg_RFU # average tech rep per biological reps for control conditions growth

avg.tech.reps.df.PEG.30.raw
avg.tech.reps.df.PEG.30.raw$avg_RFU_PEG #average tech rep per biological rep growth in NaCl 1.8% 

#Add suffixes before joining dfs
names(avg.tech.reps.df)[-1] <- paste0(names(avg.tech.reps.df)[-1], "_control")
names(avg.tech.reps.df.PEG.30.raw)[-1] <- paste0(names(avg.tech.reps.df.PEG.30.raw)[-1], "_PEG")
# Then join the df
merged_df <- left_join(avg.tech.reps.df, avg.tech.reps.df.PEG.30.raw, by = "Strain")
merged_df.filtered <- merged_df %>%filter(!Strain %in% c("KNZ2-4C", "SVR5-12G", "SVR3-9D"))  # outliers
#### Regression of control and PEG 30% (averaged technical reps)  ####
model<-lm(log10(avg_RFU_PEG)~sqrt(avg_RFU_control), data=merged_df.filtered)
plot(resid(model)~fitted(model))
summary(model)
summary(model)$coefficients
#                         Estimate   Std. Error    t value     Pr(>|t|)
#(Intercept)            6.951233e+00 1.079978e-01 64.364596 3.733158e-111
#sqrt(avg_RFU_control) -4.145774e-05 2.000818e-05 -2.072039  3.997376e-02
summary(model)$r.squared
#  0.02782587
confint(model)
#                           2.5 %        97.5 %
#(Intercept)            6.737840e+00  7.164626e+00
#sqrt(avg_RFU_control) -8.099202e-05 -1.923463e-06
anova(model)
#                       Df  Sum Sq  Mean Sq  F value  Pr(>F)    
#sqrt(avg_RFU_control)   1  0.942 0.94190  4.2933 0.03997 *
#Residuals             150 32.908 0.21938              

# Create a new column with the transformed variables for plotting
merged_df.filtered <- merged_df.filtered %>%
  mutate(
    sqrt_control = sqrt(avg_RFU_control),
    log_PEG = log10(avg_RFU_PEG) )
# Plot
ggplot(merged_df.filtered, aes(x = sqrt_control, y = log_PEG, color=Site.2_control)) +
  geom_point(size = 3, alpha = 0.9) +  # data points
  geom_smooth(method = "lm", color = "#2e294e", se = TRUE) +  # regression line + CI
  scale_color_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) + 
  #scale_fill_manual(values=c("#aa4465", "#ffa69e","#3d5a80"),limits=c("WK", "TLI", "KNZ")) +
  labs(x = "Control conditions growth (sqrt)",y = "PEG-stressed growth (log10)") +
  theme_classic()+theme(legend.position = "none")
ggsave(filename= "Control.PEG1.8.tradeoff.regression.nolegend.with.raw.PEG.growth.pdf", device="pdf", units="mm", dpi=300, width=100, height=100, path="/Users/nicholeginnan/Documents/KU/Luteibacter/plots")
