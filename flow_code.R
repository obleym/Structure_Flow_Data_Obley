# FLOW ----
# R code for looking at variance difference between length among families in control vs treatment tanks
install.packages("lme4", "lmerTest", "tidyverse", "ggplot", "ggrepel", "dplyr", "ggcorrplot")
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(ggcorrplot)
library(glmmTMB)

### reading the data into R ----
OHRC.dat.fam <- OHRC_data_forR <- read_csv("Google_Drive/Graduate/MS/FISH SCIENCE/Data Analysis/Flow/OHRC_data_forR.csv")
View(OHRC.dat.fam)

OHRC.dat.fam <- na.omit(OHRC.dat.fam)
f.zzq=as.data.frame(table(OHRC.dat.fam[c("treatment", "tank","family")]))
f.zzq=subset(f.zzq, Freq>0)
# shows how many fish were in each family in each tank
f.zzq




# Before adjusting for tank effect ----
## Subset flow and control tanks ----
flow= subset(OHRC.dat.fam, treatment == 'flow')
f.control= subset(OHRC.dat.fam, treatment == 'f.control')

## family means (NOT adjusted) ----
f.fam.mean=aggregate(OHRC.dat.fam[c("length")],
                     by= OHRC.dat.fam[c("family","treatment")], mean)

tA = subset(OHRC.dat.fam, tank == 'A')
tB = subset(OHRC.dat.fam, tank == 'B')
tC = subset(OHRC.dat.fam, tank == 'C')
tD = subset(OHRC.dat.fam, tank == 'D')
tE = subset(OHRC.dat.fam, tank == 'E')
tF = subset(OHRC.dat.fam, tank == 'F')
tH = subset(OHRC.dat.fam, tank == 'H')
tI = subset(OHRC.dat.fam, tank == 'I')

## tank ICCs ----
require(ICC)

ICCest(family, length, data = tA, CI.type = "S")
ICCest(family, length, data = tB, CI.type = "S")
ICCest(family, length, data = tC, CI.type = "S")
ICCest(family, length, data = tD, CI.type = "S")
ICCest(family, length, data = tE, CI.type = "S")
ICCest(family, length, data = tF, CI.type = "S")
ICCest(family, length, data = tH, CI.type = "S")
ICCest(family, length, data = tI, CI.type = "S")


### Create tank ICC dataframe ----
f.tank <- c('A','B','C','D','E','F','H','I')
f.treatment <- c('Flow', 'Control', 'Control', 'Flow', 'Control', 'Flow', 'Flow', 'Control')
f.ICC <- c(0.1852441,0.154041,0.1917256,0.1717507,0.1538535,0.2850669,0.1813925,0.167842)
f.varw <- c(1.741834, 2.023114, 1.912486, 1.879169, 1.994413, 1.728293, 1.472022, 1.933854)
f.vara <-c(0.3960259, 0.3683897, 0.4536486, 0.3896757, 0.3626411, 0.6891262, 0.3261804, 0.3900485)


f.ICCdf <- data.frame(f.tank, f.treatment, f.ICC, f.varw, f.vara)
f.ICCdf


### Plot tank ICCs ----
require(ggplot2)
require(ggrepel)

f.ICCplot.tank <- ggplot(data = f.ICCdf, 
                         aes(x = f.treatment, 
                             y = f.ICC)) +
  geom_point(aes(color=f.treatment),
             size=3.5, 
             position = position_dodge2(w = 0.15)) +
  scale_x_discrete(expand=c(0.1, 0.5)) + 
  theme(panel.grid = element_line()) + 
  ggtitle("ICC comparison between\n control and flow tanks") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Treatment", y = "ICC") + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))
f.ICCplot.tank

### geom_label_repel
f.ICCplot.tank + 
  geom_label_repel(aes(label = f.tank),
                   box.padding   = .6,
                   point.padding = .1,
                   min.segment.length = Inf,
                   segment.color = 'grey50') 


## perform Welch's t-test ----
t.test(f.ICC ~ f.treatment, f.ICCdf)

## Family plot ----
f.famplot.na <- ggplot(data = f.fam.mean, aes(x = treatment, y = length, label = family, group = family)) +
  geom_point() + geom_line(aes(color=family)) + scale_y_continuous(limits=c(13,18)) + theme(panel.grid = element_line()) + labs(title="Mean family length in control and flow (not adjusted)") + theme(plot.title = element_text(hjust = 0.5))
f.famplot.na

f.famplot.na + 
  geom_label_repel(data = subset(f.fam.mean, treatment %in% c("f.control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(f.fam.mean, treatment %in% c("flow")),
                   nudge_x       = .2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50')




# Adjust for tank effects ----
# create a dataframe with the average length for each tank
# find the average of each treatment and subtract to find variance of tank from treatment. 
# add this (var) to each fish length to account for tank effect

## find mean tank length ----
f.tankmean <- aggregate(OHRC.dat.fam[c("length")],
                        by= OHRC.dat.fam[c("tank","treatment")], mean)

## find mean treatment length  ----
f.treatmentmean <- aggregate(OHRC.dat.fam[c("length")],
                             by= OHRC.dat.fam[c("treatment")], mean)

### find variance  ----
f.tankmean$var <- c((15.88241-15.68562), (15.51119-15.38397), (15.38443-15.38397), (15.53215-15.68562), (15.09645-15.38397), (15.50170-15.68562), (15.85162-15.68562), (15.55328-15.38397))

### adjust length ----
tA$adjusted <-c(t(tA[, 4] - -0.28752))
tB$adjusted <-c(t(tB[, 4] - 0.19679))
tC$adjusted <-c(t(tC[, 4] - 0.12722))
tD$adjusted <-c(t(tD[, 4] - -0.18392))
tE$adjusted <-c(t(tE[, 4] - 0.00046))
tF$adjusted <-c(t(tF[, 4] - 0.16600))
tH$adjusted <-c(t(tH[, 4] - 0.16931))
tI$adjusted <-c(t(tI[, 4] - -0.16931))

### add to OHRC.dat.fam  ----
f.adjusted = c(tA$adjusted, tB$adjusted, tC$adjusted, tD$adjusted, tE$adjusted, tF$adjusted, tH$adjusted, tI$adjusted)
OHRC.dat.fam['adjusted'] <- f.adjusted
OHRC.dat.fam

## re-subset flow and control tanks to include adjusted values ----
flow= subset(OHRC.dat.fam, treatment == 'flow')
f.control= subset(OHRC.dat.fam, treatment == 'f.control')

## adjusted family means  ----
f.fam.mean.adjusted=aggregate(OHRC.dat.fam[c("adjusted")],
                              by= OHRC.dat.fam[c("family","treatment")], mean)




# After adjusting for tank effect  ----

## treatment ICCS (adjusted) ----
require(ICC)
ICCest(family, f.adjusted, data = flow, CI.type = "S")
# ICC = 0.208972
# varw = 1.780449
# vara = 0.4703549
# upper CI = 0.339061
# lower CI = 0.07888295
ICCest(family, f.adjusted, data = f.control, CI.type = "S")
# ICC = 0.160748
# varw = 2.020777
# vara = 0.3870542
# upper CI = 0.2685716
# lower CI = 0.05292453

# Create ICC dataframe
f.treatment2 <- c('flow', 'control')
f.ICC2 <- c(0.208972,0.160748)
f.vara2 <- c(0.4703549, 0.3870542)
f.varw2 <- c(1.780449, 2.020777)
f.upperCI2 <- c(0.339061, 0.2685716)
f.lowerCI2 <- c(0.07888295, 0.05292453)

f.ICCdf.treat <- data.frame(f.treatment2, f.ICC2, f.vara2, f.varw2, f.upperCI2, f.lowerCI2)
f.ICCdf.treat







### Plot treatment ICCs (adjusted)----
require(ggplot2)
# with error bars
f.ICCplot.treat <- ggplot(data = f.ICCdf.treat, aes(x = f.treatment2, y = f.ICC2)) +
  geom_point(aes(color=f.treatment2),size=3.5, position = position_dodge2(w = 0.15)) +
  scale_x_discrete(expand=c(0.1, 0.5)) + 
  scale_y_continuous(name="ICC", limits=c(0, .4)) +
  theme(panel.grid = element_line()) + 
  ggtitle("ICC comparison between\n Flow and Control (adjusted)") + 
  geom_errorbar(
    aes(x=f.treatment2, 
        ymin = f.lowerCI2, 
        ymax = f.upperCI2), width=0.1, 
    color = "red"
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Treatment", y = "ICC") + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))
f.ICCplot.treat


## CV and opportunity for selection ----
# fish
# control cv
fc.cv.fish <- sd(f.control$adjusted) / mean(f.control$adjusted) * 100
fc.cv.fish # 10.05787

# control ofs
fc.ofs.fish <- ((fc.cv.fish / 100) ^2)
fc.ofs.fish # 0.01011608

# flow cv
ff.cv.fish <- sd(flow$adjusted) / mean(flow$adjusted) * 100
ff.cv.fish # 9.475402

# flow ofs
ff.ofs.fish <- ((ff.cv.fish / 100) ^2)
ff.ofs.fish # 0.008978323

# fam mean
# subset fam means for f.control and flow
f.control.fam.means<-subset(f.fam.mean.adjusted, treatment!="flow")
flow.fam.means<-subset(f.fam.mean.adjusted, treatment!="f.control")

# control cv
fc.cv.fam <- sd(f.control.fam.means$adjusted) / mean(f.control.fam.means$adjusted) * 100
fc.cv.fam # 4.355066

# control ofs
fc.ofs.fam <- ((fc.cv.fam / 100) ^2)
fc.ofs.fam # 0.00189666

# flow cv
ff.cv.fam <- sd(flow.fam.means$adjusted) / mean(flow.fam.means$adjusted) * 100
ff.cv.fam # 4.650409

# flow ofs
ff.ofs.fam <- ((ff.cv.fam / 100) ^2)
ff.ofs.fam # 0.00216263


## Family plots (adjusted) ----
# plot family mean lengths comparing control and flow (adjusted)
f.famplot.a <- ggplot(data = f.fam.mean.adjusted, 
                      aes(x = treatment, 
                          y = adjusted, 
                          label = family, 
                          group = family)) +
  geom_point() + 
  geom_line(aes(color=family)) + 
  scale_y_continuous(limits=c(13,18)) + 
  theme(panel.grid = element_line()) + 
  labs(title="Mean family length in\n control and flow (adjusted)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y= "Length", x = "Treatment") + 
  labs(x="Treatment", y = "Length (cm)") + 
  labs(colour = "Famliy") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))
f.famplot.a

f.famplot.a + 
  geom_label_repel(data = subset(f.fam.mean.adjusted, treatment %in% c("f.control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(f.fam.mean.adjusted, treatment %in% c("flow")),
                   nudge_x       = .2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50')





# F test treat ----
f.ftest.treat <- var.test(adjusted ~ treatment, data = OHRC.dat.fam)
f.ftest.treat
# F test to compare two variances
# 
# data:  adjusted by treatment
# F = 1.0733, num df = 1540, denom df = 1557, p-value = 0.1639
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.9715313 1.1858070
# sample estimates:
#   ratio of variances 
# 1.07332 

### Fligner ----
fligner.test(adjusted ~ treatment, data = OHRC.dat.fam)
# Fligner-Killeen:med chi-squared = 0.084501, df = 1, p-value = 0.7713

#### variances ----
f.control.var.treat <- var(subset(OHRC.dat.fam, treatment == 'f.control', adjusted))
f.control.var.treat
# 2.380918

flow.var.treat <- var(subset(OHRC.dat.fam, treatment == 'flow', adjusted))
flow.var.treat
# 2.218274


# F test fam ----
f.ftest.fam <- var.test(adjusted ~ treatment, data = f.fam.mean.adjusted)
f.ftest.fam
# F test to compare two variances
# 
# data:  adjusted by treatment
# F = 0.83275, num df = 14, denom df = 14, p-value = 0.7368
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.2795782 2.4804135
# sample estimates:
#   ratio of variances 
# 0.8327482  

### Fligner ----
fligner.test(adjusted ~ treatment, data = f.fam.mean.adjusted)
# Fligner-Killeen:med chi-squared = 0.33684, df = 1, p-value = 0.5617


#### variances ----
f.control.var.fam <- var(subset(f.fam.mean.adjusted, treatment == 'f.control', adjusted))
f.control.var.fam
# 0.4432099

flow.var.fam <- var(subset(f.fam.mean.adjusted, treatment == 'flow', adjusted))
flow.var.fam
# 0.5322256








# MIXED MODEL ----
modf1 <- lmer(length ~ 1 + treatment + (1+treatment|family) + (1|treatment/tank), OHRC.dat.fam, REML = F)
summary(modf1)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: length ~ 1 + treatment + (1 + treatment | family) + (1 | treatment/tank)
#    Data: OHRC.dat.fam
# 
#      AIC      BIC   logLik deviance df.resid 
#  10778.5  10826.8  -5381.3  10762.5     3091 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -5.3113 -0.4237  0.1172  0.6392  2.4198 
# 
# Random effects:
#  Groups         Name          Variance Std.Dev. Corr
#  family         (Intercept)   0.40047  0.6328 <----- ## 4.004 is the amount of variance from mean length of fam in the control to mean length of fam in flow. BUT its just telling you the amount of total variation, not the among family. Look at the ICC values to show if the among fam variation is different between tanks and treatments. 

#                 treatmentflow 0.04317  0.2078   0.09 <----- ## 0.431 is the amount of variance of the slopes from control to treatment. So this is how much overall the family slopes vary. But that could be families going up or down. 

#  tank:treatment (Intercept)   0.03173  0.1781       
#  treatment      (Intercept)   0.00000  0.0000       
#  Residual                     1.83723  1.3554       
# Number of obs: 3099, groups:  family, 15; tank:treatment, 8; treatment, 2
# 
# Fixed effects:
#               Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)    15.3322     0.1894 19.8782  80.947   <2e-16 ***
# treatmentflow   0.3321     0.1456  9.3125   2.281   0.0475 *  <----- A clear weak treatment effect
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# treatmntflw -0.304
# optimizer (nloptwrap) convergence code: 0 (OK)

# boundary (singular) fit: see help('isSingular')
anova(mod1)
rand(modf1)

VCrandom <- VarCorr(modf1)
print(VCrandom, comp = "Variance")

modf2 <- lmer(length ~ 1 + treatment + (1+treatment|family) + (1|tank), OHRC.dat.fam, REML = F)
summary(modf2)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: length ~ 1 + treatment + (1 + treatment | family) + (1 | tank)
#    Data: OHRC.dat.fam
# 
#      AIC      BIC   logLik deviance df.resid 
#  10776.5  10818.8  -5381.3  10762.5     3092 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -5.3113 -0.4237  0.1172  0.6392  2.4198 
# 
# Random effects:
#  Groups   Name          Variance Std.Dev. Corr
#  family   (Intercept)   0.40048  0.6328       
#           treatmentflow 0.04318  0.2078   0.09
#  tank     (Intercept)   0.03173  0.1781       
#  Residual               1.83723  1.3554       
# Number of obs: 3099, groups:  family, 15; tank, 8
# 
# Fixed effects:
#               Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)    15.3322     0.1894 19.8778  80.947   <2e-16 ***
# treatmentflow   0.3321     0.1456  9.3135   2.281   0.0475 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# treatmntflw -0.304

# Compare AIC values to determine if "tank" should be nested:
# modf1 = 10778.5
# modf2 = 10776.5  <--- the better model without nested



modf3 <- lmer(length ~ 1 + treatment + family + (+ (1|tank), OHRC.dat.fam, REML = F)))
summary(modf3)
anova(modf3)

# Boxplots ----
# Create data with reordered group levels
OHRC.dat.fam_ordered <- OHRC.dat.fam  
# Order boxes by median
OHRC.dat.fam_ordered.s <- with(OHRC.dat.fam,reorder(family, length, median))
OHRC.dat.fam_ordered$family <- factor(OHRC.dat.fam_ordered$family,levels = levels(OHRC.dat.fam_ordered.s))

# not adjusted
ggplot(OHRC.dat.fam, aes(x=length, color=treatment)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + labs(x="Length(cm)", y = "Count") +
  ggtitle("Lengths for Flow and Control Distribution") + 
  expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5))

# adjusted
ggplot(OHRC.dat.fam, aes(x=adjusted, color=treatment)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + labs(x="Length (cm)", y = "Count") +
  ggtitle("Lengths for Flow and Control Distribution\n (adjusted)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Treatment") +
  expand_limits(x = 0, y = 0) +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))


# not adjusted
ggplot(OHRC.dat.fam_ordered, aes(x = family, y = length, fill = treatment)) + 
  geom_boxplot() +
  labs(title = "Variability of Flow vs. Control for Each Family", 
       x = "family", 
       y = "length") +
  theme_minimal()

# adjusted
ggplot(OHRC.dat.fam_ordered, aes(x = family, y = adjusted, fill = treatment)) + 
  geom_boxplot() +
  labs(title = "Variability of Flow vs. Control \n for Each Family (adjusted)", 
       x = "Family", 
       y = "Length (cm)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))






## control (NOT adjusted) ----

# Create data with reordered group levels
f.control_ordered <- f.control    
# Order boxes by median
fam_ordered.cf <- with(f.control,reorder(family, length, median))
f.control_ordered$family <- factor(f.control_ordered$family,
                                   levels = levels(fam_ordered.cf))
# Draw ggplot2 boxplot ordered by median
box.f.control.na <- ggplot(f.control_ordered,                              
                   aes(x = family,
                       y = length)) +
                geom_boxplot() + 
                labs(y= "Length", x = "Family") + 
                ggtitle("Lengths for Control (not adjusted)") + 
                theme(plot.title = element_text(hjust = 0.5))
box.f.control.na

## flow (NOT adjusted) ----
# Create data with reordered group levels
flow_ordered <- flow    
# Order boxes by median
fam_ordered.ff <- with(flow,reorder(family, length, median))
flow_ordered$family <- factor(flow_ordered$family,
                                   levels = levels(fam_ordered.ff))
# Draw ggplot2 boxplot ordered by median
box.flow.na <- ggplot(flow_ordered,                              
                     aes(x = family,
                         y = length)) +
                geom_boxplot() + 
                labs(y= "Length", x = "Family") + 
                ggtitle("Lengths for Flow (not adjusted)") + 
                theme(plot.title = element_text(hjust = 0.5))
box.flow.na

## control (adjusted) ----

# Create data with reordered group levels
f.control_ordered.a <- f.control    
# Order boxes by median
fam_ordered.cfa <- with(f.control,reorder(family, adjusted, median))
f.control_ordered.a$family <- factor(f.control_ordered.a$family,
                                     levels = levels(fam_ordered.cfa))

# Draw ggplot2 boxplot ordered by median
box.f.control.a <- ggplot(f.control_ordered.a,                              
                       aes(x = family,
                           y = adjusted)) +
                  geom_boxplot() + 
                  labs(y= "Length (adjusted)", x = "Family") + 
                  ggtitle("Lengths (adjusted) for Control") + 
                  theme(plot.title = element_text(hjust = 0.5))
box.f.control.a

## flow (adjusted)----
# Create data with reordered group levels
flow_ordered.a <- flow    
# Order boxes by median
fam_ordered.ffa <- with(flow,reorder(family, adjusted, median))
flow_ordered.a$family <- factor(flow_ordered.a$family,
                                     levels = levels(fam_ordered.ffa))

# Draw ggplot2 boxplot ordered by median
box.flow.a <- ggplot(flow_ordered.a,                              
                   aes(x = family,
                       y = adjusted)) +
              geom_boxplot() + 
              labs(y= "Length (adjusted)", x = "Family") + 
              ggtitle("Lengths (adjusted) for Flow") + 
              theme(plot.title = element_text(hjust = 0.5))    
box.flow.a





# Histograms ----
## control count (NOT adjusted) ----
hist.f.control.na <- ggplot(f.control_ordered, aes(x=length)) + 
                      geom_histogram(colour="black", fill="white") +
                      labs(x="Length(cm)", y = "Count") +
                      ggtitle("Lengths for Control Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5)) 
hist.f.control.na

## control density (NOT adjusted) ----
hist.f.control.den.na <- ggplot(f.control_ordered, aes(x=length)) + 
                  geom_histogram(aes(y=after_stat(density)), colour="black", fill="white") +
                  geom_density(alpha=.2, fill="#FF6666") +
                  labs(x="Length(cm)", y = "Density") +
                  ggtitle("Lengths for Control Distribution") + 
                  theme(plot.title = element_text(hjust = 0.5))
hist.f.control.den.na

## flow count (NOT adjusted) ----
hist.flow.na <- ggplot(flow_ordered, aes(x=length)) + 
                  geom_histogram(colour="black", fill="white") +
                  labs(x="Length(cm)", y = "Count") +
                  ggtitle("Lengths for Flow Distribution") + 
                  theme(plot.title = element_text(hjust = 0.5)) 
hist.flow.na

## flow density (NOT adjusted) ----
hist.flow.den.na <- ggplot(flow_ordered, aes(x=length)) + 
                      geom_histogram(aes(y=after_stat(density)), colour="black", fill="white") +
                      geom_density(alpha=.2, fill="#FF6666") +
                      labs(x="Length(cm)", y = "Density") +
                      ggtitle("Lengths for Flow Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5))
hist.flow.den.na



## control count (adjusted) ----
hist.f.control.a <- ggplot(f.control_ordered, aes(x=adjusted)) + 
                      geom_histogram(colour="black", fill="white") +
                      labs(x="Length(cm)", y = "Count") +
                      ggtitle("Lengths for Control Distribution (adjusted)") + 
                      theme(plot.title = element_text(hjust = 0.5)) 
hist.f.control.a

## control density (adjusted) ----
hist.f.control.den.a <- ggplot(f.control_ordered, aes(x=adjusted)) + 
                          geom_histogram(aes(y=after_stat(density)), colour="black", fill="white") +
                          geom_density(alpha=.2, fill="#FF6666") +
                          labs(x="Length(cm)", y = "Density") +
                          ggtitle("Lengths for Control Distribution (adjusted)") + 
                          theme(plot.title = element_text(hjust = 0.5))
hist.f.control.den.a

## flow count (adjusted) ----
hist.flow.a <- ggplot(flow_ordered, aes(x=adjusted)) + 
                geom_histogram(colour="black", fill="white") +
                labs(x="Length(cm)", y = "Count") +
                ggtitle("Lengths for Flow Distribution (adjusted)") + 
                theme(plot.title = element_text(hjust = 0.5)) 
hist.flow.a

## flow density (adjusted) ----
hist.flow.den.a <- ggplot(flow_ordered, aes(x=adjusted)) + 
                    geom_histogram(aes(y=after_stat(density)), colour="black", fill="white")+
                    geom_density(alpha=.2, fill="#FF6666") +
                    labs(x="Length(cm)", y = "Density") +
                    ggtitle("Lengths for Flow Distribution (adjusted)") + 
                    theme(plot.title = element_text(hjust = 0.5))
hist.flow.den.a


# Q-Q plots ----
## (NOT adjusted) ----
qq.f.na <- ggplot(OHRC.dat.fam, aes(sample = length)) +
            stat_qq(aes(color = treatment)) +
            scale_color_manual(values = c("darkgrey", "black"))+
            labs(y = "Length") +
            ggtitle("Q-Q plot for Control and Flow Length") + 
            theme(plot.title = element_text(hjust = 0.5))
qq.f.na

## (adjusted) ----
qq.f.a <- ggplot(OHRC.dat.fam, aes(sample = adjusted)) +
            stat_qq(aes(color = treatment)) +
            scale_color_manual(values = c("darkgrey", "black"))+
            labs(y = "Length") +
            ggtitle("Q-Q plot for Control and Flow Length (adjusted)") + 
            theme(plot.title = element_text(hjust = 0.5))
qq.f.a

# Pearson correlation matrix ----
library(dplyr)
library(tidyverse)

pcm.data <- read_csv("Google_Drive/Graduate/MS/FISH SCIENCE/Data Analysis/Structure/pcm.data.csv")
View(pcm.data)  
pcm.data <- na.omit(pcm.data)

pcm.data2 <- pcm.data %>% remove_rownames %>% column_to_rownames(var="family")
pcm <- cor(pcm.data2, method = c("spearman"))
pcm
#           f.control      flow s.control structure
# f.control 1.0000000 0.9178571 0.8642857 0.6642857
# flow      0.9178571 1.0000000 0.7607143 0.6964286
# s.control 0.8642857 0.7607143 1.0000000 0.7928571
# structure 0.6642857 0.6964286 0.7928571 1.0000000
install.packages("GGally")
library(GGally) 
data(mtcars) 

# Scatter plot matrix

ggplot(pcm.data2, aes(x=f.control, y=flow)) + 
  geom_point() +
  stat_ellipse() +
  labs(title ="f.control x flow", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))

ggplot(pcm.data2, aes(x=f.control, y=s.control)) + 
  geom_point() +
  stat_ellipse() +
  labs(title ="f.control x s.control", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))

ggplot(pcm.data2, aes(x=flow, y=s.control)) + 
  geom_point() +
  stat_ellipse() +
  labs(title ="flow x s.control", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))

ggplot(pcm.data2, aes(x=f.control, y=structure)) + 
  geom_point() +
  stat_ellipse() +
  labs(title ="f.control x structure", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))

ggplot(pcm.data2, aes(x=flow, y=structure)) + 
  geom_point() +
  stat_ellipse() +
  labs(title ="flow x structure", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))


ggplot(pcm.data2, aes(x=s.control, y=structure)) + 
  geom_point() +
  stat_ellipse() +
  labs(title="s.control x structure", x = "", y = "") +
  xlim(13, 19) +
  ylim(13, 19) +
  theme(axis.text=element_text(size=25),
        title=element_text(size=30))


install.packages("psych")
library(psych)
Scatter_Matrix <- pairs.panels(pcm.data2[, c(1:4)], main = "Scatter Plot Matrix for Structure and Flow Experiments") 

pairs.panels(pcm.data2[, c(1:4)], 
             method = "pearson", # correlation method
             smooth = FALSE,
             rug = FALSE,
             scale = TRUE,
             hist.col = FALSE,
             density = FALSE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )




ggcorrplot(pcm.data.2[,4:1],
           hc.order = TRUE,
           type = "lower",
           outline.color = "white",
           lab = TRUE)
library(data.table)
cordt <- as.data.table(pcm.data2, keep.rownames = 'col_name')
cordt <- melt(cordt, id.vars = 'col_name', variable.name = 'row_name')

# convert to factor so that rows and columns have the same order as the data
cordt[, row_name := factor(row_name, levels = rev(rownames(pcm.data2)))]
cordt[, col_name := factor(col_name, levels = rownames(pcm.data2))]

# set diagonal and the top-right half of the matrix to 0 so that those cells appears white
cordt[ncol(pcm.data2) - as.integer(row_name) < as.integer(col_name), value := 0]
# remove the last column and the bottom row (where left cells are self correlations only)
cordt <- cordt[as.integer(row_name) < ncol(pcm.data2) &
                 as.integer(col_name) < ncol(pcm.data2)]

ggplot(cordt, aes(x = col_name, y = row_name, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high = 'red') +
  labs(x = NULL, y = NULL, fill = 'Corr') +
  theme_minimal()


# Percent Variation ----
f.control.flow <- 0.9178571^2
# 0.8424617
f.control.s.control <- 0.8642857^2
# 0.7469898
f.control.structure <- 0.6642857^2
# 0.4412755
flow.s.control <- 0.7607143^2
# 0.5786862
flow.structure <- 0.6964286^2
# 0.4850128
s.control.structure <- 0.7928571^2
# 0.6286224
