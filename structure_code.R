# STRUCTURE ----
# R code for looking at variance difference between length among families in control vs treatment tanks
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readr)


### reading the data into R ----
AAHL.dat.fam <- AAHL_data_forR <- read_csv("/Users/mimiobley/Google_Drive/Graduate/MS/FISH SCIENCE/Data Analysis/Structure/AAHL_data_forR.csv")
View(AAHL_data_forR)

AAHL_data_forR <- na.omit(AAHL_data_forR)
AAHL.dat.fam <- na.omit(AAHL.dat.fam)
s.zzq=as.data.frame(table(AAHL.dat.fam[c("treatment", "tank","family")]))
s.zzq=subset(s.zzq, Freq>0)
# shows how many fish were in each family in each tank
s.zzq




# Before adjusting for tank effect ----
## Subset structure and control tanks ----
structure= subset(AAHL.dat.fam, treatment == 'structure')
s.control= subset(AAHL.dat.fam, treatment == 's.control')

## family means (NOT adjusted) ----
s.fam.mean=aggregate(AAHL.dat.fam[c("length")],
                     by= AAHL.dat.fam[c("family","treatment")], mean)

r1 = subset(AAHL.dat.fam, tank == 'R1')
r2 = subset(AAHL.dat.fam, tank == 'R2')
r3 = subset(AAHL.dat.fam, tank == 'R3')
r4 = subset(AAHL.dat.fam, tank == 'R4')
t1 = subset(AAHL.dat.fam, tank == 'T1')
t2 = subset(AAHL.dat.fam, tank == 'T2')
t3 = subset(AAHL.dat.fam, tank == 'T3')
t4 = subset(AAHL.dat.fam, tank == 'T4')

## tank ICCs ----
require(ICC)

ICCest(family, length, data = r1, CI.type = "S")
ICCest(family, length, data = r2, CI.type = "S")
ICCest(family, length, data = r3, CI.type = "S")
ICCest(family, length, data = r4, CI.type = "S")
ICCest(family, length, data = t1, CI.type = "S")
ICCest(family, length, data = t2, CI.type = "S")
ICCest(family, length, data = t3, CI.type = "S")
ICCest(family, length, data = t4, CI.type = "S")


### Create tank ICC dataframe ----
s.tank <- c('r1','r2','r3','r4','t1','t2','t3','t4')
s.treatment <- c('control', 'control', 'structure', 'control', 'structure', 'structure', 'control', 'structure')
s.ICC <- c(0.2016267,0.09377938,0.07035455,0.2200658,0.1504637,0.1338225,0.2017065,0.1007435)
s.varw <- c(5.4555,3.522996,5.738956,3.417962,2.574569,4.279555,3.015139,3.960443)
s.vara <- c(1.377769,0.3645739,0.4343179,0.9644105,0.4559889,0.6611817,0.7618414,0.4436875)


s.ICCdf <- data.frame(s.tank, s.treatment, s.ICC, s.varw, s.vara)
s.ICCdf

### Plot tank ICCs ----
require(ggplot2)
require(ggrepel)

s.ICCplot.tank <- ggplot(data = s.ICCdf, aes(x = s.treatment, y = s.ICC)) +
  geom_point(aes(color=s.treatment),size=3.5, 
             position = position_dodge2(w = 0.15)) +
  scale_x_discrete(expand=c(0.1, 0.5)) + 
  theme(panel.grid = element_line()) + 
  ggtitle("ICC comparison between control\n and structure tanks") + 
  labs(x="Treatment", y = "ICC") + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15))
s.ICCplot.tank

### geom_label_repel
s.ICCplot.tank + 
  geom_label_repel(aes(label = s.tank),
                   box.padding   = .4,
                   point.padding = .1,
                   min.segment.length = Inf,
                   segment.color = 'grey50') 


## perform Welch's t-test ----
t.test(s.ICC ~ s.treatment, s.ICCdf)

## Family plot ----
s.famplot.na <- ggplot(data = s.fam.mean, aes(x = treatment, y = length, label = family, group = family)) +
  geom_point() + geom_line(aes(color=family)) + 
  scale_y_continuous(limits=c(13,18)) + 
  theme(panel.grid = element_line()) + 
  labs(title="Mean family length in control and structure") + 
  theme(plot.title = element_text(hjust = 0.5)) 
s.famplot.na

s.famplot.na + 
  geom_label_repel(data = subset(s.fam.mean, treatment %in% c("s.control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(s.fam.mean, treatment %in% c("structure")),
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
s.tankmean <- aggregate(AAHL.dat.fam[c("length")],
                        by= AAHL.dat.fam[c("tank","treatment")], mean)

## find mean treatment length  ----
s.treatmentmean <- aggregate(AAHL.dat.fam[c("length")],
                             by= AAHL.dat.fam[c("treatment")], mean)
# treatment   length
# 1   Control 15.48447
# 2 Structure 15.75731

### find variance  ----
s.tankmean$var <- c((15.67659-15.48447), (15.39800-15.48447), (15.24010-15.48447), (15.62492-15.48447), (15.43973-15.75731), (15.95804-15.75731), (15.85069-15.75731), (15.78360-15.75731))
# tank treatment   length      var
# 1   r1   Control 15.67659  0.19212
# 2   r2   Control 15.39800 -0.08647
# 3   r4   Control 15.24010 -0.24437
# 4   t3   Control 15.62492  0.14045
# 5   r3 Structure 15.43973 -0.31758
# 6   t1 Structure 15.95804  0.20073
# 7   t2 Structure 15.85069  0.09338
# 8   t4 Structure 15.78360  0.02629

### adjust length ----
r1$adjusted <-c(t(r1[, 1] - 0.19212))
r2$adjusted <-c(t(r2[, 1] - -0.08647))
r4$adjusted <-c(t(r4[, 1] - -0.24437))
t3$adjusted <-c(t(t3[, 1] - 0.14045))
r3$adjusted <-c(t(r3[, 1] - -0.31758))
t1$adjusted <-c(t(t1[, 1] - 0.20073))
t2$adjusted <-c(t(t2[, 1] - 0.09338))
t4$adjusted <-c(t(t4[, 1] - 0.02629))

### add to AAHL.dat.fam  ----
s.adjusted = c(r1$adjusted, r2$adjusted, r3$adjusted, r4$adjusted, t1$adjusted, t2$adjusted, t3$adjusted, t4$adjusted)
AAHL.dat.fam['adjusted'] <- s.adjusted
AAHL.dat.fam

## re-subset structure and control tanks to include adjusted values ----
structure= subset(AAHL.dat.fam, treatment == 'structure')
s.control= subset(AAHL.dat.fam, treatment == 's.control')

## adjusted family means  ----
s.fam.mean.adjusted=aggregate(AAHL.dat.fam[c("adjusted")],
                              by= AAHL.dat.fam[c("family","treatment")], mean)



write.csv(AAHL.dat.fam, file = "AAHL.dat.fam.csv")
# After adjusting for tank effect  ----

## treatment ICCS (adjusted) ----
require(ICC)
ICCest(family, adjusted, data = structure, CI.type = "S")
# ICC = 0.0800958
# varw = 4.250114
# vara = 0.3700562
# upper CI = 0.1439589
# lower CI = 0.01623266
ICCest(family, adjusted, data = s.control, CI.type = "S")
# ICC = 0.172281
# varw = 3.861283
# vara = 0.8036853
# upper CI = 0.2859454
# lower CI = 0.05861654

# Create ICC dataframe
s.treatment2 <- c('structure', 'control')
s.ICC2 <- c(0.0800958,0.172281)
s.varw2 <- c(4.250114, 3.861283)
s.vara2 <- c(0.3700562, 0.8036853)
s.upperCI2 <- c(0.1439589, 0.2859454)
s.lowerCI2 <- c(0.01623266, 0.05861654)

s.ICCdf2 <- data.frame(s.treatment2, s.ICC2, s.varw2, s.vara2, s.upperCI2, s.lowerCI2)
s.ICCdf2

### Plot treatment ICCs (adjusted) ----
require(ggplot2)
# with error bars
s.ICCplot.treat <- ggplot(data = s.ICCdf2, aes(x = s.treatment2, y = s.ICC2)) +
  geom_point(aes(color=s.treatment2),size=3.5, position = position_dodge2(w = 0.15)) +
  scale_x_discrete(expand=c(0.1, 0.5)) + 
  scale_y_continuous(name="ICC", limits=c(0, .4)) +
  theme(panel.grid = element_line()) + 
  ggtitle("ICC comparison between\n Structure and Control (adjusted)") + 
  labs(x="Treatment", y = "ICC") + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15)) +
  geom_errorbar(
    aes(x=s.treatment2, 
        ymin = s.lowerCI2, 
        ymax = s.upperCI2), width=0.1, 
    color = "red"
  )
s.ICCplot.treat

## perform Welch's t-test ----
t.test(s.ICC2 ~ s.treatment2, s.ICCdf2)


## CV and opportunity for selection ----
# fish
# control cv

sc.cv.fish <- sd(s.control$adjusted) / mean(s.control$adjusted) * 100
sc.cv.fish # 13.86672

# control ofs
sc.ofs.fish <- ((sc.cv.fish / 100) ^2)
sc.ofs.fish # 0.01922861

# structure cv
ss.cv.fish <- ((.86) * 100) / mean(structure$adjusted) 
ss.cv.fish # 13.60368

# structure ofs
ss.ofs.fish <- ((ss.cv.fish / 100) ^2)
ss.ofs.fish # 0.01850602

# fam mean
# subset fam means for f.control and structure
s.control.fam.means<-subset(s.fam.mean.adjusted, treatment!="structure")
structure.fam.means<-subset(s.fam.mean.adjusted, treatment!="s.control")

# control cv
sc.cv.fam <- sd(s.control.fam.means$adjusted) / mean(s.control.fam.means$adjusted) * 100
sc.cv.fam # 5.913798

# control ofs
sc.ofs.fam <- ((sc.cv.fam / 100) ^2)
sc.ofs.fam # 0.0034973

# structure cv
ss.cv.fam <- sd(structure.fam.means$adjusted) / mean(structure.fam.means$adjusted) * 100
ss.cv.fam # 4.222503

# structure ofs
ss.ofs.fam <- ((ss.cv.fam / 100) ^2)
ss.ofs.fam # 0.001782953


## Family plots (adjusted) ----
# plot family mean lengths comparing control and structure (adjusted)
s.famplot.a <- ggplot(data = s.fam.mean.adjusted, 
                      aes(x = treatment, 
                          y = adjusted, 
                          label = family, 
                          group = family)) + 
  geom_point() + geom_line(aes(color=family)) + 
  scale_y_continuous(limits=c(13,18)) + 
  theme(panel.grid = element_line()) + 
  labs(title="Mean family length in\n Control and Structure (adjusted)")  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y= "Length", 
       x = "Treatment") +
  labs(colour = "Family") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20)) +
  scale_x_discrete(labels=c('Control', 'Structure'))
s.famplot.a

s.famplot.a + 
  geom_label_repel(data = subset(s.fam.mean.adjusted, treatment %in% c("s.control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(s.fam.mean.adjusted, treatment %in% c("structure")),
                   nudge_x       = .2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50')



          ## Example
ex.famplot <- ggplot(data = t_example, 
                      aes(x = Group, 
                          y = Length, 
                          label = Family, 
                          group = Family)) + 
  geom_point() + geom_line(aes(color=Family)) + 
  scale_y_continuous(limits=c(13,19)) + 
  theme(panel.grid = element_line()) + 
  labs(title="Example Mean family length in\n Control and Treatment")  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y= "Length", 
       x = "Treatment") +
  labs(colour = "Family") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20)) +
  scale_x_discrete(labels=c('Control', 'Treatment'))
ex.famplot

ex.famplot + 
  geom_label_repel(data = subset(t_example, Group %in% c("Control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(example, Group %in% c("Treatment")),
                   nudge_x       = .2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50')



# MIXED MODEL ----
mods1 <- lmer(length ~ 1 + treatment + (1+treatment|family) + (1|treatment/tank), AAHL.dat.fam, REML = F)
summary(mods1)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: length ~ 1 + treatment + (1 + treatment | family) + (1 | treatment/tank)
#    Data: AAHL.dat.fam
# 
#      AIC      BIC   logLik deviance df.resid 
#  10450.4  10496.8  -5217.2  10434.4     2437 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.8761 -0.4172  0.1554  0.6518  2.1712 
# 
# Random effects:
#  Groups         Name               Variance Std.Dev. Corr 
#  family         (Intercept)        0.74326  0.8621        
#                 treatmentstructure 0.29887  0.5467   -0.73
#  tank:treatment (Intercept)        0.03596  0.1896        
#  treatment      (Intercept)        0.00000  0.0000        
#  Residual                          4.05823  2.0145        
# Number of obs: 2445, groups:  family, 15; tank:treatment, 8; treatment, 2
# 
# Fixed effects:
#                    Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)         15.4739     0.2487 18.1211  62.219   <2e-16 ***
# treatmentstructure   0.2589     0.2115 14.4718   1.224    0.241    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# trtmntstrct -0.667
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')
anova(mods1)
rand(mods1)



mods2 <- lmer(length ~ 1 + treatment + (1+treatment|family) + (1|tank), AAHL.dat.fam, REML = F)
summary(mods2)
# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: length ~ 1 + treatment + (1 + treatment | family) + (1 | tank)
#    Data: AAHL.dat.fam
# 
#      AIC      BIC   logLik deviance df.resid 
#  10448.4  10489.0  -5217.2  10434.4     2438 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.8761 -0.4172  0.1554  0.6518  2.1712 
# 
# Random effects:
#  Groups   Name               Variance Std.Dev. Corr 
#  family   (Intercept)        0.74322  0.8621        
#           treatmentstructure 0.29887  0.5467   -0.73
#  tank     (Intercept)        0.03596  0.1896        
#  Residual                    4.05823  2.0145        
# Number of obs: 2445, groups:  family, 15; tank, 8
# 
# Fixed effects:
#                    Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)         15.4739     0.2487 18.1223  62.220   <2e-16 ***
# treatmentstructure   0.2589     0.2115 14.4722   1.224    0.241    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr)
# trtmntstrct -0.667

# Compare AIC values to determine if "tank" should be nested:
# modf1 = 10450.4
# modf2 = 10448.4  <--- the better model without nested



mods3 <- lmer(length ~ 1 + treatment*family + (1|tank), AAHL.dat.fam, REML = F)
summary(mods3)
anova(mods3)



### F test treat ----
s.ftest.treat <- var.test(adjusted ~ treatment, data = AAHL.dat.fam)
s.ftest.treat
# F test to compare two variances
# 
# data:  adjusted by treatment
# F = 1.0034, num df = 1261, denom df = 1182, p-value = 0.9531
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.896758 1.122513
# sample estimates:
#   ratio of variances 
# 1.003409

#### variances ----
s.control.var.treat <- var(subset(AAHL.dat.fam, treatment == 's.control', adjusted))
s.control.var.treat
# 4.61042

structure.var.treat <- var(subset(AAHL.dat.fam, treatment == 'structure', adjusted))
structure.var.treat
# 4.594755


### F test fam ----
s.ftest.fam <- var.test(adjusted ~ treatment, data = s.fam.mean.adjusted)
s.ftest.fam
# F test to compare two variances
# 
# data:  adjusted by treatment
# F = 1.8985, num df = 14, denom df = 14, p-value = 0.2426
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.6373828 5.6548494
# sample estimates:
#   ratio of variances 
# 1.8985 

#### variances ----
s.control.var.fam <- var(subset(s.fam.mean.adjusted, treatment == 's.control', adjusted))
s.control.var.fam
# 0.8372986

structure.var.fam <- var(subset(s.fam.mean.adjusted, treatment == 'structure', adjusted))
structure.var.fam
# 0.4410316


### without family E ----
no.e<-subset(AAHL.dat.fam, family!="E")
no.e

# fit nested ANOVA with interaction term (family*treatment) without family E
ANOVA_no.e <- aov(length ~ treatment + family + (treatment/factor(tank)) + (family*treatment), data = no.e)
# view summary of nested ANOVA
summary(ANOVA_no.e)



#### Family plots without family E ----
## adjusted family means from no.e dataset
s.fam.mean.a.no.e=aggregate(no.e[c("adjusted")],
                              by= no.e[c("family","treatment")], mean)

# plot family mean lengths comparing control and structure (adjusted)
s.famplot.a.no.e <- ggplot(data = s.fam.mean.adjusted.no.e, aes(x = treatment, y = adjusted, label = family, group = family)) +
  geom_point() + 
  geom_line(aes(color=family)) + 
  scale_y_continuous(limits=c(13,18)) + 
  theme(panel.grid = element_line()) + 
  labs(title="Mean family length in control and structure (adjusted).no.e") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y= "Length", x = "Treatment")
s.famplot.a.no.e

s.famplot.a.no.e + 
  geom_label_repel(data = subset(s.fam.mean.adjusted.no.e, treatment %in% c("s.control")),
                   nudge_x       = -.2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  geom_label_repel(data = subset(s.fam.mean.adjusted.no.e, treatment %in% c("structure")),
                   nudge_x       = .2,
                   direction     = "y",
                   hjust         = 1,
                   segment.size  = 0.2,
                   box.padding   = 0.2, 
                   point.padding = 0.1,
                   segment.color = 'grey50')

#### ICCtank without fam E ----
# subset tanks from dataset without fam E
r1.no.e = subset(no.e, tank == 'R1')
r2.no.e = subset(no.e, tank == 'R2')
r3.no.e = subset(no.e, tank == 'R3')
r4.no.e = subset(no.e, tank == 'R4')
t1.no.e = subset(no.e, tank == 'T1')
t2.no.e = subset(no.e, tank == 'T2')
t3.no.e = subset(no.e, tank == 'T3')
t4.no.e = subset(no.e, tank == 'T4')

## tank ICCs 
require(ICC)

ICCest(family, length, data = r1.no.e, CI.type = "S")
ICCest(family, length, data = r2.no.e, CI.type = "S")
ICCest(family, length, data = r3.no.e, CI.type = "S")
ICCest(family, length, data = r4.no.e, CI.type = "S")
ICCest(family, length, data = t1.no.e, CI.type = "S")
ICCest(family, length, data = t2.no.e, CI.type = "S")
ICCest(family, length, data = t3.no.e, CI.type = "S")
ICCest(family, length, data = t4.no.e, CI.type = "S")


### Create tank ICC dataframe 
s.tank.no.e <- c('r1','r2','r3','r4','t1','t2','t3','t4')
s.treatment.no.e <- c('control', 'control', 'structure', 'control', 'structure', 'structure', 'control', 'structure')
s.ICC.no.e <- c(0.104768, 0.07822021, 0.07276745, 0.1753236, 0.1504637, 0.1375568, 0.1984511, 0.1099726)
s.varw.no.e <- c(4.953575, 3.421729, 5.553356, 2.970968, 2.374316, 4.425261, 2.753218, 3.907608)
s.vara.no.e <- c(0.5797113, 0.2903604, 0.4358169, 0.6415218, 0.5047719, 0.7058143, 0.6816539, 0.4828278)


s.ICCdf.no.e <- data.frame(s.tank.no.e, s.treatment.no.e, s.ICC.no.e, s.varw.no.e, s.vara.no.e)
s.ICCdf.no.e

### Plot tank ICCs 
require(ggplot2)
require(ggrepel)

s.ICCplot.tank.no.e <- ggplot(data = s.ICCdf.no.e, aes(x = s.treatment.no.e, y = s.ICC.no.e)) +
  geom_point(aes(color=s.treatment.no.e),size=3.5, 
             position = position_dodge2(w = 0.15))+
  scale_x_discrete(expand=c(0.1, 0.5)) + theme(panel.grid = element_line()) + ggtitle("ICC comparison between control and structure tanks.no.e")
s.ICCplot.tank.no.e

### geom_label_repel
s.ICCplot.tank.no.e + 
  geom_label_repel(aes(label = s.tank.no.e),
                   box.padding   = .4,
                   point.padding = .1,
                   min.segment.length = Inf,
                   segment.color = 'grey50') 


## perform Welch's t-test 
t.test(s.ICC.no.e ~ s.treatment.no.e, s.ICCdf.no.e)




#### ICCtreatment without fam E ----
# subset treatments from dataset without fam E
structure.no.e= subset(no.e, treatment == 'structure')
s.control.no.e= subset(no.e, treatment == 's.control')

## treatment ICCS (adjusted)
require(ICC)
ICCest(family, adjusted, data = structure.no.e, CI.type = "S")
# ICC = 0.0864179
# varw = 4.179328
# vara = 0.3953325
# upper CI = 0.1567756
# lower CI = 0.01606021
ICCest(family, adjusted, data = s.control.no.e, CI.type = "S")
# ICC = 0.1274566
# varw = 3.516351
# vara = 0.5136501
# upper CI = 0.2215751
# lower CI = 0.033338

# Create ICC dataframe
s.treatment2.no.e <- c('structure', 'control')
s.ICC2.no.e <- c(0.0864179,0.1274566)
s.varw2.no.e <- c(4.179328, 3.516351)
s.vara2.no.e <- c(0.3953325, 0.5136501)
s.upperCI2.no.e <- c(0.1567756, 0.2215751)
s.lowerCI2.no.e <- c(0.01606021, 0.033338)

s.ICCdf2.no.e <- data.frame(s.treatment2.no.e, s.ICC2.no.e, s.varw2.no.e, s.vara2.no.e, s.upperCI2.no.e, s.lowerCI2.no.e)
s.ICCdf2.no.e

### Plot treatment ICCs (adjusted) 
require(ggplot2)
# with error bars
s.ICCplot.treat.no.e <- ggplot(data = s.ICCdf2.no.e, aes(x = s.treatment2.no.e, y = s.ICC2.no.e)) +
  geom_point(aes(color=s.treatment2.no.e),size=3.5, position = position_dodge2(w = 0.15)) +
  scale_x_discrete(expand=c(0.1, 0.5)) + 
  scale_y_continuous(name="ICC", limits=c(0, .4)) +
  theme(panel.grid = element_line()) + 
  ggtitle("ICC comparison between Structure and Control (adjusted).no.e") + 
  geom_errorbar(
    aes(x=s.treatment2, 
        ymin = s.lowerCI2.no.e, 
        ymax = s.upperCI2.no.e), width=0.1, 
    color = "red"
  )
s.ICCplot.treat.no.e




# Boxplots ----
# Create data with reordered group levels
AAHL.dat.fam_ordered <- AAHL.dat.fam  
# Order boxes by median
AAHL.dat.fam_ordered.s <- with(AAHL.dat.fam,reorder(family, length, median))
AAHL.dat.fam_ordered$family <- factor(AAHL.dat.fam_ordered$family,levels = levels(AAHL.dat.fam_ordered.s))

# histogram combined not adjusted
ggplot(AAHL.dat.fam, aes(x=length, color=treatment)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + labs(x="Length (cm)", y = "Count") +
  ggtitle("Lengths for Structure and Control\n Distribution") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15)) 

# histogram combined adjusted
ggplot(AAHL.dat.fam, aes(x=adjusted, color=treatment)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + labs(x="Length (cm)", y = "Count") +
  ggtitle("Lengths for Structure and Control\n Distribution (adjusted)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(colour = "Treatment") +
  expand_limits(x = 0, y = 0) +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15)) 


# not adjusted
ggplot(AAHL.dat.fam_ordered, aes(x = family, y = length, fill = treatment)) + 
  geom_boxplot() +
  labs(title = "Variability of Structure vs. Control for Each Family", 
       x = "family", 
       y = "length") +
  theme_minimal()

# adjusted
ggplot(AAHL.dat.fam_ordered, aes(x = family, y = adjusted, fill = treatment)) + 
  geom_boxplot() +
  labs(title = "Variability of Structure vs. Control\n for Each Family (adjusted)", 
       x = "Family", 
       y = "Length (cm)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() + 
  labs(colour = "Treatment") +
  theme(axis.text=element_text(size=15),
        title=element_text(size=20),
        legend.text=element_text(size=15)) 



## control (NOT adjusted) ----

# Create data with reordered group levels
s.control_ordered <- s.control    
# Order boxes by median
fam_ordered.cs <- with(s.control,reorder(family, length, median))
s.control_ordered$family <- factor(s.control_ordered$family,
                                   levels = levels(fam_ordered.cs))

# Draw ggplot2 boxplot ordered by median
box.s.control.na <-  ggplot(s.control_ordered,                    
                           aes(x = family,
                               y = length)) +
                      geom_boxplot() + 
                      labs(y= "Length", x = "Family") + 
                      ggtitle("Lengths for Control") + 
                      theme(plot.title = element_text(hjust = 0.5))
box.s.control.na

## structure (NOT adjusted) ----
# Create data with reordered group levels
structure_ordered <- structure    
# Order boxes by median
fam_ordered.ss <- with(structure,reorder(family, length, median))
structure_ordered$family <- factor(structure_ordered$family,
                                   levels = levels(fam_ordered.ss))

# Draw ggplot2 boxplot ordered by median
box.strucutre.na <-  ggplot(structure_ordered,                          
                             aes(x = family,
                                 y = length)) +
                        geom_boxplot() + 
                        labs(y= "Length", x = "Family") + 
                        ggtitle("Lengths for Structure") + 
                        theme(plot.title = element_text(hjust = 0.5))
box.strucutre.na

## control (adjusted) ----

# Create data with reordered group levels
s.control_ordered.a <- s.control    
# Order boxes by median
fam_ordered.csa <- with(s.control,reorder(family, adjusted, median))
s.control_ordered.a$family <- factor(s.control_ordered.a$family,
                                     levels = levels(fam_ordered.csa))

# Draw ggplot2 boxplot ordered by median
box.s.control.a <- ggplot(s.control_ordered.a,                   
                               aes(x = family,
                                   y = adjusted)) +
                          geom_boxplot() + 
                          labs(y= "Length (adjusted)", x = "Family") + 
                          ggtitle("Lengths (adjusted) for Control") + 
                          theme(plot.title = element_text(hjust = 0.5))
box.s.control.a

## structure (adjusted)----
# Create data with reordered group levels
structure_ordered.a <- structure    
# Order boxes by median
fam_ordered.ssa <- with(structure,reorder(family, adjusted, median))
structure_ordered.a$family <- factor(structure_ordered.a$family,
                                     levels = levels(fam_ordered.ssa))

# Draw ggplot2 boxplot ordered by median
box.structure.a <- ggplot(structure_ordered.a,                     
                               aes(x = family,
                                   y = adjusted)) +
                          geom_boxplot() + 
                          labs(y= "Length (adjusted)", x = "Family") + 
                          ggtitle("Lengths (adjusted) for Structure") + 
                          theme(plot.title = element_text(hjust = 0.5)) 
box.structure.a


# Histograms ----
## control count (NOT adjusted) ----
hist.s.control.na <- ggplot(s.control_ordered, aes(x=length)) + 
                      geom_histogram(colour="black", fill="white") +
                      labs(x="Length(cm)", y = "Count") +
                      ggtitle("Lengths for Control Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5)) 
hist.s.control.na

## control density (NOT adjusted) ----
hist.s.control.den.na <-  ggplot(s.control_ordered, aes(x=length)) + 
                      geom_histogram(aes(y=after_stat(density)),
                                     colour="black", fill="white") +
                      geom_density(alpha=.2, fill="#FF6666") +
                      labs(x="Length(cm)", y = "Density") +
                      ggtitle("Lengths for Control Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5))
hist.s.control.den.na

## structure count (NOT adjusted) ----
hist.structure.na <-  ggplot(structure_ordered, aes(x=length)) + 
                      geom_histogram(colour="black", fill="white") +
                      labs(x="Length(cm)", y = "Count") +
                      ggtitle("Lengths for Structure Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5)) 
hist.structure.na

## structure density (NOT adjusted) ----
hist.structure.den.na <-  ggplot(structure_ordered, aes(x=length)) + 
                      geom_histogram(aes(y=after_stat(density)), 
                                     colour="black", fill="white")+
                      geom_density(alpha=.2, fill="#FF6666") +
                      labs(x="Length(cm)", y = "Density") +
                      ggtitle("Lengths for Structure Distribution") + 
                      theme(plot.title = element_text(hjust = 0.5))
hist.structure.den.na



## control count (adjusted) ----
hist.s.control.a <- ggplot(s.control_ordered, aes(x=adjusted)) + 
                    geom_histogram(colour="black", fill="white") +
                    labs(x="Length(cm)", y = "Count") +
                    ggtitle("Lengths for Control Distribution (adjusted)") +
                    theme(plot.title = element_text(hjust = 0.5)) 
hist.s.control.a

## control density (adjusted) ----
hist.s.control.den.a <- ggplot(s.control_ordered, aes(x=adjusted)) + 
                    geom_histogram(aes(y=after_stat(density)), 
                                   colour="black", fill="white")+
                    geom_density(alpha=.2, fill="#FF6666") +
                    labs(x="Length(cm)", y = "Density") +
                    ggtitle("Lengths for Control Distribution (adjusted)") +
                    theme(plot.title = element_text(hjust = 0.5))
hist.s.control.den.a

## structure count (adjusted) ----
hist.structure.a <- ggplot(structure_ordered, aes(x=adjusted)) + 
                    geom_histogram(colour="black", fill="white") +
                    labs(x="Length(cm)", y = "Count") +
                    ggtitle("Lengths for Structure Distribution (adjusted)") +
                    theme(plot.title = element_text(hjust = 0.5)) 
hist.structure.a

## structure density (adjusted) ----
hist.structure.den.a <- ggplot(structure_ordered, aes(x=adjusted)) + 
                    geom_histogram(aes(y=after_stat(density)), 
                                   colour="black", fill="white")+
                    geom_density(alpha=.2, fill="#FF6666") +
                    labs(x="Length(cm)", y = "Density") +
                    ggtitle("Lengths for Structure Distribution (adjusted)") +
                    theme(plot.title = element_text(hjust = 0.5))
hist.structure.den.a


# Q-Q plots ----
## (NOT adjusted) ----
qq.s.na <- ggplot(AAHL.dat.fam, aes(sample = length)) +
            stat_qq(aes(color = treatment)) +
            scale_color_manual(values = c("darkgrey", "black"))+
            labs(y = "Length") +
            ggtitle("Q-Q plot for Control and Structure Length") + 
            theme(plot.title = element_text(hjust = 0.5))
qq.s.na

## (adjusted) ----
qq.s.a <- ggplot(AAHL.dat.fam, aes(sample = adjusted)) +
            stat_qq(aes(color = treatment)) +
            scale_color_manual(values = c("darkgrey", "black"))+
            labs(y = "Length") +
            ggtitle("Q-Q plot for Control and Structure Length (adjusted)") + 
            theme(plot.title = element_text(hjust = 0.5))
qq.s.a

# pearson correlation matrix (in flow code) ----

