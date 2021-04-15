
# New Western data sent by Rahel on 4/1/21: These data are normalized to a total protein stain instead of actin and uses a new program for the analysis.
# data: 2021_04_01 empiria TIDES Snyder Subset Signal Intensity.xlsx

# ICCs for BP+CF and CS
# Correlations of the normalized protein only 
#     protein x protein
#     protein x normalized transcripts 
# LME associations with 
#     p x t
#     hospital record data
#     phthalates

library(PerformanceAnalytics)
library(ggpubr)
library(dplyr)
library(stargazer)
library(tidyverse)

m <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)

m$Sample <-as.factor(m$Sample)
m$Replicate <-as.factor(m$Replicate)
m$Duplicate <-as.factor(m$Duplicate)

mo <- aggregate(PPARg ~ Duplicate + Sample + Tissue, m, mean) # repeat for all 3 proteins

write.csv(mo, file = "mo.csv")

# subset data by tissue type
mobp <- mo[which(mo$Tissue == "BP"), c(1:6)]
mocf <- mo[which(mo$Tissue == "CF"), c(1:6)]
mocs <- mo[which(mo$Tissue == "CS"), c(1:6)]

# ICCs
# sample var/(residual var + sample var)
n <- lmer(hCGb~Replicate + (1|Sample), data = mo)
summary(n)

# PPARg
# BP: 11357/(6343+11357) = 0.6416384
# CF: 23220/(23220+15492) = 0.599814
# CS: 501560/(501560+88244) = 0.8503842

# hCGb
# BP: 26279268/(26279268+1489387) = 0.9463645
# CF: 10491846/(10491846+814042) = 0.9279984
# CS: 927628/(927628+352606) = 0.7245769

# hCGa
# BP: 5325246/(5325246+528472) = 0.9097203
# CF: 2015876/(2015876+213437) = 0.9042588
# CS: 0/(0+64185) = NA

# collapse BP and CF replicates to give n of 24 instead of stacking and getting n of 48 
mobpcf <- rbind(mobp, mocf)
em <- aggregate(hCGb ~ Replicate + Sample, mobpcf, mean) # repeat for all 3 proteins

# 15612/(15612+3002) = 0.8387235 # PPARg
# 3265518/(3265518+296794) = 0.916685 # hCGb
# 16376928/(16376928+772025) = 0.9549812 # hCGa

# matrix showing bivariate scatterplots, distributions, and correlation values with significance levels
# but doesn't give CIs
c <- mo[which(mo$Tissue == 'BP'), c(4:6)]
chart.Correlation(c, histogram = TRUE, pch = 19, method = "spearman")
mtext("Tides BP", side = 3, line = 3)

# function to obtain spearman correlations with CIs instead
spearman.test <- function(x, y, conf.level = 0.95) {
  RIN <- function(x){qnorm((rank(x) - 0.5)/(length(rank(x))))}
  x_rin <- RIN(x)
  y_rin <- RIN(y)
  list(cor.test(x,y, method='spearman'),
       'RIN corrected CI'= cor.test(x_rin,y_rin)$conf.int[1:2]
  )
}

# apply function (repeat for all tissue types)
spearman.test(emobp$PPARg, emobp$hCGb, conf.level = 0.95)
spearman.test(emobp$PPARg, emobp$hCGa, conf.level = 0.95)
spearman.test(emobp$hCGa, emobp$hCGb, conf.level = 0.95)

# protein x protein correlation results
# PPARg x hCGb (bp, cf, cs): 
# PPARg x hCGa (bp, cf, cs):
# hCGa x hCGb (bp, cf, cs):

# -0.05 (-0.39, 0.42) | -0.34 (-0.69, 0.0009)	| -0.24 (-0.5, 0.3)
# -0.06 (-0.39, 0.42) |	-0.26 (-0.65, 0.08)	  | 0.07 (-0.25, 0.53)
# 0.94 (0.84, 0.97)	  | 0.91 (0.84, 0.97)	    | 0.31 (-0.02, 0.68)

# protein x transcript correlations

moo <- merge(mo, em, by = c("Sample", "Tissue", "Replicate"))
moo <- moo %>% dplyr::rename(PPARg.m = PPARg.y, hCGb.m = hCGb, hCGa.m = hCGa) # new name = old name

spearman.test(cs$CGA.n, cs$hCGa.m, conf.level = 0.95)
spearman.test(cs$CGA.n, cs$hCGb.m, conf.level = 0.95)
spearman.test(cs$CGA.n, cs$PPARg.m, conf.level = 0.95)

spearman.test(cs$CGB.n, cs$hCGa.m, conf.level = 0.95)
spearman.test(cs$CGB.n, cs$hCGb.m, conf.level = 0.95)
spearman.test(cs$CGB.n, cs$PPARg.m, conf.level = 0.95)

spearman.test(cs$PPARg.SDHA, cs$hCGa.m, conf.level = 0.95)
spearman.test(cs$PPARg.SDHA, cs$hCGb.m, conf.level = 0.95)
spearman.test(cs$PPARg.SDHA, cs$PPARg.m, conf.level = 0.95)

spearman.test(cs$SNORA46.SDHA, cs$hCGa.m, conf.level = 0.95)
spearman.test(cs$SNORA46.SDHA, cs$hCGb.m, conf.level = 0.95)
spearman.test(cs$SNORA46.SDHA, cs$PPARg.m, conf.level = 0.95)

# BP	CF	CS
# CGA x hCGa	0.02 (-0.32, 0.48)	0.26 ( -0.19, 0.58)	0.02 (-0.37, 0.47)
# CGA x hCGb	0.05 (-0.28, 0.51)	0.3 (-0.09, 0.64)	0.21 (-0.26, 0.56)
# CGA x PPARg	0.11 (-0.3, 0.5)	0.05 (-0.42, 0.39)	0.22 (-0.21, 0.6)
# CGB x hCGa	0.40 (0.09, 0.74)	0.54 (0.086, 0.73)	0.08 (-0.34, 0.5)
# CGB x hCGb 	0.42 (0.1, 0.74)	0.56 (0.17, 0.77)	0.23 (-0.25, 0.56)
# CGB x PPARg 	0.31 (-0.12, 0.63)	0.07 (-0.4, 0.41)	0.14 (-0.38, 0.46)
# PPARg x hCGa	0.13 (-0.33, 0.47)	-0.21 (-0.53, 0.25)	0.24 (-0.13, 0.64)
# PPARg x hCGb	0.07 (-0.34, 0.46)	-0.26 (-0.57, 0.2)	0.19 (-0.19, 0.61)
# PPARg x PPARg	0.23 (-0.23, 0.55)	0.31 (-0.05, 0.67)	0.58 (0.31, 0.84)
# SNORA46 x hCGa	0.17 (-0.27, 0.52)	-0.1 (-0.48, 0.32)	0.42 (-0.05, 0.69)
# SNORA46 x hCGb 	0.18 (-0.26, 0.53)	0.04 (-0.41, 0.4)	0.08 (-0.34, 0.49)
# SNORA46 x PPARg	-0.19 (-0.56, 0.21)	-0.18 (-0.53, 0.26)	0.34 (0.04, 0.73)

# linear mixed effects model protein x transcript 

hist(moo$PPARg.m)
shapiro.test(moo$hCGb.m)

moo$hCGb.ml = log(moo$hCGb.m)
moo$hCGa.ml = log(moo$hCGa.m)
moo$PPARg.ml = log(moo$PPARg.m)

# check for tissue type interactions with global p-value
p <- lmer(hCGa.ml ~ SNORA46.SDHA:Tissue + (1|Sample) + (1|Tissue), data = moo)

tbl_regression(p, exponentiate = F) %>% 
  add_global_p() %>%
  mold_p () %>%
  mold_labels() %>%
  italicize_levels()

# loop with remaining proteins

a <- lmer(PPARg.ml ~ CGA.n + CGA.n:Tissue + (1|Sample) + (1|Tissue), moo)
b <- lmer(PPARg.ml ~ CGB.n + CGB.n:Tissue + (1|Sample) + (1|Tissue), moo)
c <- lmer(PPARg.ml ~ PPARg.SDHA + PPARg.SDHA:Tissue + (1|Sample) + (1|Tissue), moo)
d <- lmer(PPARg.ml ~ SNORA46.SDHA + SNORA46.SDHA:Tissue + (1|Sample) + (1|Tissue), moo)

class(a) <- "lmerMod"
class(b) <- "lmerMod"
class(c) <- "lmerMod"
class(d) <- "lmerMod"

stargazer(a,b,c,d,
          type="html",
          out="tp.html",
          intercept.mottom = F,
          intercept.top = T,
          ci = T, digits=2,
          model.names = T,
          single.row = T)

# repeat below for phthalates/hr data. Saved in TIDES_emp_4.07.Rmd

# outcome
out_start=7
out_end= 15
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=21
exp_end=58
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)
number=1

for (i in out_start:out_end){
  outcome = colnames(boo)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(boo)[j]
    model <- lmer(get(outcome) ~ get(exposure) + (1|Tissue) + (1|Sample),
                  na.action = na.exclude,
                  data=boo)
    v <- v(model, useScale = FALSE)
    beta <- fixef(model)
    se <- sqrt(diag(v))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    out_beta[number] = as.numeric(beta[2])
    out_se[number] = as.numeric(se[2])
    out_pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    number = number + 1
    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }
}

outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)

outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue,
    obs = out_nobs
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue,
    obs = exp_nobs
  )
all = rbind(outcome, exposure)
all = na.omit(all)

head(all)

data = all %>% 
  mutate(
    type = substr(variable, 1, 2)
  ) %>% 
  spread(type, variable) %>% 
  rename(
    d = dx,
    i = ix
  ) %>% 
  mutate (
    beta = round(beta, 5),
    se = round(se, 5),
    pvalue = round(pvalue, 5)
  ) %>% 
  select(d, i, beta, se, pvalue)

