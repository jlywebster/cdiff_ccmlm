
#-------------------------------------------------------------#
# Code used for manuscript:                                   # 
#                                                             #
# Estimating the impact of patient-level risk factors and     #
# time-varying hospital unit on healthcare-associated         #
# Clostridioides difficile infection using cross-classified   #
# multilevel models                                           #
#                                                             #
# Author: Jessica L. Webster                                  #
# jlywebster@gmail.com                                        #
#                                                             #
# Updated 06.12.2024                                          #
#-------------------------------------------------------------#


# Calling in libraries
library(tidyverse)
library(lme4)
library(broom)
library(broom.mixed)

# Pulling in long dataset
chris_long_a <- read.csv("christiana_long_analytic.csv", header=T, na.strings = c("","NA"))
datlong <- chris_long_a
# Pulling in wide dataset
chris_wide_a <- read.csv("christiana_wide_analytic.csv", header=T, na.strings = c("","NA"))
datwide <- chris_wide_a

# Subsetting to variables I want to include in my model: age (cont), insurance (cat), AHRQ (cont),
# any medications (cat), any antibiotics (cat), any procedures (cat), length of stay (cont)
vars_wide <- c("CdiffPositive","PatientPersonID","AnyAntibiotics","AnyMedications","AnyProcedures","age","insur_cat",
               "AHRQ","los","TotalUnits","longest_unit","season_admission")

vars_long <- c("CdiffPositive","PatientPersonID","AnyDailyAntibiotics","AnyDailyMedications","AnyDailyProcedures",
               "AgeInYears","Insurance_cat","AHRQComorbidityCnt","hosp_day0","cum_units","unit_dept","season_hosp_date")

datlong_vars <- datlong[,vars_long]
datwide_vars <- datwide[,vars_wide]

# some cleaning
datlong_vars <- within(datlong_vars, {
  Insurance_cat <- factor(Insurance_cat)
  AnyDailyAntibiotics <- factor(AnyDailyAntibiotics)
  AnyDailyMedications <- factor(AnyDailyMedications)
  AnyDailyProcedures <- factor(AnyDailyProcedures)
  CdiffPositive <- factor(CdiffPositive)
  season_hosp_date <- factor(season_hosp_date)
})
datwide_vars <- within(datwide_vars, {
  CdiffPositive <- factor(CdiffPositive)
  insur_cat <- factor(insur_cat)
  AnyMedications <- factor(AnyMedications)
  AnyAntibiotics <- factor(AnyAntibiotics)
  AnyProcedures <- factor(AnyProcedures)
  longest_unit <- factor(longest_unit)
  longest_unit <- relevel(longest_unit, ref = "C_medicine")
  season_admission <- factor(season_admission)
})

# rescaling and reshaping numeric variables
datlong$log_hosp_day0 <- log(datlong$hosp_day0)
numcols <- c("AgeInYears", "AHRQComorbidityCnt", "log_hosp_day0","cum_units")
datlong[,numcols] <- scale(datlong[,numcols])

#### Regression analyses ####

## 3.1. Hierarchical ML models
# patients nested within longest unit
model1.empty <- glmer(CdiffPositive ~ (1|longest_unit), 
                      data=datwide_vars, family=binomial(link="logit"));summary(model1.empty)

model1.adj <- glmer(CdiffPositive ~ AnyAntibiotics + AnyMedications + AnyProcedures
                    + age + insur_cat + AHRQ + los + TotalUnits
                    + (1|longest_unit),
                    data=datwide_vars, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model1.adj)
tidy(model1.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

## 3.2. Cross-classified MLM using patient day and patient, unit as covariate
# patient-day-level vars: any antibiotic treatment, any non-antibiotic medications, any procedures, length of stay, and number of wards visited/transfers
# patient-level vars: age at admission, insurance, and AHRQ comorbidities
# ward: time-varying ward indicator 

## empty model
model2.empty <- glmer(CdiffPositive ~ (1|PatientPersonID) + (1|unit_dept), 
                      data=datlong, family=binomial(link="logit"));summary(model2.empty)

model2.adj <- glmer(CdiffPositive ~ AnyDailyAntibiotics + AnyDailyMedications + AnyDailyProcedures 
                    + AgeInYears + Insurance_cat + AHRQComorbidityCnt + cum_units + hosp_day0
                    + (1|PatientPersonID) + (1|unit_dept), 
                    data=datlong, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model2.adj)
tidy(model2.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")


## calculating median odds ratios and confidence intervals
model.list<-list(model1.empty, model1.adj, model2.empty, model2.adj)

model.list_output<-map_dfr(model.list, function(model){
  coefs<-tidy(model) %>% filter(effect=="ran_pars")%>%
    mutate(tau2=estimate^2,
           MOR=exp(0.95*sqrt(tau2))) %>%
    select(group, tau2, MOR)
  coefs
})
model.list_output

## bootstrapped confidence intervals 
# empty hierarchical MLM unit
t1a<- ranef(model1.empty,condVar = TRUE)
est<-as.numeric(unlist(t1a$longest_unit))
var<- as.numeric(unlist(attr(t1a$longest_unit,"postVar")))
### bootstrap 
### create empty output collection:
mor_boot<- c()
#### iterate over replicates
for(i in 1:1000){
  ### draw vector of area random effects from normal
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  ### create data frame with all possible pairs 
  s<- combn(drw,2)
  #### estimate MOR and save in output
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
### bootstrap median and 95% CI
quantile(mor_boot, c(.025,.5,.975))
# adjusted hierarchical MLM unit
t1b<- ranef(model1.adj,condVar = TRUE)
est<-as.numeric(unlist(t1b$longest_unit))
var<- as.numeric(unlist(attr(t1b$longest_unit,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
# empty CC MLM patient
t2a_p<- ranef(model2.empty,condVar = TRUE)
est<-as.numeric(unlist(t2a_p$PatientPersonID))
var<- as.numeric(unlist(attr(t2a_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
# empty CC MLM unit
t2a_u<- ranef(model2.empty,condVar = TRUE)
est<-as.numeric(unlist(t2a_u$unit_dept))
var<- as.numeric(unlist(attr(t2a_u$unit_dept,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
# adjusted CC MLM patient
t2b_p<- ranef(model2.adj,condVar = TRUE)
est<-as.numeric(unlist(t2b_p$PatientPersonID))
var<- as.numeric(unlist(attr(t2b_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
# adjusted CC MLM unit
t2b_u<- ranef(model2.adj,condVar = TRUE)
est<-as.numeric(unlist(t2b_u$unit_dept))
var<- as.numeric(unlist(attr(t2b_u$unit_dept,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))


## 3.3. sensitivity analysis #1 (clustering by season)
# hierarchical MLM
model3.empty <- glmer(CdiffPositive ~ (1|longest_unit) + (1|season_admission), 
                      data=datwide, family=binomial(link="logit"));summary(model3.empty)

model3.adj <- glmer(CdiffPositive ~ AnyAntibiotics + AnyMedications + AnyProcedures
                    + age + insur_cat + AHRQ + los + TotalUnits
                    + (1|longest_unit) + (1|season_admission),
                    data=datwide, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model3.adj)
tidy(model3.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

# CC MLM
model4.empty <- glmer(CdiffPositive ~ (1|PatientPersonID) + (1|unit_dept) + (1|season_hosp_date), 
                      data=datlong, family=binomial(link="logit"));summary(model4.empty)

model4.adj <- glmer(CdiffPositive ~ AnyDailyAntibiotics + AnyDailyMedications + AnyDailyProcedures 
                    + AgeInYears + Insurance_cat + AHRQComorbidityCnt + hosp_day0 + cum_units
                    
                    + (1|PatientPersonID) + (1|unit_dept) + (1|season_hosp_date), 
                    data=datlong, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model4.adj)
tidy(model4.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")


model.list2<-list(model3.empty, model3.adj, model4.empty, model4.adj)

## Median odds ratio
model.list2_output<-map_dfr(model.list2, function(model){
  coefs<-tidy(model) %>% filter(effect=="ran_pars")%>%
    mutate(tau2=estimate^2,
           MOR=exp(0.95*sqrt(tau2))) %>%
    select(group, tau2, MOR)
  coefs
})
model.list2_output

## bootstrapping MOR confidence intervals
t3a_u<- ranef(model3.empty,condVar = TRUE)
est<-as.numeric(unlist(t3a_u$longest_unit))
var<- as.numeric(unlist(attr(t3a_u$longest_unit,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t3a_s<- ranef(model3.empty,condVar = TRUE)
est<-as.numeric(unlist(t3a_s$season_admission))
var<- as.numeric(unlist(attr(t3a_s$season_admission,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))

t3b_u<- ranef(model3.adj,condVar = TRUE)
est<-as.numeric(unlist(t3b_u$longest_unit))
var<- as.numeric(unlist(attr(t3b_u$longest_unit,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t3b_s<- ranef(model3.adj,condVar = TRUE)
est<-as.numeric(unlist(t3b_s$season_admission))
var<- as.numeric(unlist(attr(t3b_s$season_admission,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4a_p<- ranef(model4.empty,condVar = TRUE)
est<-as.numeric(unlist(t4a_p$PatientPersonID))
var<- as.numeric(unlist(attr(t4a_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4a_u<- ranef(model4.empty,condVar = TRUE)
est<-as.numeric(unlist(t4a_u$unit_dept))
var<- as.numeric(unlist(attr(t4a_u$unit_dept,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4a_s<- ranef(model4.empty,condVar = TRUE)
est<-as.numeric(unlist(t4a_s$season_hosp_date))
var<- as.numeric(unlist(attr(t4a_s$season_hosp_date,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4b_p<- ranef(model4.adj,condVar = TRUE)
est<-as.numeric(unlist(t4b_p$PatientPersonID))
var<- as.numeric(unlist(attr(t4b_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4b_u<- ranef(model4.adj,condVar = TRUE)
est<-as.numeric(unlist(t4b_u$unit_dept))
var<- as.numeric(unlist(attr(t4b_u$unit_dept,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t4b_s<- ranef(model4.adj,condVar = TRUE)
est<-as.numeric(unlist(t4b_s$season_hosp_date))
var<- as.numeric(unlist(attr(t4b_s$season_hosp_date,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))



## 3.4. sensitivity analysis #2
# creating a 1 and 2 day lag in long dataset for unit and time-varying covars
# unit, any antibiotics, any medication, any procedures
datlong_vars <- 
  datlong_vars %>%
  group_by(PatientPersonID) %>%
  mutate(lag.unit_dept1 = dplyr::lag(unit_dept, n = 1, default = NA)) %>%
  mutate(lag.unit_dept2 = dplyr::lag(unit_dept, n = 2, default = NA)) %>%
  mutate(lag.hosp_day01 = dplyr::lag(hosp_day0, n = 1, default = NA)) %>%
  mutate(lag.hosp_day02 = dplyr::lag(hosp_day0, n = 2, default = NA)) %>%
  mutate(lag.any_abx1 = dplyr::lag(AnyDailyAntibiotics, n = 1, default = NA)) %>%
  mutate(lag.any_abx2 = dplyr::lag(AnyDailyAntibiotics, n = 2, default = NA)) %>%
  mutate(lag.any_med1 = dplyr::lag(AnyDailyMedications, n = 1, default = NA)) %>%
  mutate(lag.any_med2 = dplyr::lag(AnyDailyMedications, n = 2, default = NA)) %>%
  mutate(lag.any_proc1 = dplyr::lag(AnyDailyProcedures, n = 1, default = NA)) %>%
  mutate(lag.any_proc2 = dplyr::lag(AnyDailyProcedures, n = 2, default = NA)) 
datlong_lag1 <- datlong_vars[!is.na(datlong_vars$lag.unit_dept1),]
datlong_lag2 <- datlong_vars[!is.na(datlong_vars$lag.unit_dept2),]

## 1 day lag
model5.empty <- glmer(CdiffPositive ~ (1|PatientPersonID) + (1|lag.unit_dept1), 
                      data=datlong_lag1, family=binomial(link="logit"),
                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model5.empty)

model5.adj <- glmer(CdiffPositive ~ lag.any_abx1 + lag.any_med1 + lag.any_proc1 
                    + AgeInYears + Insurance_cat + AHRQComorbidityCnt + lag.hosp_day01 + cum_units
                    + (1|PatientPersonID) + (1|lag.unit_dept1), 
                    data=datlong_lag1, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model5.adj)
tidy(model5.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
## 2 day lag
model6.empty <- glmer(CdiffPositive ~ (1|PatientPersonID) + (1|lag.unit_dept2), 
                      data=datlong_lag2, family=binomial(link="logit"),
                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model6.empty)

model6.adj <- glmer(CdiffPositive ~ lag.any_abx2 + lag.any_med2 + lag.any_proc2 
                    + AgeInYears + Insurance_cat + AHRQComorbidityCnt + lag.hosp_day02 + cum_units
                    
                    + (1|PatientPersonID) + (1|lag.unit_dept2), 
                    data=datlong_lag2, family=binomial(link="logit"),
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 20000)));summary(model6.adj)
tidy(model6.adj,conf.int=TRUE,exponentiate=TRUE,effects="fixed")


model.list<-list(model5.empty, model5.adj, model6.empty, model6.adj)
## Median odds ratio
model.list_output<-map_dfr(model.list, function(model){
  coefs<-tidy(model) %>% filter(effect=="ran_pars")%>%
    mutate(tau2=estimate^2,
           MOR=exp(0.95*sqrt(tau2))) %>%
    select(group, tau2, MOR)
  coefs
})
model.list_output

# bootstrapped MOR confidence intervals
t5a_p<- ranef(model5.empty,condVar = TRUE)
est<-as.numeric(unlist(t5a_p$PatientPersonID))
var<- as.numeric(unlist(attr(t5a_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t5a_u<- ranef(model5.empty,condVar = TRUE)
est<-as.numeric(unlist(t5a_u$lag.unit_dept1))
var<- as.numeric(unlist(attr(t5a_u$lag.unit_dept1,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t5b_p<- ranef(model5.adj,condVar = TRUE)
est<-as.numeric(unlist(t5b_p$PatientPersonID))
var<- as.numeric(unlist(attr(t5b_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t5b_u<- ranef(model5.adj,condVar = TRUE)
est<-as.numeric(unlist(t5b_u$lag.unit_dept1))
var<- as.numeric(unlist(attr(t5b_u$lag.unit_dept1,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t6a_p<- ranef(model6.empty,condVar = TRUE)
est<-as.numeric(unlist(t6a_p$PatientPersonID))
var<- as.numeric(unlist(attr(t6a_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t6a_u<- ranef(model6.empty,condVar = TRUE)
est<-as.numeric(unlist(t6a_u$lag.unit_dept2))
var<- as.numeric(unlist(attr(t6a_u$lag.unit_dept2,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t6b_p<- ranef(model6.adj,condVar = TRUE)
est<-as.numeric(unlist(t6b_p$PatientPersonID))
var<- as.numeric(unlist(attr(t6b_p$PatientPersonID,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))
#
t6b_u<- ranef(model6.adj,condVar = TRUE)
est<-as.numeric(unlist(t6b_u$lag.unit_dept2))
var<- as.numeric(unlist(attr(t6b_u$lag.unit_dept2,"postVar")))
mor_boot<- c()
for(i in 1:1000){
  drw<- rnorm(n = length(est),mean = est,sd = sqrt(var))
  s<- combn(drw,2)
  mor_boot<- c(mor_boot, median(exp(abs(s[1,]- s[2,]))))
}
quantile(mor_boot, c(.025,.5,.975))

