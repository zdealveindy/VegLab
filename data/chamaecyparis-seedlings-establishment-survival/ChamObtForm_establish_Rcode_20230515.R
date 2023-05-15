library(tidyverse)
library(dplyr)
library(lme4)
library(survival)
library(lmerTest)
library(DHARMa)
library(ggplot2)
library(geepack)
library(mclogit)
library(gnm)
library(jtools)
library(janitor)
library (readr) 
#### Data preparation ####

data <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_establish_plot_record_20220805.txt')
data <- data[-c(seq(from = 6, to = 72, by =6)),] # remove the extra one random plot from TMS, TSS and WDY
seedling <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_establish_seedling_record_20220805.txt')
GIS <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_establish_GIS_20220805.txt')

data$Sub_type <- factor(data$Sub_type, levels = c("Soil","CWD_D","CWD_halfD", "CWD_nonD","Mat","Living"))
data <- left_join(data, GIS, by = "Plot_ID")


# Species cover conversion
cover_change <- function(data = data){
  # transform Sd_cv_1 scale into percentage scale
  ifelse(is.na(data), 0, data)
  data [data == "0"] <- 0 # 0
  data [data == "3"] <- 37.5 # 3
  data [data == "1"] <- 3 # 1
  data [data == "r"] <- 1 # r 
  data [data == "+"] <- 2 # +
  data [data == "2a"] <- 10 # 2a
  data [data == "2b"] <- 20 # 2b
  data [data == "4"] <- 62.5 # 4
  data [data == "5"] <- 87.5 # 5
  data [data == "2m"] <- 5 # 2m
  data  <- as.numeric (data)
  ifelse(is.na(data), 0, data)
}
data$Herb_cv <- cover_change(data = data$Herb_cv)
data$Sd_cv <- cover_change(data = data$Sd_cv)
data$Sp_cv <- cover_change(data = data$Sp_cv)
data$`Dead Yushania cover` <- cover_change(data = data$`Dead Yushania cover`)
data$Herb_cv[is.na(data$Herb_cv)] <- 0
data$Sd_cv[is.na(data$Sd_cv)] <- 0
data$Sp_cv[is.na(data$Sp_cv)] <- 0
data$`Dead Yushania cover`[is.na(data$`Dead Yushania cover`)] <- 0
seedling %>% filter(Species == "Chamaecyparis_obtusa_formosana") %>% 
  group_by(Plot_ID) %>% 
  summarise(sd_num = n()) -> sd
sd$Sd_class <- 0

sd$Sd_class [sd$sd_num <= 10] <- 1 
sd$Sd_class [sd$sd_num > 10 &sd$sd_num <=15 ] <- 2 
sd$Sd_class [sd$sd_num > 15 ] <- 3


data <- left_join(data, sd, by = "Plot_ID")
data$sd_num[is.na(data$sd_num)] <- 0
data$Sd_class[is.na(data$Sd_class)] <- 0

data$Herb_ht[data$Herb_ht == "NA"] <- 0


#### Data analysis of conditional logistic regression ####
# full model (all variable without interaction and polynomial terms)
m1 <- clogit(formula = Case ~ Sub_type + scale(M_cv) + scale(M_th) + scale(L_th)+scale(L_cv)  + scale(Sd_cv) + scale(Sp_cv) + scale(Herb_cv) + strata(ID.x) + scale(Canopy_cv), data = data, method = "efron")

m0 <- clogit(formula = Case ~ Sub_type + strata(ID.x) , data = data , method = "efron")

AIC <- step(m0,
               scope=list(lower=m0, upper=m1), "both")

BIC <- step(m0,
            scope=list(lower=m0, upper=m1), "both", k = log(370))  

summary(AIC) #Sub_type + strata(ID.x) + scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th) + scale(L_cv)
summary(BIC) # Sub_type + strata(ID.x) + scale(L_th) + scale(M_cv) + scale(Canopy_cv)

car::Anova(AIC)
car::Anova(BIC)
anova(AIC,BIC) # AIC model has better performance than BIC model by LRT (p = 0.028)
rms::vif(AIC)
rms::vif(BIC)

# BIC model
m2.0 <- clogit(formula = Case ~ strata(ID.x) +
                 Sub_type + 
                 scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th) + scale(L_cv), 
                 data = data, method = "efron")

# Identify which variables contribute to better model of AIC (M_th or L_cv)
m2.1.1 <- clogit(formula = Case ~  strata(ID.x) +
                 Sub_type +
                 scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(L_cv) ,
                 data = data, method = "efron")
anova(m2.0, m2.1.1) # omit M_th is bad (LRT p < 0.001) -> retain L_th
m2.1.2 <- clogit(formula = Case ~  strata(ID.x) +
                   Sub_type +
                   scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th) , 
                 data = data , method = "efron")
anova(m2.0, m2.1.2) # omit L_cv do not influence the model performance (LRT p = 0.152) -> byebye L_cv

# Reduce the model from AIC to m2.1.2, and testing whether there should be a polynomial terms
# polynomial for M_th
m2.2.1 <- clogit(formula = Case ~  strata(ID.x) + 
                 Sub_type + scale(M_th) +
                 scale(L_th) + scale(M_cv) + scale(Canopy_cv)+ I(scale(M_th)^2) , data = data , method = "efron")
# polynomial for L_th
m2.2.2 <- clogit(formula = Case ~  strata(ID.x) + 
                   Sub_type + scale(M_th) +
                   scale(L_th) + scale(M_cv) + scale(Canopy_cv)+ I(scale(L_th)^2) , data = data , method = "efron")
# polynomial for M_cv
m2.2.3 <- clogit(formula = Case ~  strata(ID.x) + 
                   Sub_type + scale(M_th) +
                   scale(L_th) + scale(M_cv) + scale(Canopy_cv)+ I(scale(M_cv)^2) , data = data , method = "efron")
# polynomial for Canuopy_cv
m2.2.4 <- clogit(formula = Case ~  strata(ID.x) + 
                   Sub_type + scale(M_th) +
                   scale(L_th) + scale(M_cv) + scale(Canopy_cv)+ I(scale(Canopy_cv)^2) , data = data , method = "efron")
anova(m2.1.2, m2.2.1)
anova(m2.1.2, m2.2.2)
anova(m2.1.2, m2.2.3)
anova(m2.1.2, m2.2.4)
# all the polynomial model is not significant by LRT -> no need for polynomial terms

# test whether mother tree distance will affect the result
mother.effect <- coxph(formula = Surv(rep(1, 370L), Case) ~ strata(ID.x) +
                         Sub_type +
                         scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th) + factor(Mother),
                       data = data , method = "efron")
car::Anova(mother.effect)
summary(mother.effect)
anova(m2.1.2, mother.effect) # Mother does not significantly improve the model by LRT (p = 0.098)

#### FINAL MODEL
m2.1.2 <- clogit(formula = Case ~  strata(ID.x) +
                   Sub_type +
                   scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th) , data = data , method = "efron")
plot(resid(m2.1.2))
car::Anova(m2.1.2)
summary(m2.1.2)

1-((m2.1.2$loglik[2]-9)/m2.1.2$loglik[1]) # adjusted R2 for parameters = 77.3%
1-(m2.1.2$loglik[2]/m2.1.2$loglik[1]) # R2 84.8%


##### create table of regression summary ####
m.final.table <- m2.1.2 %>% tidy(exponentiate = T, conf.int = TRUE) %>%       
  mutate(across(where(is.numeric), round, digits = 3)) %>% 
  flextable::flextable()%>% flextable::autofit()
flextable::save_as_docx("m.final.table" = m.final.table, path = "m3table.docx")

data %>% tabyl(Sub_type, Case, show_na= T) %>% 
  adorn_totals("both") %>% 
  adorn_percentages(denominator = "all") %>%  # convert to proportions
  adorn_pct_formatting() %>%                  # convert to percents
  adorn_ns(position = "front") %>%            # display as: "count (percent)"
  adorn_title(                                # adjust titles
    row_name = "Substrate types",
    placement = "combined") %>% 
  flextable::flextable()%>% flextable::autofit() -> table
flextable::save_as_docx("Case/Control" = table, path = "file.docx")


#### understanding which sub type is different from each other  ####
#    -> CWD_D and CWD_halfD cannot be seperate from each other, but both are different from Soil and CWD_nonD. Mat and Living cannot be differentiated from any othe sub type.
# Soil is different from CWD_D and CWD_halfD
data$Sub_type <- factor(data$Sub_type, levels = c("CWD_D","CWD_halfD", "CWD_nonD","Mat","Living","Soil"))
# CWD_nonD and Soil is different from CWD_D, and CWD_halfD cannot be distinct from CWD_D
data$Sub_type <- factor(data$Sub_type, levels = c("CWD_halfD","CWD_D","Living", "CWD_nonD","Mat","Soil"))
data$Sub_type <- factor(data$Sub_type, levels = c("Living","CWD_D","CWD_halfD", "CWD_nonD","Mat","Soil"))
data$Sub_type <- factor(data$Sub_type, levels = c("CWD_nonD","CWD_D","Living", "CWD_halfD","Mat","Soil"))
data$Sub_type <- factor(data$Sub_type, levels = c("Mat","CWD_nonD","CWD_D","Living", "CWD_halfD","Soil"))

summary(clogit(formula = Case ~ Sub_type +
               scale(L_th) + scale(M_cv) + scale(Canopy_cv) + scale(M_th),
               data = data, method = "efron"))
