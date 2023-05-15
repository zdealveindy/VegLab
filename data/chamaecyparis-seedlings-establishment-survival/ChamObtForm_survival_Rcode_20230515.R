library(readxl)
library(survival)
library(dplyr)
library(survminer)
library(StepReg)
library(rms)
library(coxme)
library(janitor)
library(ggplot2)
library(broom)
library(CoxR2)
library(relimp)
library(etm)
library(cmprskcoxmsm)
setwd("/Users/weishuo/Desktop")

#### Data preparation ####
# Import data from survival survey
sur <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_survival_seedling_level_20230311.txt')[,c(-20)]
# Import data for plot- and individual- level variables
env <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_survival_plot_level_20230311.txt')
# Import data for region- level variables
cli <- read_delim ('https://raw.githubusercontent.com/zdealveindy/VegLab/main/data/chamaecyparis-seedlings-establishment-survival/ChamObtForm_survival_regional_level_20230311.txt')
# Species cover conversion
data <- left_join(sur, env, by = c("Plot_ID", "Subtype"))

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
data$herb_cv <- cover_change(data = data$herb_cv)
data$Sd_cv <- cover_change(data = data$Sd_cv)
data$Sp1_cv <- cover_change(data = data$Sp1_cv)
data$Sp2_cv <- cover_change(data = data$Sp2_cv)
data$herb_cv[is.na(data$herb_cv)] <- 0
data$Sd_cv[is.na(data$Sd_cv)] <- 0
data$Sp1_cv[is.na(data$Sp1_cv)] <- 0
data$Sp2_cv[is.na(data$Sp2_cv)] <- 0
data$Sp_cv <- data$Sp1_cv + data$Sp2_cv # merge cover of sp1 with sp2

# calculate the number of Chamaecyparis obtusa var. formosana seedling per plot
CHA_num <- data %>% filter(Species == "Chamaecyparis_obtusa_formosana") %>% group_by(Plot_ID, Subtype) %>% 
  summarise(CHA_num = n())
data <- left_join(data,CHA_num, by = c("Plot_ID", "Subtype"))
data$CHA_num[is.na(data$CHA_num)] <- 0
data$herbivory <- as.factor(data$herbivory)
data$cover <- as.factor(data$cover)
data$M_cv <- data$cv*data$moss_cv/100
data$L_cv <- data$cv*data$litter_cv/100
data$CHA_num <- data$CHA_num/data$cv
data$Subtype <- factor(data$Subtype, levels = c("soil", "CWD_nonD","CWD_D","CWD_halfD","mat"))
data$M_L <- factor(data$M_L, levels = c("L","S","M"))
data$status <- factor(data$status, levels = c("L","D","H","M"))
data$time2 <- data$time_no7-1
data <- data %>% filter(Species == "Chamaecyparis_obtusa_formosana")
data <- data %>% filter(time2 != 0)
data.nomissing <- data #%>% filter(status4 !="M") 
data.cli <- merge(data.nomissing, cli,  by = c("Plot_ID"))

# Table for showing the fate of the seedlings in different regions
data %>% tabyl(status, Region.x, show_na= T) %>% 
  adorn_totals("both") %>% 
  adorn_percentages(denominator = "all") %>%  # convert to proportions
  adorn_pct_formatting() %>%                  # convert to percents
  adorn_ns(position = "front") %>%            # display as: "count (percent)"
  adorn_title(                                # adjust titles
    row_name = "Substrate types",
    placement = "combined") %>% 
  flextable::flextable()%>% flextable::autofit()

# Fig. 4 CIF from competing risk
CIF.fit <- survfit(Surv(time2, status, type="mstate")~1, data=data.nomissing)
png("CIF.png",  res = 600,width = 8, height = 6,units = "in")
plot(CIF.fit, ylim=c(0,0.15), xlim=c(0,4), lty=1, lwd=2,
     col=c("#d8b365" , "#A3B18A", "grey"), xlab="Time (month)", ylab="CIF Estimates",xaxt = "n")
axis(1, at = 0:4, labels = c("Sep","Dec", "Mar", "Jun", "Sep"))
legend("topleft", lwd=2, bty="n",
       col=c("#d8b365" , "#A3B18A", "grey"), legend=c("Environment", "Herbivory", "Missing"))
dev.off()

cif.abortion <- etmCIF(Surv(time2, status != "L")~1, data=data.nomissing, etype = status)
summary(cif.abortion)
plot(cif.abortion, which.cif = c("D", "H", "M"), ylim = c(0, 0.2), xlim = c(0, 5), ci.type = "pointwise", col = c("#d8b365" , "#A3B18A", "grey"), lwd = 1.5, lty = 1, cex = 1.3)

D_mean <- list(0, 0.01150, 0.01150, 0.04948,0.11737, 0.11737)
D_upper <- list(0, 0.02128, 0.02128, 0.06614,0.14068, 0.14068)
D_lower <- list(0, 0.00620, 0.00620, 0.03693,0.09770, 0.09770)
H_mean <- list(0, 0.00920, 0.037974, 0.040276, 0.057537, 0.057537)
H_upper <- list(0, 0.01832, 0.05300, 0.05564, 0.07521, 0.07521)
H_lower <- list(0, 0.00461, 0.02714, 0.02908, 0.04391, 0.04391)
M_mean <- list(0, 0.034522, 0.06098, 0.066743, 0.09321, 0.09321)
timstp <- list(0, 1, 2, 3, 4, 5)
timstp_H <- list(0.005, 1.005, 2.005, 3.005, 4.005, 5)
timstp_D <- list(-0.005, 0.995, 1.995, 2.995, 3.995, 5)

png("CIF.png",  res = 600,width = 8, height = 6,units = "in")
plot(as.integer(timstp), D_mean, type = 'n', ylim = c(0, 0.15), xlim = c(1,5), col = "#d8b365", lwd = 1.5, xlab = "Season (month)", ylab = "CIF estimates", xaxt = "n")
axis(1, at = 1:5, labels = c("Sep","Dec", "Mar", "Jun", "Sep"))
axis(1, at = 1.5:4.5, labels = c("Autumn","Winter","Spring", "Summer"), tick = F, padj = 0.9)
polygon(c(1, 1, 3 ,3),c(0.02128, 0.00620, 0.00620, 0.02128), col = rgb(216, 179, 101, max = 255, alpha = 75), border=NA)
polygon(c(3 ,3, 4 ,4),c(0.06614, 0.03693, 0.03693, 0.06614), col = rgb(216, 179, 101, max = 255, alpha = 75), border=NA)
polygon(c(4 ,4, 5, 5),c(0.14068, 0.09770, 0.09770, 0.14068), col = rgb(216, 179, 101, max = 255, alpha = 75), border=NA)

polygon(c(1, 1, 2 ,2),c(0.01832, 0.00461, 0.00461, 0.01832), col = rgb(163, 177, 138, max = 255, alpha = 75), border=NA)
polygon(c(2, 2, 3 ,3),c(0.05300, 0.02714, 0.02714, 0.05300), col = rgb(163, 177, 138, max = 255, alpha = 75), border=NA)
polygon(c(3 ,3, 4, 4),c(0.05564, 0.02908, 0.02908, 0.05564), col = rgb(163, 177, 138, max = 255, alpha = 75), border=NA)
polygon(c(4 ,4, 5, 5),c(0.07521, 0.04391, 0.04391, 0.07521), col = rgb(163, 177, 138, max = 255, alpha = 75), border=NA)
lines(as.integer(timstp), M_mean, type = 's', col = rgb(190, 190, 190, max = 255, alpha = 100), lwd = 2)
lines((timstp_H), H_mean, type = 's', col = "#A3B18A", lwd = 2)
lines((timstp_D), D_mean, type = 's', col = "#d8b365", lwd = 2)
lines(timstp_D, D_upper, type = 's', col = "#d8b365", lty = 3, lwd = 1.5)
lines(timstp_D, D_lower, type = 's', col = "#d8b365", lty = 3, lwd = 1.5)
lines(timstp_H, H_upper, type = 's', col = "#A3B18A", lty = 3, lwd = 1.5)
lines(timstp_H, H_lower, type = 's', col = "#A3B18A", lty = 3, lwd = 1.5)

legend("topleft", lwd=1.5, bty="n",
       col=c("#d8b365" , "#A3B18A", "grey"), legend=c("Environment", "Herbivory", "Missing"))
dev.off()

#### Data analysis ####
#### ANALYSIS: for environmental-caused mortality (remove death from missing or herbivory) ####
# construct the full model contain all the individual and plot-level variables
# status = most conservative
coxph.plot.m0 <- coxph(Surv(time2, status == "D") ~  scale(Dsurf) + scale(ht) + herbivory + cover + M_L + scale(Canopy_cv) + scale(CHA_num) + scale(moss_th) + scale(M_cv) + scale(litter_th) + scale(L_cv) + scale(herb_cv) + scale(Sp_cv) + scale(Sd_cv) + Subtype + strata(Region.x), data = data.nomissing, ties = "efron", outer.max = 1000, iter.max = 10000, cluster = Plot_ID)

# stepwise selection (retain strata(Region.x) in the selection process)
stats::step(coxph.plot.m0, scope = list(upper = coxph.plot.m0,lower = Surv(time2, status == "D") ~ strata(Region.x)), direction = "both")

coxph.plot.1 <- coxph(Surv(time2, status == "D") ~ cover + scale(moss_th) + scale(M_cv) + scale(herb_cv) + strata(Region.x) , data = data.nomissing, ties = "efron", outer.max=1000, iter.max=10000, cluster = Plot_ID)
cox.zph(coxph.plot.1, transform = "rank") # PH-assumption testing -> M_cv and herb_cv violate PH-assumption

# Final model
coxph.plot.t1 <- coxph(Surv(time2, status == "D")~ cover + scale(moss_th) +  scale(M_cv) + scale(herb_cv) + tt(scale(moss_th))+ tt(scale(herb_cv)) + strata(Region.x), data = data.nomissing, ties = "efron", outer.max = 1000, iter.max = 10000, cluster = Plot_ID, tt = function(x,t,...) x*t)

# Examine whether time-varying variable should be differnt structures
coxph.plot.logt <-  coxph(Surv(time2, status == "D")~ cover + scale(moss_th) +  scale(M_cv) + scale(herb_cv) + tt(scale(moss_th))+ tt(scale(herb_cv)) + strata(Region.x), data = data.nomissing, ties = "efron",outer.max=1000, iter.max=10000, cluster = Plot_ID,
                          tt = function(x,t,...) x*log1p(t))
coxph.plot.logtt <-  coxph(Surv(time2, status == "D")~ cover + scale(moss_th) +  scale(M_cv) + scale(herb_cv) + tt(scale(moss_th))+ tt(scale(herb_cv))+ strata(Region.x), data = data.nomissing, ties = "efron",outer.max=1000, iter.max=10000, cluster = Plot_ID,
                           tt = function(x,t,...) cbind(x*log1p(t), x*t))

AIC(coxph.plot.t1, coxph.plot.logt, coxph.plot.logtt) # time-varing model (x*t) has the lowest AIC among other time-varing structures

coxr2(coxph.plot.t1) # R2 = 0.392
AIC(coxph.plot.t1, coxph.plot.1) # model with time-varying variables has lower AIC
car::Anova(coxph.plot.t1)
rms::vif(coxph.plot.t1)
summary(coxph.plot.t1)
# Examine whether there is strong pattern in residuals
ggcoxdiagnostics(coxph.plot.t1, type = "dfbeta",
                 linear.predictions = F, ggtheme = theme_bw())

broom::glance(coxph.plot.t1) # glimpse on model output

# Summary table
coxph.env.table <- coxph.plot.t1 %>% broom::tidy(exponentiate = F, conf.int = TRUE) %>%       
  mutate(across(where(is.numeric), round, digits = 3)) %>% 
  flextable::flextable()%>% flextable::autofit()
flextable::save_as_docx("coxph.env.table" = coxph.env.table, path = "coxph.env.table.docx")

# Fig. S2 Forest plot for the model of environmental-filtering survival analysis
a <- coxph.plot.t1 %>% tidy(exponentiate = F, conf.int = TRUE) 
a$term <- factor(a$term, levels = rev(a$term))
ggplot(a, aes(x = estimate, y = term)) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = conf.high, xmin = conf.low), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")+
  scale_x_continuous(breaks = (seq(-4, 4, 1)), labels = seq(-4, 4, 1),
                     limits = (c(-4,4)))+
  scale_y_discrete(labels = c("tt(Herb cover)", "tt(Bryophyte thickness)", "Herb cover", "Bryophyte cover", "Bryophyte thickness", "Litter covered"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  ylab("") +
  xlab("Hazard ratio (log scale)") +
  annotate(geom = "text", y =1.2, x = 1.5, 
           label = "LR~test~p<0.001\nAdjusted~pseudo~R^{2}==0.39", parse = T, size = 4, hjust = 0)+
  annotate(geom = "text", y =0.8, x = 1, 
           label = "Adjusted~pseudo~R^{2}==0.39", parse = T, size = 4, hjust = 0)
ggsave("survival_env_model.png", dpi = 600,width = 8, height = 6)

#### ANALYSIS: for herbivory-environmental-caused mortality ####
coxph.herb.m0 <- coxph(Surv(time2, status == "H")~  scale(Dsurf) + scale(ht) + M_L + scale(Canopy_cv) + scale(CHA_num) + scale(moss_th) + scale(M_cv) + scale(litter_th) + scale(L_cv) + scale(herb_cv) + scale(Sp_cv) + scale(Sd_cv) + Subtype + strata(Region.x), data = data.nomissing, ties = "efron", outer.max = 1000, iter.max = 10000, cluster = Plot_ID)

# stepwise selection (retain strata(Region.x) in the selection process)
stats::step(coxph.herb.m0, scope = list(upper = coxph.herb.m0,lower = Surv(time2, status == "H") ~ strata(Region.x), direction = "both"))

coxph.herb.1 <- coxph(Surv(time2, status == "H")~ scale(Dsurf) + scale(ht) +  scale(Canopy_cv) + scale(moss_th) + scale(L_cv) + scale(herb_cv) + strata(Region.x) , data = data.nomissing, ties = "efron",outer.max=1000, iter.max=10000, cluster = Plot_ID)
cox.zph(coxph.herb.1, transform = "rank") # PH-assumption testing ->  L_cv violate PH-assumption

coxph.herb.t1 <- coxph(Surv(time2, status == "H")~ scale(Dsurf) + scale(ht) +  scale(Canopy_cv) + scale(moss_th) + scale(L_cv) + scale(herb_cv) + tt(scale(L_cv))+ strata(Region.x) , data = data.nomissing, ties = "efron",outer.max=1000, iter.max=10000, cluster = Plot_ID,tt = function(x,t,...) x*t)

coxr2(coxph.herb.t1) # R2 = 0.613
AIC(coxph.herb.1, coxph.herb.t1)
summary(coxph.herb.t1)
broom::glance(coxph.herb.t1) # glimpse on model output

# Summary table
coxph.herb.t.table <- coxph.herb.t1 %>% broom::tidy(exponentiate = F, conf.int = TRUE) %>%       
  mutate(across(where(is.numeric), round, digits = 3)) %>% 
  flextable::flextable()%>% flextable::autofit()
flextable::save_as_docx("coxph.herb.status1.table" = coxph.herb.t.table, path = "coxph.herb.status1.table.docx")

# Fig. S3 Forest plot for the model of herbivory survival analysis
a <- coxph.herb.t1 %>% tidy(exponentiate = F, conf.int = TRUE) 
a$term <- factor(a$term, levels = rev(a$term))
ggplot(a, aes(x = estimate, y = term)) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = conf.high, xmin = conf.low), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange")+
  scale_x_continuous(breaks = (seq(-2, 2, 1)), labels = seq(-2, 2, 1),
                     limits = (c(-2,2)))+
  scale_y_discrete(labels = c("tt(Litter cover)", "Herb cover", "Litter cover", "Bryophyte thickness", "Canopy cover", "Height", "Dist. to soil"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  ylab("") +
  xlab("Hazard ratio (log scale)") +
  annotate(geom = "text", y =1.2, x = 0.7, 
           label = "LR~test~p<0.001\nAdjusted~pseudo~R^{2}==0.61", parse = T, size = 4, hjust = 0)+
  annotate(geom = "text", y =0.8, x = 0.4, 
           label = "Adjusted~pseudo~R^{2}==0.61", parse = T, size = 4, hjust = 0)
ggsave("survival_herb_model.png", dpi = 600,width = 8, height = 6)

#### ANALYSIS: for regional climatic variables ####
data.cli <- merge(data.nomissing, cli,  by = c("Plot_ID"))
data.cli <- data.cli %>% mutate(
  Fog_spring = (Fog_3+Fog_4+Fog_5)/3,
  Fog_summer = (Fog_6+Fog_7+Fog_8)/3,
  Fog_autumn = (Fog_9+Fog_10+Fog_11)/3,
  Fog_winter = (Fog_12+Fog_1+Fog_2)/3)
cli <- cli %>% mutate(
  Fog_spring = (Fog_3+Fog_4+Fog_5)/3,
  Fog_summer = (Fog_6+Fog_7+Fog_8)/3,
  Fog_autumn = (Fog_9+Fog_10+Fog_11)/3,
  Fog_winter = (Fog_12+Fog_1+Fog_2)/3)

# Non time-dependent analysis for cause-specific Cox PH models from survival monitoring on regional climatic variables.
coxme.cli.m0 <- coxph(Surv(time2, status == "D") ~ scale(Fog_year) + scale(Total_prec) + scale(Annual_tme), data = data.cli, ties = "efron")
coxme.cli.m1 <- coxme(Surv(time2, status == "D") ~ scale(Fog_year) + scale(Total_prec) + scale(Annual_tme) + (1|Region.x/Plot_ID), data = data.cli, ties = "efron")
anova(coxme.cli.m0, coxme.cli.m1)
summary(coxme.cli.m1)
MuMIn::r.squaredLR(coxme.cli.m1)
rms::vif(coxme.cli.m0)

coxme.herb.m0 <- coxph(Surv(time2, status == "H")~  scale(Fog_year)+ scale(Total_prec)+ scale(Annual_tme) , data = data.cli, ties = "efron")
coxme.herb.m1 <- coxme(Surv(time2, status == "H")~  scale(Fog_year)+ scale(Total_prec)+ scale(Annual_tme)  + (1|Region.x/Plot_ID), data = data.cli, ties = "efron")
anova(coxme.herb.m0, coxme.herb.m1)
summary(coxme.herb.m1)
MuMIn::r.squaredLR(coxme.herb.m1)
rms::vif(coxme.herb.m1)

# Time-dependent analysis for cause-specific Cox PH models from survival monitoring on regional climatic variables.
x <- cli %>% pivot_longer(cols = c("Fog_spring", "Fog_summer", "Fog_autumn", "Fog_winter"),
             names_to = "Season",
             values_to = "Fog") %>% 
  group_by(Region, Season) %>% summarise(meanFog = mean(Fog))
x$Season[which(x$Season== "Fog_spring")] <- "1"
x$Season[which(x$Season== "Fog_summer")] <- "2"
x$Season[which(x$Season== "Fog_autumn")] <- "3"
x$Season[which(x$Season== "Fog_winter")] <- "4"

y <- cli %>% pivot_longer(cols = c("Spr_prec", "Sum_prec", "Aut_prec", "Win_prec"),
                          names_to = "Season",
                          values_to = "Prec") %>% 
  group_by(Region, Season) %>% summarise(meanPrec = mean(Prec))
y$Season[which(y$Season== "Spr_prec")] <- "1"
y$Season[which(y$Season== "Sum_prec")] <- "2"
y$Season[which(y$Season== "Aut_prec")] <- "3"
y$Season[which(y$Season== "Win_prec")] <- "4"

z <- cli %>% pivot_longer(cols = c("Spr_tmean_", "Sum_tmean_", "Aut_tmean_", "Win_tmean_"),
                          names_to = "Season",
                          values_to = "Temp") %>% 
  group_by(Region, Season) %>% summarise(meanTemp = mean(Temp))
z$Season[which(z$Season== "Spr_tmean_")] <- "1"
z$Season[which(z$Season== "Sum_tmean_")] <- "2"
z$Season[which(z$Season== "Aut_tmean_")] <- "3"
z$Season[which(z$Season== "Win_tmean_")] <- "4"
x$Region[x$Region == "TMM"] <-"TMS"
y$Region[y$Region == "TMM"] <-"TMS"
z$Region[which(z$Region == "TMM")]<-"TMS"

data.cli$event <- data.cli$status == "D"
try <- survSplit(data.cli,
                         cut = 1:4,
                         end = "time2",
                         event = "event")
try$Season <- "3"
try$Season[try$time2 == 1] <- "3" 
try$Season[try$time2 == 2] <- "4" 
try$Season[try$time2 == 3] <- "1" 
try$Season[try$time2 == 4] <- "2"
try2 <- left_join(try, x, by = c("Region.x" ="Region", "Season"))
try2 <- left_join(try2, y, by = c("Region.x" ="Region", "Season"))
try2 <- left_join(try2, z, by = c("Region.x" ="Region", "Season"))

# Environmental-filtering
cli.0 <- coxme(Surv(tstart,time2,event=="1") ~scale(meanFog)+scale(meanPrec)+scale(meanTemp)+(1|Region.x/Plot_ID), data=try2, ties = "efron" )
summary(cli.0)
MuMIn::r.squaredLR(cli.0)
rms::vif(cli.0)
car::Anova(cli.0, type = "3")

# Herbivory
data.cli$event <- data.cli$status == "H"
try1 <- survSplit(data.cli,
                 cut = 1:4,
                 end = "time2",
                 event = "event")
try1$Season <- "3"
try1$Season[try1$time2 == 2] <- "4" 
try1$Season[try1$time2 == 3] <- "1" 
try1$Season[try1$time2 == 4] <- "2"
try1.2 <- left_join(try1, x, by = c("Region.x" = "Region",  "Season"))
try1.2 <- left_join(try1.2, y, by = c("Region.x" ="Region", "Season"))
try1.2 <- left_join(try1.2, z, by = c("Region.x" ="Region", "Season"))

cli.herb <- coxme(Surv(tstart,time2,event == "1") ~scale(meanFog)+scale(meanPrec)+scale(meanTemp)+(1|Region.x/Plot_ID), data=try1.2, ties = "efron" ) 
#anova(cli.herb.0,cli.herb)
summary(cli.herb)
MuMIn::r.squaredLR(cli.herb)
rms::vif(cli.herb)
car::Anova(cli.herb)
