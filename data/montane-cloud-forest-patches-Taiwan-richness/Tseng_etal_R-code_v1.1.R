########################
#Packages and functions#
########################
#Packages
library(tidyr)
library(dplyr)
library(car)
library(ggplot2)
library(patchwork)


#Function for drawing correlation plot----
#Modified chart.Correlation() from R package "PerformanceAnalytics"
cor_plot <- function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                                      "spearman"), ...) 
{
  x = PerformanceAnalytics::checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 3
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex)
    text(0.8, 0.8, Signif, cex = 2, col = "black")
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "grey", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel, ...)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
             ...)
}




###############
#read csv file#
###############
land <- read.csv('Tseng_etal_data.csv')

##############################################################
#Sampling for the plots to prevent pseudo-replication problem#
##############################################################
#Extract geographical coordinates of the plot
cor <- land %>% filter(buffer == 1000) %>% select("releve_num", "Tmx", "Tmy")
rownames(cor) <- cor[, "releve_num"]
cor <- cor[, -1]

set.seed(500)
shuffled_cor <- cor[sample(1:nrow(cor)), ]
shu_dist <- as.matrix(dist(shuffled_cor, method = "euclidean", diag = T))
shu_dist[upper.tri(shu_dist, diag = T)] <- NA

i <- 0
repeat{
  i <- i + 1
  if (length(which(shu_dist[, i] < 2000)) > 0){
    shu_dist <- shu_dist[-which(shu_dist[, i] < 2000), -which(shu_dist[, i] < 2000)]
  }
  
  if (length(which(!is.na(shu_dist[, i]))) == 0){
    break
  }
  if (length(which(!is.na(shu_dist[, (i + 1)]))) == 0){
    break
  }
}
#107 plots are selected

###############
#Data Analysis#
###############
#Prepare data for the correlation plot with specialist richness and landscape variables
land_shu <- land %>% filter(buffer == 1000) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))
land_shu_075 <- land %>% filter(buffer == 750) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))
land_shu_05 <- land %>% filter(buffer == 500) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))

dat_transform <- function(dat) {
  row.names(dat) <- dat[, 1]
  dat <- dat[, -1]
  dat <- scale(dat)
  dat <- dat[, c("specialist_richness", "SMCF_area_p", "con_len", "tri_mean_p", "veg_di", "near_dist", "edge_len")]
  colnames(dat) <- c("Species richness of specialist", "Area of SMCFs", "Connectivity", "Topographical heterogeneity", "Vegetation diversity", "Nearest distance to the edge", "Edge length")
  return(dat)
}

dat_shu <- dat_transform(land_shu)
dat_shu_075 <- dat_transform(land_shu_075)
dat_shu_05 <- dat_transform(land_shu_05)

#correlation plots
cor_plot(dat_shu, histogram = TRUE, method = "kendall")
cor_plot(dat_shu_075, histogram = TRUE, method = "kendall")
cor_plot(dat_shu_05, histogram = TRUE, method = "kendall")

#multiple linear model and VIF for the SMCF patches in each buffer zone
for(mod in c('land_shu', 'land_shu_075', 'land_shu_05')){
  fit.mod <- lm(specialist_richness ~ SMCF_area_p + veg_di + tri_mean_p + near_dist + edge_len , data = as.data.frame(get(mod)))
  print(summary(fit.mod))
  print(car::vif(fit.mod))
}

# Fig. 3 Partial added-plots of the multiple linear models with the species richness of specialists in SMCFs as the response variable, and the landscape variables as the independent variables within 1 km radius of the buffer zone. 
dat <- as.data.frame(land_shu)
model <- lm(specialist_richness ~ SMCF_area_p + veg_di + tri_mean_p + near_dist + edge_len, data = dat)
model_wo_area <- lm(specialist_richness ~ veg_di + tri_mean_p + near_dist + edge_len, data = dat)
model_wo_veg <- lm(specialist_richness ~ SMCF_area_p + tri_mean_p + near_dist + edge_len, data = dat)
model_wo_tri <- lm(specialist_richness ~ veg_di + SMCF_area_p + near_dist + edge_len, data = dat)
model_wo_nd <- lm(specialist_richness ~ veg_di + tri_mean_p + SMCF_area_p + edge_len, data = dat)
model_wo_edge <- lm(specialist_richness ~ veg_di + tri_mean_p + near_dist + SMCF_area_p, data = dat)

dat$resid_noarea <- residuals(model_wo_area)
dat$resid_noveg <- residuals(model_wo_veg)
dat$resid_notri <- residuals(model_wo_tri)
dat$resid_nond <- residuals(model_wo_nd)
dat$resid_noedge <- residuals(model_wo_edge)

model_area <- lm(SMCF_area_p ~ veg_di + tri_mean_p + near_dist + edge_len, data = dat)
model_veg <- lm(veg_di ~ SMCF_area_p + tri_mean_p + near_dist + edge_len, data = dat)
model_tri <- lm(tri_mean_p ~ SMCF_area_p + veg_di + near_dist + edge_len, data = dat)
model_nd <- lm(near_dist ~ SMCF_area_p + veg_di + tri_mean_p + edge_len, data = dat)
model_edge <- lm(edge_len ~ SMCF_area_p + veg_di + tri_mean_p + near_dist, data = dat)

dat$resid_area <- residuals(model_area)
dat$resid_veg <- residuals(model_veg)
dat$resid_tri <- residuals(model_tri)
dat$resid_nd <- residuals(model_nd)
dat$resid_edge <- residuals(model_edge)

txt_size <- 18
area <- ggplot(dat, aes(x = resid_area, y = resid_noarea)) +
  geom_point() + geom_smooth(method = 'lm', color = 'black') +
  xlab("The area of SMCFs | Others") +
  ylab("Species richness of specialists in SMCFs | Others)") +
  theme(text = element_text(size = txt_size)) 
veg <- ggplot(dat, aes(x = resid_area, y = resid_noveg)) +
  geom_point() +
  xlab("Vegetation diversity | Others") +
  ylab("Species richness of specialists in SMCFs | Others)")+
  theme(text = element_text(size = txt_size))
tri <- ggplot(dat, aes(x = resid_tri, y = resid_notri)) +
  geom_point() +
  xlab("Topographical heterogeneity | Others") +
  ylab("Species richness of specialists in SMCFs | Others)")+
  theme(text = element_text(size = txt_size))
nd <- ggplot(dat, aes(x = resid_nd, y = resid_nond)) +
  geom_point() +
  xlab("Nearest distance to the edge | Others") +
  ylab("Species richness of specialists in SMCFs | Others)")+
  theme(text = element_text(size = txt_size))
edge <- ggplot(dat, aes(x = resid_edge, y = resid_noedge)) +
  geom_point() +
  xlab("Edge length | Others") +
  ylab("Species richness of specialists in SMCFs | Others)")+
  theme(text = element_text(size = txt_size))
# to examine the results' robustness
#save plot
png("avplot.png",
    width=1200, height=900)
area + veg + tri + nd + edge
dev.off()


##############################################
#100-time re-sampling for results' robustness#
##############################################
set.seed(500)
ce <- c() #correlation coefficient
p <- c() #p-value
r2 <- c() #r-square
n <- c() 

for(r in 1:100)
{
  shuffled_cor <- cor[sample(1:nrow(cor)), ]
  shu_dist <- as.matrix(dist(shuffled_cor, method = "euclidean", diag = T))
  shu_dist[upper.tri(shu_dist, diag = T)] <- NA
  
  i <- 0
  repeat{
    i <- i + 1
    if (length(which(shu_dist[, i] < 2000)) > 0){
      shu_dist <- shu_dist[-which(shu_dist[, i] < 2000), -which(shu_dist[, i] < 2000)]
    }
    
    if (length(which(!is.na(shu_dist[, i]))) == 0){
      break
    }
    if (length(which(!is.na(shu_dist[, (i + 1)]))) == 0){
      break
    }
  }
  n <- c(n, nrow(shu_dist))
  
  land_shu <- land %>% filter(buffer == 1000) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))
  land_shu_075 <- land %>% filter(buffer == 750) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))
  land_shu_05 <- land %>% filter(buffer == 500) %>% filter(releve_num %in% as.numeric(rownames(shu_dist)))
  
  dat_transform_only <- function(dat) {
    row.names(dat) <- dat[, 1]
    dat <- dat[, -1]
    dat <- scale(dat)
    return(dat)
  }
  
  dat_shu <- dat_transform_only(land_shu)
  dat_shu_075 <- dat_transform_only(land_shu_075)
  dat_shu_05 <- dat_transform_only(land_shu_05)
  
  for(mod in c('dat_shu', 'dat_shu_075', 'dat_shu_05')){
    fit.mod <- lm(specialist_richness ~ SMCF_area_p + veg_di + tri_mean_p + near_dist + edge_len, data = as.data.frame(get(mod)))
    ce <- c(ce, summary(fit.mod)$coefficients[2:6, 1])
    p <- c(p, summary(fit.mod)$coefficients[2:6, 4])
    r2 <- c(r2, summary(fit.mod)$adj.r.squared)
  }
}

ce_mat <- matrix(ce, ncol = 15, nrow = 100, byrow = T)
p_mat <- matrix(p, ncol = 15, nrow = 100, byrow = T)
r2_mat <- matrix(r2, ncol = 3, nrow = 100, byrow = T)

#data frame for the mean and standard deviation of correlation coefficients as well as percentage of significant results for different landscape variables with 100-time re-sampling
cor_coe <- matrix(c(apply(ce_mat, 2, mean), apply(ce_mat, 2, sd), apply(p_mat, 2, function(x){length(which(x<0.05))/100})), nrow = 15, ncol = 3, byrow = F)
cor_coe <- as.data.frame(cor_coe)
cor_coe$variable <- rep(c("Area of SMCFs", "Vegetation diversity", "Topographical heterogeneity", "Nearest distance to the edge", "Edge length"), 3)
cor_coe$buffer <- rep(c(1000, 750, 500), each = 5)
colnames(cor_coe)[1:3] <- c("correlation coefficient", "sd", "percentage of significant results")
write.csv(cor_coe, 'cor_coe.csv')

#data frame for the mean and standard deviation of adjusted R-square for 100-time re-sampling models
r2_st_mat <- matrix(c(apply(r2_mat, 2, mean), apply(r2_mat, 2, sd)), nrow = 2, ncol = 3, byrow = T)
colnames(r2_st_mat) <- c(1000, 750 ,500)
rownames(r2_st_mat) <- c("mean", "sd")
write.csv(r2_st_mat, 'r2.csv')

