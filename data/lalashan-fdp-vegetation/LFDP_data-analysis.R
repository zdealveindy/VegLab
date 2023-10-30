# Lalashan Forest Dynamics Plot - woody species analysis, Chen et al. ----

## Load libraries ----
library (readxl)
library (tidyverse)
library (interp)  # for interpolation of elevation values
library (plot3D)
#devtools::install_github("nsj3/riojaExtra")  # needs to be installed before twinspanR
#devtools::install_github('zdealveindy/twinspanR')
library (twinspanR)  # install from github
library (vegan)
#devtools::install_github('sdray/adespatial')
library (adespatial)

## Definition of functions ----
tb2df <- function (tb)
{
  df <- as.data.frame (tb)
  row.names(df) <- df[,1]
  df <- df[,-1]
  return (df)
}

## Load data ----
#setwd ('h:\\Dropbox_select\\CLANKY\\LFDP\\forest vegetation\\')
setwd ('c:\\Users\\zeleny\\Dropbox\\CLANKY\\LFDP\\forest vegetation\\')
ind <- read_excel ('LFDP_data_v1.0.xlsx', sheet = 'individual_data')
check <- read_excel ('LFDP_data_v1.0.xlsx', sheet = 'checklist')
pile <- read_excel ('LFDP_data_v1.0.xlsx', sheet = 'pile_elevations')
env_100 <- read_excel ('LFDP_data_v1.0.xlsx', sheet = 'env_100')
env_025 <- read_excel ('LFDP_data_v1.0.xlsx', sheet = 'env_25')
env <- left_join (env_100, env_025[,-5], by = 'subplot') %>% tb2df()

# Transform soil variables ----
env$C <- log10(env$C)
env$tN <- log10(env$tN)
env$CN_ratio <- log10(env$CN_ratio)
env$eN <- log10(env$eN)
env$P <- log10(env$P)
env$K <- log10(env$K)
env$Ca <- log10(env$Ca)
env$Mg <- log10(env$Mg)
env$Cu <- log10(env$Cu*100)
env$Zn <- log10(env$Zn)


# aggregate branches
ind_agg <- ind %>% mutate (BA = pi*(DBH/2)^2) %>%
  group_by (tag_no) %>%
  summarise (subplot = unique (subplot),
             x = unique (x),
             y = unique (y),
             qx = (x + as.numeric (str_sub (subplot, start = 2, end = 2))*100)/10,
             qy = (y + as.numeric (str_sub (subplot, start = 4, end = 4))*100)/10,
             BA = sum (BA),
             DBH = sqrt (BA/pi)*2,
             DBH_class = ifelse (DBH <= 5, 1, ifelse (DBH > 5 & DBH <= 30, 2, 3)),
             no_branches = n()-1,
             species = unique (species)) %>%
  left_join (check[, c('species', 'family', 'leaf_type')], by = c('species' = 'species'))

env$conif_BA <- ind_agg %>% group_by (subplot) %>% summarise (BA = sum (BA[leaf_type == 'Conifer'])) %>% select (BA) %>% sqrt () %>% .$BA
env$decid_BA <- ind_agg %>% group_by (subplot) %>% summarise (BA = sum (BA[leaf_type == 'Deciduous'])) %>% select (BA) %>% sqrt () %>% .$BA 
env$evergr_BA <- ind_agg %>% group_by (subplot) %>% summarise (BA = sum (BA[leaf_type == 'Evergreen'])) %>% select (BA) %>% sqrt () %>% .$BA

# interpolate the piles' elevation for smooth drawing of contours
pile_sub <- pile[1:121,] 
IL <- interp (pile_sub$x, pile_sub$y, pile_sub$elevation, nx = 110, ny = 110, method="akima", kernel="gaussian",solver="QR", input = 'points', output = 'grid')
IL_df <- IL
eg <- expand.grid (list (x = IL$x, y = IL$y))
IL_df$x <- eg$x
IL_df$y <- eg$y
IL_df$z <- as.vector (IL$z)
IL_df <- as.data.frame (IL_df)

# Preparing IVI data ----
com_long <- ind_agg %>% group_by (subplot, species) %>%
  summarise (BA = sum (BA), dens = n()) %>%
  group_by (subplot) %>%
  mutate (IVI = ((BA/sum (BA)*100) + (dens/sum (dens)*100))/2)

com <- com_long %>% select (-BA, -dens) %>% 
  pivot_wider (names_from = species, values_from = IVI, values_fill = 0) %>%
  tb2df ()

com_BA <- com_long %>% select (-dens, -IVI) %>%
  pivot_wider (names_from = species, values_from = BA, values_fill = 0) %>%
  tb2df ()

# Vegetation classification (modified TWINSPAN) ----
tw <- twinspan (com, modif = T, clusters = 3)
veg_types <- cut (tw)
table (veg_types)

# tw_BA <- twinspan (log1p (com_BA), modif = T, clusters = 3, cut.levels = c(0, 2, 4, 6, 8, 10))
# veg_types_BA <- cut (tw_BA)
# table (veg_types_BA)

# Detrended correspondence analysis ----
DCA <- decorana (log10(com+1))
env_topo <- env[, c('elevation', 'convexity', 'slope', 'northeasterness', 'heatload_1', 'heatload_2', 'windwardness', 'soil_depth', 'rockiness', 'pH')]
ef.topo <- envfit (DCA, env_topo, permutations = how (within = Within (type = 'grid', ncol = 10, nrow = 10, mirror = TRUE), complete = TRUE, nperm = 999))
ef.soil <- envfit (DCA, env[, c(25:43)], na.rm = TRUE, permutations = how (within = Within (type = 'grid', ncol = 5, nrow = 5, mirror = TRUE), complete = TRUE))

# Analysis of soil decomposition data against env. vars ----
env_sub <- env[!is.na(env$s_buried),]

adespatial::forward.sel (Y = env_sub$k_combust, X = env_sub[,c(-1, -40:-43)])  # Fe and evergr_BA

adespatial::forward.sel (Y = env_sub$s_combust, X = env_sub[,c(-1, -40:-43)])  # northeasterness and conifer_BA

lm_k <- lm (scale (k_combust) ~ scale (Fe) + scale (evergr_BA), data = env_sub)
summary (lm_k)

lm_s <- lm (scale (s_combust) ~ scale (northeasterness) + scale (conif_BA), data = env_sub)
summary (lm_s)



com_sub <- com[!is.na (env$s_combust),]
com_sub <- com_sub[,colSums (com_sub)>0]

# Soil decomposition as expl. var for species composition ----
RDA_raw <- rda (com_sub ~ s_combust, data = env_sub)
summary (RDA_raw)
RsquareAdj (RDA_raw) # 5.75%

com_BA_sub <- com_BA[!is.na (env$s_combust),]
com_BA_sub <- com_BA_sub[,colSums (com_BA_sub)>0]
RDA_BA_raw <- rda (com_BA_sub ~ s_combust, data = env_sub)
summary (RDA_BA_raw)
RsquareAdj (RDA_BA_raw) # 7.81%
anova (RDA_BA_raw)

DCA_sub <- decorana (sqrt (com_BA[!is.na (env$s_combust),]))
envfit (DCA_sub, env_sub)
# Soil properties are more related to the analysis done on 
# species data without transformation (IVI, BA), which means
# they are more related to the pattern formed by dominant species
