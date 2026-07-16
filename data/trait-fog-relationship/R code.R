# 1.	Library----
library(dplyr)
library(tidyr)
library(rio)
library(permute)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggpp)
# 2. Customized function----
# 2.2. Analysis----
# 2.2.1. Restricted-permutation-based Pearson's correlation test----
cor.test_HT <- function(x,
                        y,
                        n_perm = 999,
                        hierar = gl(n = 9, k = 3),
                        set_seed = 1211,
                        use = "everything",
                        na.action = "na.exclude"){
  # f.1. origin test----
  p_ori <- round(cor.test(x, y, na.action = na.action)$p.value, 3)
  cor_ori <- round(cor(x, y, use = use), 3)
  r2 <- round(cor(x, y, use = use) ^ 2, 3)
  perm_cor <- 1
  # f.2. permutation----
  HT_strcuture <- how(within = Within(type = "free"),
                      plots = Plots(strata = hierar, 
                                    type = "free"))
  
  # Create shuffle data frame
  set.seed(set_seed)
  shuffle_set <- shuffleSet(n = length(hierar), 
                            control = HT_strcuture,
                            nset = n_perm)
  for(i in 1:nrow(shuffle_set)){
    y_shuffle <- y[shuffle_set[i, ]]
    if(abs(cor(x, y_shuffle, use = use)) >= 
       abs(cor(x, y, use = use))){
      perm_cor <- perm_cor + 1
    }
  }
  p_shuff <- round(perm_cor / (n_perm + 1),
                   3)
  return(data.frame(n_perm = n_perm,
                    r = cor_ori,
                    r2 = r2,
                    p_ori = p_ori,
                    p_shuff = p_shuff))
}
# 2.2.2. Restricted-permutation-based simple linear regression----
slr_HT_fval <- function(formula1,
                        data1, # should include "y"
                        n_perm = 999,
                        hierar = gl(n = 9, k = 3),
                        set_seed = 1211){
  HT_strcuture <- 
    how(within =  Within(type = "free"),
        plots = Plots(strata = hierar, 
                      type = "free"))
  
  formula1 <- as.formula(formula1)
  # Create shuffle data frame
  set.seed(set_seed)
  shuffle_set <- shuffleSet(n = length(hierar), 
                            control = HT_strcuture,
                            nset = n_perm)
  # f.1. Calculated original p value----
  lm_ori <- lm(data1[, as.character(formula1[[2]])] ~ .,
               data = 
                 as.data.frame(data1[, gsub(pattern = " ", replacement = "",
                                            strsplit(as.character(formula1)[3], split = "[+]")[[1]])]))
  F_ori <- summary(lm_ori)$fstatistic[1]
  perm_F <- 1
  
  slope <- summary(lm_ori)$coefficients[2, 1]
  
  # f.2. Start permutation----
  for(i in 1:nrow(shuffle_set)){
    lm_shuff <- 
      lm(data1[shuffle_set[i, ], as.character(formula1[[2]])] ~ .,
         data = as.data.frame(data1[, gsub(pattern = " ", replacement = "",
                                           strsplit(as.character(formula1)[3], split = "[+]")[[1]])]))
    # extract p-value
    F_shuff <- summary(lm_shuff)$fstatistic[1]
    perm_F <- perm_F + ifelse(F_shuff > F_ori, 1, 0)
  }
  final_p <- perm_F / (n_perm + 1)
  
  info <- paste("Permutation time:", n_perm)
  return_summary <- list(summary(lm_ori), 
                         info, 
                         slope = slope,
                         final_p)
  
  names(return_summary) <- 
    c("lm_ori", "info", "slope", "final_p")
  return(return_summary)
  
}
# 2.2.3. Restricted-permutation-based multiple linear regression----
lmp_HT <- function(x1, x2, y,
                   n_perm = 999,
                   hierar = gl(n = 9, k = 3),
                   set_seed = 1211){
  # f.1. shuffle structure ----
  HT_strcuture <- how(within =  Within(type = "free"),
                      plots = Plots(strata = hierar, 
                                    type = "free"))
  set.seed(set_seed)
  shuffle_set <- shuffleSet(n = length(hierar), 
                            control = HT_strcuture,
                            nset = n_perm)
  # f.2. origin statistics----
  # f.2.1. whole model----
  F_ori <- summary(lm(y ~ x1 + x2))$fstatistic[1]
  slope_x1 <- 
    round(summary(lm(y ~ x1 + x2))$coefficients[2, 1],
          3)
  slope_x2 <- 
    round(summary(lm(y ~ x1 + x2))$coefficients[3, 1],
          3)
  vif <- round(car::vif(lm(y ~ x1 + x2)), 3)[1]
  F_p <- 1
  # f.2.2. beta1----
  resid_yx2 <- lm(y ~ x2)$residuals
  resid_x1x2 <- lm(x1 ~ x2)$residuals
  r2_x1_ori <- 
    sum((resid_yx2 * resid_x1x2) ^ 2) /
    (sum(resid_yx2 ^ 2) * sum(resid_x1x2 ^ 2))
  r2_x1_p <- 1
  # f.2.3. beta2----
  resid_yx1 <- lm(y ~ x1)$residuals
  resid_x2x1 <- lm(x2 ~ x1)$residuals
  r2_x2_ori <- 
    sum((resid_yx1 * resid_x2x1) ^ 2) /
    (sum(resid_yx1 ^ 2) * sum(resid_x2x1 ^ 2))
  r2_x2_p <- 1
  # f.3. shuffle----
  for(i in 1:nrow(shuffle_set)){
    # f.3.1. whole model----
    y_shuffle <- y[shuffle_set[i, ]]
    F_shuffle <- 
      summary(lm(y_shuffle ~ x1 + x2))$fstatistic[1]
    if(F_shuffle >= F_ori){
      F_p <- F_p + 1
    } # if
    # f.3.2. beta1----
    y_resid.shf <- predict(lm(y ~ x2)) + resid_yx2[shuffle_set[i, ]]
    resid_yx2_shuffle <- lm(y_resid.shf ~ x2)$residuals
    r2_x1_shuffle <- 
      sum((resid_yx2_shuffle * resid_x1x2) ^ 2) /
      (sum(resid_yx2_shuffle ^ 2) * sum(resid_x1x2 ^ 2))
    if(r2_x1_shuffle >= r2_x1_ori){
      r2_x1_p <- r2_x1_p + 1
    } # if
    # f.3.3. beta2----
    y_resid.shf <- predict(lm(y ~ x1)) + resid_yx1[shuffle_set[i, ]]
    resid_yx1_shuffle <- lm(y_resid.shf ~ x1)$residuals
    r2_x2_shuffle <- 
      sum((resid_yx1_shuffle * resid_x2x1) ^ 2) /
      (sum(resid_yx1_shuffle ^ 2) * sum(resid_x2x1 ^ 2))
    if(r2_x2_shuffle >= r2_x2_ori){
      r2_x2_p <- r2_x2_p + 1
    } # if
  } # for
  # f.4. return----
  beta1_p <- round(r2_x1_p / (n_perm + 1), 3)
  beta2_p <- round(r2_x2_p / (n_perm + 1), 3)
  model_p <- round(F_p / (n_perm + 1), 3)
  
  return(list(ori_summary = summary(lm(y ~ x1 + x2)),
              summary = 
                data.frame(
                  n_perm = n_perm,
                  slope_x1 = slope_x1,
                  p_x1 = beta1_p,
                  slope_x2 = slope_x2,
                  p_x2 = beta2_p,
                  model_p = model_p,
                  vif = vif
                )))
} # function
# 2.2.4. Restricted-permutation-based p_max test----
permute_test <- 
  function(trait_df,
           weight_df,
           env_vec, 
           type = NULL,
           n_perm = 999,
           hierar = gl(n = 9, k = 3),
           scale = TRUE,
           set_seed = 1211){
    set.seed(set_seed)
    # f.1. ensure there is no empty column----
    if(prod(!is.na(env_vec)) != 1){
      warning("There is NA in the environmental variable, but the analysis will keep going.")
    }
    if(!is.numeric(env_vec)){
      warning("There is character in the environmental variable.")
      stop()
    }
    
    if(length(which(colSums(is.na(trait_df)) == nrow(trait_df))) != 0){
      warning("Delete the column with all NA in trait_df")
      weight_df <- weight_df[, -which(colSums(is.na(trait_df)) == nrow(trait_df))]
      trait_df <- trait_df[, -which(colSums(is.na(trait_df)) == nrow(trait_df))]
      # to ensure weight_wide share same columns
    }
    if(length(which(colSums(is.na(weight_df)) == nrow(weight_df))) != 0){
      warning("Delete the column with all NA in weight_df")
      weight_df <- weight_df[, -which(colSums(is.na(weight_df)) == nrow(weight_df))]
      trait_df <- trait_df[, -which(colSums(is.na(weight_df)) == nrow(weight_df))]
      # to ensure weight_wide share same columns
    } # if loop
    
    weight_df <- sweep(weight_df, 
                       1, 
                       rowSums(weight_df, na.rm = TRUE),
                       FUN = "/")
    
    
    if(ncol(trait_df) != ncol(weight_df)){
      warning("After cleaning NA column, trait_df and weight_df don't have same column")
      stop()
    }
    
    # f.2. generate trait_df and weight_df----
    if(type == "spe"){
      trait_df <- trait_df
    } else if(type == "fixed"){
      trait_df <- 
        apply(trait_df,
              2,
              function(x){
                x[!is.na(x)] <- mean(x, na.rm = TRUE);
                x
              })
    } else if(type == "intra"){
      trait_df <- 
        apply(trait_df,
              2,
              function(x){
                x[!is.na(x)] <- 
                  scale(x[!is.na(x)], scale = FALSE);
                x
              })
    } else {
      warning("wrong type input")
    }
    
    # f.3. p_row----
    if(as.character(sum(weight_df, na.rm = TRUE)) != 
       as.character(nrow(weight_df))){
      warning("Sum of weight value is not equal to 1")
      print(paste(as.character(sum(weight_df, na.rm = TRUE)), "!=", as.character(nrow(weight_df))))
    }
    CWM <- rowSums(trait_df * weight_df, na.rm = TRUE)
    if(scale == TRUE){
      CWM_env_df <- 
        data.frame(CWM = scale(CWM),
                   env = scale(env_vec))
    } else {
      CWM_env_df <- 
        data.frame(CWM = CWM,
                   env = env_vec)
      
    }
    slope <- 
      round(summary(lm(CWM ~ env, 
                       data = CWM_env_df))$coefficients[2, 1], 
            3)
    Rsquare <- round(summary(lm(CWM ~ env, 
                                data = CWM_env_df))$r.square,
                     3)
    p_ori <- round(anova(lm(CWM ~ env, 
                            data = CWM_env_df))[1, 5],
                   3)
    # f.4. p_row ----
    # set.seed(set_seed)
    p_row <- slr_HT_fval(formula1 = CWM ~ env, 
                         data1 = CWM_env_df,
                         n_perm = n_perm,
                         hierar = hierar)$final_p
    if(type == "intra"){
      print("intra-specific CWM don't need to do p-max test")
      p_col <- "No"
      p_max <- p_row
    } else if(type == "fixed" | type == "spe"){
      # f.5. p_col ----
      F_ori <- summary(lm(CWM ~ env, 
                          data = CWM_env_df))$fstatistic[1]
      perm_F <- 1
      # set.seed(set_seed)
      for(i in 1:n_perm){
        shuff_trait_df <- 
          sweep(sweep(trait_df %>%
                        apply(2,
                              function(x){
                                x[!is.na(x)] <- 
                                  sample(x[!is.na(x)]);
                                x
                              }), 
                      2, 
                      colMeans(trait_df, 
                               na.rm = TRUE),
                      FUN = "-"),
                2,
                sample(colMeans(trait_df, 
                                na.rm = TRUE)),
                FUN = "+") %>%
          as.data.frame()
        # print("19 warnings is for the species which only appears once", )
        shuff_CWM <- 
          rowSums(shuff_trait_df * weight_df, na.rm = TRUE)
        shuff_CWM_env_df <- 
          data.frame(CWM = shuff_CWM,
                     env = env_vec)
        F_shuff <- 
          summary(lm(CWM ~ env, 
                     data = shuff_CWM_env_df))$fstatistic[1]
        perm_F <- perm_F + ifelse(F_shuff >= F_ori, 1, 0)
      } # shuffle for
      p_col <- perm_F/ (n_perm + 1)
      p_max <- max(p_row, p_col)
    }
    
    return(data.frame(permute = paste0("permute times: ", n_perm),
                      p_ori = p_ori,
                      p_col = p_col,
                      p_row = p_row,
                      type = type,
                      p_max = p_max,
                      Rsquare = Rsquare,
                      Coef = slope,
                      scale = scale)
    ) # return
  }

# 2.3. Pairwise correlation of monthly variables----
pairwise_cor_table <- function(major_df = NULL,
                               major_name = "major",
                               minor_df = NULL,
                               minor_name = "minor",
                               n_to_1 = TRUE){
  # colnames of major_df and minor_df should be consistent
  cor_table <- 
    data.frame(env1 = NULL,
               env2 = NULL,
               r = NULL,
               p = NULL)
  
  for(col in colnames(major_df)){
    if(n_to_1 == TRUE){
      cor.test_result <- cor.test_HT(major_df[, col],
                                     minor_df,
                                     n_perm = 999,
                                     set_seed = 1211)
      cor_table <- 
        rbind(cor_table,
              data.frame(env1 = paste0(major_name, "_", col),
                         env2 = minor_name,
                         r = round(cor.test_result$r, 3),
                         p = round(cor.test_result$p_shuff, 3)))
    } else {
      cor.test_result <- cor.test_HT(major_df[, col],
                                     minor_df[, col],
                                     n_perm = 999,
                                     set_seed = 1211)
      cor_table <- 
        rbind(cor_table,
              data.frame(env1 = paste0(major_name, "_", col),
                         env2 = paste0(minor_name, "_", col),
                         r = round(cor.test_result$r, 3),
                         p = round(cor.test_result$p_shuff, 3)))
    } # else 
  } # for
  return(cor_table)
}

# 2.4. Community-weighted mean calculation----
CWM_generator <- 
  function(trait_df,
           weight_df,
           type = NULL){
    # f.0. ensure there is no empty column----
    if(length(which(colSums(is.na(trait_df)) == nrow(trait_df))) != 0){
      warning("Delete the column with all NA in trait_df")
      weight_df <- weight_df[, -which(colSums(is.na(trait_df)) == nrow(trait_df))]
      trait_df <- trait_df[, -which(colSums(is.na(trait_df)) == nrow(trait_df))]
      # to ensure weight_wide share same columns
    }
    if(length(which(colSums(is.na(weight_df)) == nrow(weight_df))) != 0){
      warning("Delete the column with all NA in weight_df")
      weight_df <- weight_df[, -which(colSums(is.na(weight_df)) == nrow(weight_df))]
      trait_df <- trait_df[, -which(colSums(is.na(weight_df)) == nrow(weight_df))]
      # to ensure weight_wide share same columns
    } # if loop
    
    weight_df_rescale <- sweep(weight_df, 
                               1, 
                               rowSums(weight_df, na.rm = TRUE),
                               FUN = "/")
    if(ncol(trait_df) != ncol(weight_df_rescale)){
      warning("After cleaning NA column, trait_df and weight_df don't have same number of column")
      stop()
    }
    
    # f.1. generate trait_df and weight_df----
    if(type == "spe"){
      trait_df <- trait_df
    } else if(type == "fixed"){
      trait_df <- 
        apply(trait_df,
              2,
              function(x){
                x[!is.na(x)] <- mean(x, na.rm = TRUE);
                x
              })
    } else if(type == "intra"){
      trait_df <- 
        apply(trait_df,
              2,
              function(x){
                x[!is.na(x)] <- 
                  scale(x[!is.na(x)], scale = FALSE);
                x
              })
    } else {
      warning("wrong type input")
    }
    
    if(as.character(sum(weight_df_rescale, na.rm = TRUE)) != 
       as.character(nrow(weight_df_rescale))){
      warning("Sum of weight value is not equal to 1")
      print(paste(as.character(sum(weight_df_rescale, na.rm = TRUE)), "!=", as.character(nrow(weight_df))))
    }
    CWM <- rowSums(trait_df * weight_df_rescale, na.rm = TRUE)
    return(CWM)
  }

# 2.5. CWM_generator_table----
CWM_generator_table <- function(type = NULL){
  trait_CWM <- 
    matrix(ncol = 7, nrow = 27) %>%
    as.data.frame() %>%
    `colnames<-`(c("THICK","CHL", "log_LA", "log_SLA", "LDMC",
                   "EWT", "WOOD_DEN"))
  
  for(trait in names(trait_cube)){
    trait_df <- trait_cube[[trait]]
    if(trait == "WOOD_DEN"){
      weight_df <- weight_IVI_wood
    } else {
      weight_df <- weight_IVI_leaf
    }
    trait_CWM[, trait] <- 
      CWM_generator(trait_df = trait_df,
                    weight_df = weight_df,
                    type = type)
  } # trait loop
  return(trait_CWM)
}

# 2.6. gg_align_plot----
gg_align_plot <- function(gg_plotlist = NULL,
                          layout){
  GG_list_grob <- list()
  grob_widths <- list()
  grob_text <- ""
  grob_width_text <- ""
  for(i in 1:length(gg_plotlist)){
    gg_label <- 
      ggdraw(gg_plotlist[[i]]) +
      draw_text(paste0("(", letters[i], ")"), x = 0.10, y = 0.95, hjust = 1, vjust = 1, size = 24, fontface = "bold")
    GG_list_grob[[i]] <- ggplotGrob(gg_label)
    grob_widths[[i]] <- GG_list_grob[[i]]$widths
    grob_text <- paste0(grob_text, "GG_list_grob[[", i, "]], ")
    grob_width_text <- paste0(grob_width_text, "grob_widths[[", i, "]], ")
  } # for loop
  
  maxWidth <- eval(parse(text = paste0("grid::unit.pmax(",
                                       substr(grob_width_text, 1, nchar(grob_width_text) - 2),
                                       ")")))
  
  for(i in 1:length(gg_plotlist)){
    GG_list_grob[[i]]$widths <- maxWidth
  } # for loop
  
  final_plot <- 
    eval(parse(text = paste0("grid.arrange(", grob_text, "layout_matrix = layout)")))
  return(final_plot)
} # function

# 2.7. Abbreviation generator----
# 2.7.1. Environmental variables----
env_abbr <- function(x){
  abbr <- 
    switch(x,
           "freq_jan" = "expression(FF[Jan])",
           "freq_jul" = "expression(FF[Jul])",
           "elevation" = "expression(Elevation)",
           "temp_jan" = "expression(T[Jan])",
           "temp_jul" = "expression(T[Jul])",
           "RH_jan" = "expression(RH[Jan])",
           "RH_jul" = "expression(RH[Jul])",
           "precip_jan" = "expression(PP[Jan])",
           "precip_jul" = "expression(PP[Jul])",
           "pH" = "expression(pH)",
           "N" = "expression(N)",
           "log_P" = "expression(log[10](P))",
           "log_K" = "expression(log[10](K))",
           "freezing" = "expression(Subfreezing)")
  return(abbr)
}

# 2.7.2. trait variables----
trait_abbr <- function(trait, type){
  trait_type <- paste0(type, trait)
  abbr <- 
    switch(trait_type,
           "speTHICK" = "expression(CWM[S](Lth))",
           "speCHL" = "expression(CWM[S](CHL))",
           "spelog_LA" = "expression(CWM[S](LA[log10]))",
           "spelog_SLA" = "expression(CWM[S](SLA[log10]))",
           "speLDMC" = "expression(CWM[S](LDMC))",
           "speEWT" = "expression(CWM[S](EWT))",
           "speWOOD_DEN" = "expression(CWM[S](WD))",
           "fixedTHICK" = "expression(CWM[F](Lth))",
           "fixedCHL" = "expression(CWM[F](CHL))",
           "fixedlog_LA" = "expression(CWM[F](LA[log10]))",
           "fixedlog_SLA" = "expression(CWM[F](SLA[log10]))",
           "fixedLDMC" = "expression(CWM[F](LDMC))",
           "fixedEWT" = "expression(CWM[F](EWT))",
           "fixedWOOD_DEN" = "expression(CWM[F](WD))",
           "intraTHICK" = "expression(CWM[ITV](Lth))",
           "intraCHL" = "expression(CWM[ITV](CHL))",
           "intralog_LA" = "expression(CWM[ITV](LA[log10]))",
           "intralog_SLA" = "expression(CWM[ITV](SLA[log10]))",
           "intraLDMC" = "expression(CWM[ITV](LDMC))",
           "intraEWT" = "expression(CWM[ITV](EWT))",
           "intraWOOD_DEN" = "expression(CWM[ITV](WD))")
  return(abbr)
}
# 2.7.3. plant type dominance----
plant_domin_abbr <- function(x){
  abbr <- 
    switch(x,
           "C_BA" = "expression(BA[C])",
           "D_BA" = "expression(BA[D])",
           "E_BA" = "expression(BA[E])",
           "D_IVI_broad" = "expression(IVI[D/D+E])")
  return(abbr)
}


# 3.	Species trait table ----
{
# 3.1. Import trait----
raw_trait <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Trait.xlsx", 
              sheet = "Trait") %>%
  mutate(latin = gsub(pattern = "[ ]", 
                      replacement = ".",
                      .$latin)) %>%
  group_by(plot, latin) %>%
  summarise(thick_ave = mean(as.numeric(thick_ave)),
            chl_ave = mean(as.numeric(chl_ave)),
            logLA = mean(log10(as.numeric(LA))),
            logSLA = mean(log10(as.numeric(SLA))),
            LDMC = mean(as.numeric(LDMC)),
            EWT = mean(as.numeric(EWT)),
            WOOD_DEN = mean(as.numeric(wood_density), na.rm = TRUE)) 
# The warning message is caused by the NA in wood density.

# Calculate Q table of each trait
THICK <- 
  raw_trait[, c("plot", "latin", "thick_ave")] %>% 
  pivot_wider(names_from = latin, 
              values_from = thick_ave) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()

CHL <- 
  raw_trait[, c("plot", "latin", "chl_ave")] %>% 
  pivot_wider(names_from = latin,
              values_from = chl_ave) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()

LDMC <- 
  raw_trait[, c("plot", "latin", "LDMC")] %>% 
  pivot_wider(names_from = latin, 
              values_from = LDMC) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()

log_LA <-
  raw_trait[, c("plot", "latin", "logLA")] %>% 
  pivot_wider(names_from = latin, 
              values_from = logLA) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()

log_SLA <- 
  raw_trait[, c("plot", "latin", "logSLA")] %>% 
  pivot_wider(names_from = latin, 
              values_from = logSLA) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()

EWT <-
  raw_trait[, c("plot", "latin", "EWT")] %>% 
  pivot_wider(names_from = latin, 
              values_from = EWT) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()
WOOD_DEN <- 
  raw_trait[, c("plot", "latin", "WOOD_DEN")] %>% 
  pivot_wider(names_from = latin, 
              values_from = WOOD_DEN) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(sort(unique(raw_trait$latin))) %>%
  dplyr::select(-which(colSums(is.na(.)) == 27)) %>%
  as.data.frame()
# 3.2. trait cube----
trait_cube <- 
  list(THICK = THICK,
       CHL = CHL,
       log_LA = log_LA,
       log_SLA = log_SLA,
       LDMC = LDMC,
       EWT = EWT,
       WOOD_DEN = WOOD_DEN)
}
# 4.	Species composition table ----
{
raw_composition <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Species%20composition.xlsx", 
              which = "Species_composition")
plant_type <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Species%20composition.xlsx", 
              which = "Plant_type")
}
# 4.1. produce IVI table (for leaf)----
weight_IVI_leaf <- 
  raw_composition %>%
  filter(latin %in% 
           plant_type$latin[plant_type$plant_type != "C"]) %>%
  filter(latin != "Viburnum.urceolatum") %>% # The individual from this species was too small that we would kill it for collection.
  group_by(plot, latin) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = sqrt(sum((dbh ^ 2)/4)) * 2,
          count = n()) %>% 
  unique() %>%
  ungroup() %>%
  group_by(plot) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = dbh,
          count = count,
          relative_BA = dbh ^ 2 / sum(dbh ^ 2),
          relative_abun = count / sum(count),
          IVI = (relative_abun + relative_BA) / 2) %>%
  dplyr::select(plot, latin, IVI) %>%
  pivot_wider(names_from = latin, values_from = IVI) %>%
  as.data.frame() %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(sort(colnames(.)))
# 4.2. produce IVI table (for wood)----
weight_IVI_wood <- 
  raw_composition %>%
  filter(latin %in% 
           plant_type$latin[plant_type$plant_type != "C"]) %>%
  filter(latin != "Viburnum.urceolatum",
         latin != "Photinia.serratifolia") %>% # The individuals from these species were too small that we would kill them for collection. 
  group_by(plot, latin) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = sqrt(sum((dbh ^ 2)/4)) * 2,
          count = n()) %>% 
  unique() %>%
  ungroup() %>%
  group_by(plot) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = dbh,
          count = count,
          relative_BA = dbh ^ 2 / sum(dbh ^ 2),
          relative_abun = count / sum(count),
          IVI = (relative_abun + relative_BA) / 2) %>%
  dplyr::select(plot, latin, IVI) %>%
  pivot_wider(names_from = latin, values_from = IVI) %>%
  as.data.frame() %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(sort(colnames(.)))
# 4.3. produce IVI table (for all woody species)----
weight_IVI_woody <- 
  raw_composition %>%
  group_by(plot, latin) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = sqrt(sum((dbh ^ 2)/4)) * 2,
          count = n()) %>% 
  unique() %>%
  ungroup() %>%
  group_by(plot) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = dbh,
          count = count,
          relative_BA = dbh ^ 2 / sum(dbh ^ 2),
          relative_abun = count / sum(count),
          IVI = (relative_abun + relative_BA) / 2) %>%
  dplyr::select(plot, latin, IVI) %>%
  pivot_wider(names_from = latin, values_from = IVI) %>%
  as.data.frame() %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(sort(colnames(.)))

# 4.4. produce reBA table (for all woody species)----
weight_reBA_woody <- 
  raw_composition %>%
  group_by(plot, latin) %>%
  reframe(plot = plot,
          latin = latin,
          dbh = sqrt(sum((dbh ^ 2)/4)) * 2) %>% 
  unique() %>%
  ungroup() %>%
  group_by(plot) %>%
  reframe(plot = plot,
          latin = latin,
          reBA = dbh ^ 2 / sum(dbh ^ 2)) %>%
  dplyr::select(plot, latin, reBA) %>%
  pivot_wider(names_from = latin, values_from = reBA) %>%
  as.data.frame() %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(sort(colnames(.)))

# 4.6. plant type dominance----
plant_dominance_df <- 
  data.frame(C_BA = rowSums(weight_reBA_woody[, plant_type$latin[which(plant_type$plant_type == "C")]], na.rm = TRUE),
             E_BA = rowSums(weight_reBA_woody[, plant_type$latin[which(plant_type$plant_type == "E")]], na.rm = TRUE),
             D_BA = rowSums(weight_reBA_woody[, plant_type$latin[which(plant_type$plant_type == "D")]], na.rm = TRUE),
             D_IVI_broad = rowSums(weight_IVI_woody[, plant_type$latin[which(plant_type$plant_type == "D")]], na.rm = TRUE) / rowSums(weight_IVI_woody[, plant_type$latin[which(plant_type$plant_type != "C")]], na.rm = TRUE)
             )

# 5.	Environmental variables table ----
# 5.1. selected environmental variables----
select_ENV <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Environmental%20variable.xlsx", 
              which = "Environmental_variables") %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot, -x, -y, -survey_date)
# the soil phosphorus in F04A and F09A is much more higher than other plots, we decided to remove them from the analysis.

# 5.2. monthly temperature----
monthly_temp <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Environmental%20variable.xlsx", which = "Monthly_temp_and_fog") %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(contains("temp"))

monthly_temp[, 1:ncol(monthly_temp)] <-
  monthly_temp[, 1:ncol(monthly_temp)] %>%
  apply(2, function(x)as.numeric(x))

monthly_fog <- 
  rio::import("https://github.com/r09b44030/trait-fog-relationship/raw/main/Environmental%20variable.xlsx", which = "Monthly_temp_and_fog") %>%
  `rownames<-`(.$plot) %>%
  dplyr::select(-plot) %>%
  dplyr::select(contains("fog"))

monthly_fog[, 1:ncol(monthly_fog)] <-
  monthly_fog[, 1:ncol(monthly_fog)] %>%
  apply(2, function(x)as.numeric(x))

# 6. Correlations between environmental factors----
# 6.1. Multiple linear regression: temp ~ fog + elevation ----
temp_fog.ele_result_scale <- 
  data.frame(month = NULL,
             fog_slope = NULL,
             fog_p = NULL,
             ele_slope = NULL,
             ele_p = NULL,
             model_p = NULL,
             vif = NULL)
for(i in names(monthly_temp)){
  month <- substr(i, 6, nchar(i))
  print(month)
  df_temp <- 
    data.frame(fog = scale(monthly_fog %>%
                             dplyr::select(contains(month))),
               temp = scale(monthly_temp %>%
                              dplyr::select(contains(month))),
               elevation = scale(select_ENV$elevation
               )) %>%
      `colnames<-`(c("fog", "temp", "elevation"))
  
  lmp_HT_result <- 
    lmp_HT(x1 = df_temp$fog,
           x2 = df_temp$elevation,
           y = df_temp$temp,
           n_perm = 999,
           set_seed = 1211)
  temp_fog.ele_result_scale <- 
    rbind(temp_fog.ele_result_scale, 
          c(month = month,
            fog_slope = 
              lmp_HT_result$summary$slope_x1,
            fog_p = 
              lmp_HT_result$summary$p_x1,
            elevation_slope = 
              lmp_HT_result$summary$slope_x2,
            elevation_p = 
              lmp_HT_result$summary$p_x2,
            model_p = 
              lmp_HT_result$summary$model_p,
            vif = 
              lmp_HT_result$summary$vif))
}
colnames(temp_fog.ele_result_scale) <- 
  c("month", "fog_slope", "fog_p", 
    "elevation_slope", "elevation_p", "model_p", "vif")

# ggplot (Figure 2)
temp_fog.ele_result_scale$month <-
  c("Aug", "Sep", "Oct", "Nov", "Dec",
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul")

temp_fog.ele_result_scale <- 
  temp_fog.ele_result_scale %>%
  mutate(month_fac = factor(month,
                            levels = c("Aug", "Sep", "Oct", "Nov", "Dec",
                                       "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul")),
         fog_slope = as.numeric(fog_slope),
         elevation_slope = as.numeric(elevation_slope),
         fog_p = as.numeric(fog_p),
         elevation_p = as.numeric(elevation_p),
         model_p = as.numeric(model_p),
         vif = as.numeric(vif))

gg_fog_ele_coef <- # the code fail to project the month on x label.
  ggplot(aes(x = month_fac, y = fog_slope, group = 1),
         data = temp_fog.ele_result_scale) + 
  geom_line(aes(x = month_fac, y = fog_slope, group = 1),
            color = "dodgerblue") + 
  geom_line(aes(x = month_fac, y = elevation_slope, group = 1),
            color = "black") +
  geom_point(aes(x = month_fac, y = fog_slope, group = 1),
             color = "dodgerblue",
             pch = 21,
             fill = c("white", "white", "dodgerblue",
                      "white", "dodgerblue", "white",
                      "dodgerblue", "white", "dodgerblue",
                      "dodgerblue", "white", "white"),
             cex = 3) + 
  geom_point(aes(x = month_fac, y = elevation_slope, group = 1),
             color = "black",
             pch = 21,
             fill = c("black", "black", "white",
                      "white", "black", "white",
                      "white", "black", "black",
                      "black", "black", "black"),
             cex = 3) +
  geom_hline(yintercept = 0,
             lty = 2,
             color = "grey") + 
  xlab("Month") +
  ylab("Coefficients") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        panel.background = 
          element_rect(fill = "white",
                       colour = "black",
                       linewidth = 0.5, 
                       linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.8, 0.8, 0.5, 0.5), "cm")
  ) # theme

gg_fog_ele_coef

ggsave("Figure_4.png",
       gg_fog_ele_coef,
       height = 5.78 * 1,
       width = 6.71 * 1,
       unit = "in",
       dpi = 300)

# 6.2. Pairwise correlation for selected variables----
GG_env_env <- list()
# Table S2
select_env_env_result <- 
  data.frame(env1 = NULL,
             env2 = NULL,
             r = NULL,
             p = NULL,
             p_ori = NULL)
{ # the soil phosphorus in F04A and F09A is much more higher than other plots, we decided to remove them from the analysis.
  select_ENV_Pout <- select_ENV
  select_ENV_Pout$log_P[which(rownames(select_ENV_Pout) %in% c("F04A", "F09A"))] <- NA
  select_ENV_Pout[, 2:ncol(select_ENV_Pout)] <-
    select_ENV_Pout[, 2:ncol(select_ENV_Pout)] %>%
    apply(2, function(x)as.numeric(x))
  }
m <- 1
for(i in colnames(select_ENV_Pout)){
  for(j in colnames(select_ENV_Pout)){
    print(paste(m, ":", j, "~", i))
    cor_env_env <- 
      cor.test_HT(scale(select_ENV_Pout[, j]),
                  scale(select_ENV_Pout[, i]),
                  n_perm = 999,
                  hierar = gl(n = 9,
                              k = 3),
                  use = "complete.obs",
                  na.action = "na.omit",
                  set_seed = 1211)
    select_env_env_result <- 
      rbind(select_env_env_result,
            data.frame(env1 = i,
                       env2 = j,
                       r = cor_env_env$r,
                       p = cor_env_env$p_shuff,
                       p_ori = round(cor_env_env$p_ori, 3)),
            make.row.names = FALSE)
    
    # Generate ggplot
    df1 <- data.frame(x = select_ENV_Pout[, i],
                      y = select_ENV_Pout[, j])
    
    ggtitle_string <- 
      eval(parse(text = paste0("expression(paste(r, \" = \", \"", 
                               format(cor_env_env$r, nsmall = 3, trim = FALSE), 
                               "\", \", \", p, \" = \",\"", 
                               format(cor_env_env$p_shuff, nsmall = 3, trim = FALSE), 
                               "\"))")))
    xlab_string <-
      eval(parse(text = env_abbr(i)))
    ylab_string <-
      eval(parse(text = env_abbr(j)))
    
    gg_env_env <- 
      ggplot(aes(x = x, 
                 y = y),
             data = df1) + 
      geom_point(size = 3) +
      ggpmisc::stat_ma_line(na.rm = TRUE,
                            se = FALSE,
                            lwd = ifelse(cor_env_env$p_shuff <= 0.10, 1, 0),
                            linetype = ifelse(cor_env_env$p_shuff <= 0.05, 1, 2)) +
      xlab(xlab_string) +
      ylab(ylab_string) +
      ggtitle(ggtitle_string) +
      ylim(c(min(df1$y) - sd(df1$y),
             max(df1$y) + sd(df1$y))) +
      theme(plot.title = element_text(size = 20,
                                      hjust = 1),
            axis.title = element_text(size = 24),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            panel.background = 
              element_rect(fill = "white",
                           colour = "black",
                           linewidth = 0.5, 
                           linetype = "solid"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) # theme
    if(i == "log_P"){
      df_Pout <- 
        data.frame(x = select_ENV[c(10, 25), i],
                   y = select_ENV[c(10, 25), j])
      gg_env_env <- 
        gg_env_env + 
        geom_point(aes(x = x, y = y),
                   data = df_Pout,
                   color = "grey",
                   size = 3)
    }
    if(j == "log_P"){
      df_Pout <- 
        data.frame(x = select_ENV[c(10, 25), i],
                   y = select_ENV[c(10, 25), j])
      gg_env_env <- 
        gg_env_env + 
        geom_point(aes(x = x, y = y),
                   data = df_Pout,
                   color = "grey",
                   size = 3)
    }
    GG_env_env[[m]] <- gg_env_env
    m <- m + 1
  } # j loop
} # i loop
# Group found with zero error variance.

# ggplot (Figure 3)
gg_align_plot(gg_plotlist = 
                GG_env_env[c(2, 3, 4,
                             19, 6, 21,
                             8, 23, 10,
                             11, 12, 13,
                             14)],
              layout = rbind(c(1, 2, 3),
                             c(4, 5, 6),
                             c(7, 8, 9),
                             c(10, 11, 12),
                             c(13, NA, NA)))
ggsave("Figure_3.png",
       gg_align_plot(gg_plotlist =
                       GG_env_env[c(2, 3, 4,
                                    19, 6, 21,
                                    8, 23, 10,
                                    11, 12, 13,
                                    14)],
                     layout = rbind(c(1, 2, 3),
                                    c(4, 5, 6),
                                    c(7, 8, 9),
                                    c(10, 11, 12),
                                    c(13, NA, NA))),
       height = 5.78 * 5,
       width = 6.71 * 3,
       unit = "in",
       dpi = 300)
# 6.3. Pairwise correlation for environmental variables and plant form dominance----
# Table S3
select_plant_env_result <- 
  data.frame(plant_form = NULL,
             env = NULL,
             r = NULL,
             p = NULL,
             p_ori = NULL)
GG_env_plant <- list()
m <- 1
for(i in colnames(plant_dominance_df)){
  for(j in colnames(select_ENV_Pout)){
    print(paste(m, ":", i, "~", j))
    cor_env_plant <- 
      cor.test_HT(scale(plant_dominance_df[, i]),
                  scale(select_ENV_Pout[, j]),
                  n_perm = 999,
                  hierar = gl(n = 9,
                              k = 3),
                  use = "complete.obs",
                  na.action = "na.omit",
                  set_seed = 1211)
    select_plant_env_result <- 
      rbind(select_plant_env_result,
            data.frame(plant_dominance_df = i,
                       env = j,
                       r = cor_env_plant$r,
                       p = cor_env_plant$p_shuff,
                       p_ori = round(cor_env_plant$p_ori, 3)),
            make.row.names = FALSE)
    
    # Generate ggplot
    df1 <- data.frame(x = select_ENV_Pout[, j],
                      y = plant_dominance_df[, i])
    
    ggtitle_string <- 
      eval(parse(text = paste0("expression(paste(r, \" = \", \"", 
                               format(cor_env_plant$r, nsmall = 3, trim = FALSE), 
                               "\", \", \", p, \" = \",\"", 
                               format(cor_env_plant$p_shuff, nsmall = 3, trim = FALSE), 
                               "\"))")))
  xlab_string <-
    eval(parse(text = env_abbr(j)))
  ylab_string <-
    eval(parse(text = plant_domin_abbr(i)))
    
    
    gg_env_plant <- 
      ggplot(aes(x = x, 
                 y = y),
             data = df1) + 
      geom_point(size = 3) +
      ggpmisc::stat_ma_line(na.rm = TRUE,
                            se = FALSE,
                            lwd = ifelse(cor_env_plant$p_shuff <= 0.10, 1, 0),
                            linetype = ifelse(cor_env_plant$p_shuff <= 0.05, 1, 2)) +
      xlab(xlab_string) +
      ylab(ylab_string) +
      ggtitle(ggtitle_string) +
      theme(plot.title = element_text(size = 20,
                                      hjust = 1),
            axis.title = element_text(size = 24),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            panel.background = 
              element_rect(fill = "white",
                           colour = "black",
                           linewidth = 0.5, 
                           linetype = "solid"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) # theme
    
    GG_env_plant[[m]] <- gg_env_plant
    m <- m + 1
  } # j loop
} # i loop

# 6.3.2. Simple linear regression----
spe_trait_CWM <- 
  CWM_generator_table(type = "spe") %>%
  `rownames<-`(rownames(trait_cube[[1]]))
fixed_trait_CWM <- 
  CWM_generator_table(type = "fixed") %>%
  `rownames<-`(rownames(trait_cube[[1]]))
intra_trait_CWM <- 
  CWM_generator_table(type = "intra") %>%
  `rownames<-`(rownames(trait_cube[[1]]))

df_SLA_D <- 
  data.frame(SLA = fixed_trait_CWM$log_SLA,
             D_broad = plant_dominance_df$D_IVI_broad)

slr_result_SLA_D <- 
  slr_HT_fval(SLA ~ D_broad, data1 = df_SLA_D, set_seed = 1211)

gg_SLA_D <- 
  ggplot(aes(x = D_broad, 
             y = SLA),
         data = df_SLA_D) + 
  geom_point(size = 3) +
  geom_smooth(method = "lm",
              se = ifelse(slr_result_SLA_D$final_p <= 0.1,
                          TRUE,
                          FALSE),
              lwd = ifelse(slr_result_SLA_D$final_p <= 0.1,
                           1,
                           0),
              linetype = ifelse(slr_result_SLA_D$final_p <= 0.05,
                                1,
                                2)) + 
  xlab(expression(IVI[D/D+E])) +
  ylab(expression(CWM[f](SLA[log10]))) +
  ggtitle(eval(parse(text = paste0("expression(paste(R^2, \" = \", \"", 
                                   format(round(slr_result_SLA_D$lm_ori$adj.r.squared, 3), nsmall = 3, trim = FALSE), 
                                   "\", \", \", p, \" = \",\"", 
                                   format(slr_result_SLA_D$final_p, nsmall = 3, trim = FALSE), 
                                   "\"))")))) + 
  theme(plot.title = element_text(size = 20,
                                  hjust = 1),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.background = 
          element_rect(fill = "white",
                       colour = "black",
                       linewidth = 0.5, 
                       linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) # theme

gg_align_plot(gg_plotlist = 
                list(GG_env_plant[[56]],
                     GG_env_plant[[43]],
                     gg_SLA_D),
              layout = rbind(c(1, 2, 3)))
# Figure 5
ggsave("Figure_5.png",
       gg_align_plot(gg_plotlist =
                       list(GG_env_plant[[56]],
                            GG_env_plant[[43]],
                            gg_SLA_D),
                     layout = rbind(c(1, 2, 3))),
       height = 5.78 * 1,
       width = 6.71 * 3,
       unit = "in",
       dpi = 300)

# 7. Relationships between traits and environmental factors----
# 7.1. Simple linear regression: trait ~ env----
# Table S1
select_trait_env <- 
  data.frame(trait = NULL,
             env = NULL,
             type = NULL,
             slope = NULL,
             p = NULL,
             p_ori = NULL,
             r2 = NULL,
             scale = NULL)
m <- 1
for(i in names(trait_cube)){
  for(j in names(select_ENV_Pout)){
    for(k in c("spe", "fixed", "intra")){
      print(paste(i, j, k, m, "/",
                  length(trait_cube) * ncol(select_ENV) * 3
      ))
      
      if(i == "WOOD_DEN"){
        weight_df <- weight_IVI_wood
      } else {
        weight_df <- weight_IVI_leaf
      }
      trait_df <- trait_cube[[i]]
      env_vec <- select_ENV_Pout[, j]
      
      p_max_table <- 
        permute_test(trait_df = trait_df,
                     weight_df = weight_df,
                     env_vec = env_vec,
                     type = k,
                     n_perm = 999,
                     hierar = gl(n = 9, k = 3),
                     scale = TRUE,
                     set_seed = 1211)
      select_trait_env <- 
        rbind(select_trait_env,
              data.frame(trait = i,
                         env = j,
                         type = k,
                         slope = p_max_table$Coef,
                         p = p_max_table$p_max,
                         p_ori = p_max_table$p_ori,
                         r2 = p_max_table$Rsquare,
                         scale = p_max_table$scale),
              make.row.names = FALSE)
      m <- m + 1
    } # k
  } # j
} # i
colnames(select_trait_env) <- c("trait", "env",
                                "type", "slope",
                                "p", "p_ori",
                                "r2", "scale")
## Figure 6, 7, S1, S2, S3----
GG_CWM_select_ENV <- list()
m <- 1
for(ev in unique(select_trait_env$env)){
  for(tp in unique(select_trait_env$type)){
    for(tr in unique(select_trait_env$trait)){
      print(paste(tp, tr, ev, m))
      if(tp == "spe"){
        df_gg_CWM_fog <-
          data.frame(CWM = spe_trait_CWM[, tr],
                     env = select_ENV_Pout[, ev])
      } else if(tp == "fixed"){
        df_gg_CWM_fog <-
          data.frame(CWM = fixed_trait_CWM[, tr],
                     env = select_ENV_Pout[, ev])
      } else if(tp == "intra"){
        df_gg_CWM_fog <-
          data.frame(CWM = intra_trait_CWM[, tr],
                     env = select_ENV_Pout[, ev])
      }
      sig_p <- 
        select_trait_env %>%
        filter(type == tp & 
                 trait == tr &
                 env == ev) %>%
        dplyr::select(p) %>%
        as.numeric()
      R2 <- 
        select_trait_env %>%
        filter(type == tp & 
                 trait == tr &
                 env == ev) %>%
        dplyr::select(r2) %>%
        as.numeric()
      
      if(tp == "intra"){
        ggtitle_string <- 
          eval(parse(text = paste0("expression(paste(R^2, \" = \", \"", 
                                   format(R2, nsmall = 3, trim = FALSE), 
                                   "\", \", \", p, \" = \",\"", 
                                   format(sig_p, nsmall = 3, trim = FALSE), 
                                   "\"))")))
      } else {
        ggtitle_string <- 
          eval(parse(text = paste0("expression(paste(R^2, \" = \", \"", 
                                   format(R2, nsmall = 3, trim = FALSE), 
                                   "\", \", \", p[max], \" = \",\"", 
                                   format(sig_p, nsmall = 3, trim = FALSE), 
                                   "\"))")))
      }
      xlab_string <-
        eval(parse(text = env_abbr(ev)))
      ylab_string <-
        eval(parse(text = trait_abbr(trait = tr,
                                     type = tp)))
      
      gg_CWM_select_ENV <- 
        ggplot(aes(x = env, y = CWM),
               data = df_gg_CWM_fog) + 
        geom_point(size = 4) +
        geom_smooth(method = "lm",
                    se = ifelse(sig_p <= 0.1,
                                TRUE,
                                FALSE),
                    lwd = ifelse(sig_p <= 0.1,
                                 1,
                                 0),
                    linetype = ifelse(sig_p <= 0.05,
                                      1,
                                      2)) + 
        xlab(xlab_string) +
        ylab(ylab_string) +
        ggtitle(ggtitle_string) +
        theme(plot.title = element_text(size = 20,
                                        hjust = 1),
              axis.title = element_text(size = 24),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              panel.background = 
                element_rect(fill = "white",
                             colour = "black",
                             linewidth = 0.5, 
                             linetype = "solid"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.margin = unit(c(1, 1, 1, 1), "cm")) # theme
      
      if(ev == "log_P"){
        if(tp == "spe"){
          df_gg_CWM_fog_Pout <-
            data.frame(CWM = spe_trait_CWM[c(10, 25), tr],
                       env = select_ENV[c(10, 25), ev])
        } else if(tp == "fixed"){
          df_gg_CWM_fog_Pout <-
            data.frame(CWM = fixed_trait_CWM[c(10, 25), tr],
                       env = select_ENV[c(10, 25), ev])
        } else if(tp == "intra"){
          df_gg_CWM_fog_Pout <-
            data.frame(CWM = intra_trait_CWM[c(10, 25), tr],
                       env = select_ENV[c(10, 25), ev])
        }
        
        gg_CWM_select_ENV <- 
          gg_CWM_select_ENV + 
          geom_point(aes(x = env, y = CWM),
                     data = df_gg_CWM_fog_Pout,
                     color = "grey",
                     size = 4)
      } # if
      GG_CWM_select_ENV[[m]] <- gg_CWM_select_ENV
      m <- m + 1
    }
  }
}

# Figure
gg_align_plot(gg_plotlist = 
                GG_CWM_select_ENV[c(53, 74, 95, 116, 
                                    137, 158, 179, 200, 
                                    221, 242, 263, 284)],
              layout = rbind(c(1, 2, 3, 4),
                             c(5, 6, 7, 8),
                             c(9, 10, 11, 12))
              )

ggsave("Figure_6.png",
       gg_align_plot(gg_plotlist =
                       GG_CWM_select_ENV[c(11, 18 ,19)],
                     layout = rbind(c(1, 2, 3))),
       height = 5.78 * 1,
       width = 6.71 * 3,
       unit = "in",
       dpi = 300)

ggsave("Figure_7.png",
       gg_align_plot(gg_plotlist =
                       GG_CWM_select_ENV[c(284, 249, 291,
                                           285, 250, 292)],
                     layout = rbind(c(1, 2, 3),
                                    c(4, 5, 6))),
       height = 5.78 * 2,
       width = 6.71 * 3,
       unit = "in",
       dpi = 300)

ggsave("Figure_S1.png",
       gg_align_plot(gg_plotlist =
                       GG_CWM_select_ENV[c(4, 5, 25, 26,
                                           11, 12, 32, 33,
                                           18, 19, 39, 40)],
                     layout = rbind(c(1, 2, 3, 4),
                                    c(5, 6, 7, 8),
                                    c(9, 10, 11, 12))
       ),
       height = 5.78 * 3,
       width = 6.71 * 4,
       unit = "in",
       dpi = 300)

ggsave("Figure_S2.png",
       gg_align_plot(gg_plotlist =
                       GG_CWM_select_ENV[c(60, 81, 61, 82,
                                           102, 123, 103, 124,
                                           144, 165, 145, 166,
                                           186, 207, 187, 208,
                                           228, 249, 229, 250,
                                           270, 291, 271, 292)],
                     layout = rbind(c(1, 2, 3, 4),
                                    c(5, 6, 7, 8),
                                    c(9, 10, 11, 12),
                                    c(13, 14, 15, 16),
                                    c(17, 18, 19, 20),
                                    c(21, 22, 23, 24))
       ),
       height = 5.78 * 6,
       width = 6.71 * 4,
       unit = "in",
       dpi = 300)

ggsave("Figure_S3.png",
       gg_align_plot(gg_plotlist =
                       GG_CWM_select_ENV[c(53, 74, 95, 116,
                                           137, 158, 179, 200,
                                           221, 242, 263, 284)],
                     layout = rbind(c(1, 2, 3, 4),
                                    c(5, 6, 7, 8),
                                    c(9, 10, 11, 12))
       ),
       height = 5.78 * 3,
       width = 6.71 * 4,
       unit = "in",
       dpi = 300)


