# LFDP: Plotting figures ----

## 3D view of LFDP plot ----
tiff ('LFDP_3D_view.tiff', width = 16, height = 9, units = 'cm', res = 600, compression = 'lzw', pointsize = 6)
require (plot3D)
par (mfrow = c(1,2))
persp3D (x = IL$x, y = IL$y, z = IL$z, col = ramp.col(c("white", "gray10")), shade = .6,
         phi = 30, theta = 30, resfac = 1, 
         scale = TRUE, expand = 1.3, border = 'gray', zlim = c(1700, 1790),
         contour = list(nlevels = 12, col = "black", side = c("zmin", "z"), lwd = .5),
         image = TRUE, box = TRUE, clab = "elevation\n[m asl]", colkey = list(length = 0.3, shift = -0.2),
         xlab = expression ("x (W-E)"), ylab = "y (S-N)", zlab = 'elevation',
         main = 'South-eastern view', cex.main = 2)
title (main = 'A', adj = 0, cex.main = 2)
persp3D (x = IL$x, y = IL$y, z = IL$z, col = ramp.col(c("white", "gray10")), shade = .6,
         phi = 30, theta = -60, resfac = 1, 
         scale = TRUE, expand = 1.3, border = 'gray', zlim = c(1700, 1790),
         contour = list(nlevels = 12, col = "black", side = c("zmin", "z"), lwd = 0.5),
         image = TRUE, box = TRUE, clab = "elevation\n[m asl]", colkey = list(length = 0.3, shift = -0.2), 
         xlab = expression ("x (W-E)"), ylab = "y (S-N)", zlab = 'elevation',
         main = 'South-western view', cex.main = 2)
title (main = 'B', adj = 0, cex.main = 2)
dev.off ()

## Plot description (contours etc) ----
rectangles <- data.frame (x1 = rep (c(0,20,40,60,80), times = 5),
                          x2 = rep (c(10,30,50,70,90), times = 5),
                          y1 = rep (c(0,20,40,60,80), each = 5),
                          y2 = rep (c(10,30,50,70,90), each = 5))

library (metR)
tiff ('LFDP-plot-map.tiff', width = 8, height = 8, units = 'cm', res = 600, compression = 'lzw')
ggplot () + 
  geom_contour (data = IL_df, mapping = aes (x = x*10, y = y*10, z = z), show.legend = F, color = 'grey', linewidth = 0.3) + 
  geom_text_contour(data = IL_df, aes(x = x*10, y = y*10, z = z), size = 1.5, skip = 1, color = 'grey20') + 
  geom_rect(data = rectangles, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.3)  +
  theme_bw () +
  theme (panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.03, 0.03), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.03, 0.03), breaks = seq(0, 100, by = 10)) + 
  xlab ('x [m]') + ylab ('y [m]')
dev.off ()


## plot contours - comparison of original (measured in corners of 10 x 10 m subplots) and smoothed ----
tiff ('LFDP-contours-orig-smoothed.tiff', width = 12, height = 8, units = 'cm', res = 600, compression = 'lzw')
ggplot (IL_df, aes (x = x*10, y = y*10, z = z)) + 
  geom_contour (data = pile_sub, aes (x = x*10, y = y*10, z = elevation, color = 'original')) +
  geom_contour (aes (color = 'smoothed')) + 
  theme_bw () +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) + 
  xlab ('x [m]') + ylab ('y [m]') + 
  scale_color_manual (name = 'Elevation contours',
                      values = c('original' = 'grey', 'smoothed' = 'navy'))
dev.off ()


## Map of individuals according to leaf-types ----
tiff ('LFDP-indiv-leaftypes.tiff', width = 16, height = 8, units = 'cm', res = 600, compression = 'lzw')
ggplot (ind_agg, aes (x = qx, y = qy)) +
  #geom_contour (data = pile[1:121,], aes (x = I(x*10), y = I(y*10), z = elevation), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (colour = factor (DBH_class), shape = factor (DBH_class), size = factor (DBH_class)), show.legend = TRUE) + 
  labs (x = 'x [m]', y = 'y [m]') +
  facet_wrap (~leaf_type) +
  theme_bw () + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(1, "lines"), legend.direction = 'horizontal', legend.position = 'bottom') +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq (0, 100, by = 20)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq (0, 100, by = 20)) + 
  scale_color_manual(name = 'DBH class', labels = c('1-5 cm', '5-30 cm', '30-103 cm'), values = c("red", "navy", "darkgreen")) + 
  scale_shape_manual(name = 'DBH class', labels = c('1-5 cm', '5-30 cm', '30-103 cm'), values = c(16, 1, 1)) + 
  scale_size_manual (name = 'DBH class', labels = c('1-5 cm', '5-30 cm', '30-103 cm'), values = c(.5,1,2))
  
dev.off ()

## Veg types and DCA1/DCA2 axes scores ----
xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_veg <- cbind (xy, z = veg_types)
xy_dca1 <- cbind (xy, z = scores (DCA, display = 'si', choices = 1)[,1])
xy_dca2 <- cbind (xy, z = scores (DCA, display = 'si', choices = 2)[,1])

g_veg <- ggplot (xy_veg, aes (x = x, y = y)) +
  #geom_contour (data = pile[1:121,], aes (x = I(x*10), y = I(y*10), z = elevation), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (shape = factor (z), colour = factor (z)), show.legend = TRUE, size = 5) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_color_manual(name = 'Vegetation type', labels = c('ridge', 'east-facing slope', 'valley'), values = RColorBrewer::brewer.pal (n = 3, name = 'Dark2')) +
  scale_shape_manual(name = 'Vegetation type', labels = c('ridge', 'east-facing slope', 'valley'), values = c(16, 17, 15))

g_dca1 <- ggplot (xy_dca1, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = abs (z), fill = factor (sign (z)), shape = factor (sign (z))), show.legend = FALSE) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_color_manual (values = c("gray15", "gray15"), guide = 'none') + 
  scale_fill_manual (values = c("gray15", "white"), guide = 'none') + 
  scale_shape_manual (values = c(21, 21), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

g_dca2 <- ggplot (xy_dca2, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = abs (z), fill = factor (sign (z)), shape = factor (sign (z))), show.legend = FALSE) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_color_manual (values = c("gray15", "gray15"), guide = FALSE) + 
  scale_fill_manual (values = c("gray15", "white"), guide = FALSE) + 
  scale_shape_manual (values = c(21, 21), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

tiff ('LFDP-vegtypes.tiff', width = 16, height = 8, units = 'cm', res = 600, compression = 'lzw')
g_veg
dev.off ()


library (patchwork)
tiff ('LFDP-DCA1-DCA2.tiff', width = 16, height = 8, units = 'cm', res = 600, compression = 'lzw')
g_dca1 + g_dca2 + plot_annotation(tag_levels = 'A')
dev.off ()

## Ordination diagrams - DCA ----
tiff ('LFDP-DCAs.tiff', width = 16, height = 8, units = 'cm', res = 600, compression = 'lzw', pointsize = 8)
#win.metafile ('LFDP-DCAs.wmf', width = 6.3, height = 3.15, pointsize = 8)
par (mfrow = c(1, 2))
ordiplot (DCA, display = 'si', type = 'n')
points (DCA, display = 'si', col = RColorBrewer::brewer.pal (n = 3, name = 'Dark2')[veg_types], pch = c(16, 17, 15)[veg_types]) 
legend ('top', legend = c('ridge', 'east-facing', 'valley'), title = 'Vegetation type', pch = c(16, 17, 15), col = RColorBrewer::brewer.pal (n = 3, name = 'Dark2'), cex = .9, horiz = TRUE, inset = c(0, -.2), xpd = T, bty = 'n')
plot (ef.topo, p.max = 0.05, arrow.mull = 1, cex = .7)
plot (ef.soil, p.max = 0.1, col = 'red', arrow.mul = 1, cex = .7)

title (main = 'A', adj = 0)

spe_names <- make.cepnames (names (com))
spe_names[c(1, 56)] <- c('Chamobtu', 'Chamform')

ordiplot (DCA, display = 'sp', type = 'n')
orditorp (DCA, display = 'sp', labels = spe_names, col = 'grey20', cex = .6, priority = colSums (com), pch = 3, pcex = .3, pcol = 'grey70', font = 3)
title (main = 'B', adj = 0)
dev.off ()

## Plot statistics (richness, BA etc.) ----
bubble_colour <- 'grey15'
xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_specnum <- cbind (xy, `Species\nrichness` = rowSums (com>0))

g_specnum <- ggplot (xy_specnum, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = `Species\nrichness`), show.legend = TRUE, colour = bubble_colour) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_shape_manual (values = c(16), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_noind <- cbind (xy, `Number of\nindividuals` = ind_agg %>% group_by (subplot) %>% summarize (no_ind = n()) %>% select (no_ind) %>% unlist)

g_noind <- ggplot (xy_noind, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = `Number of\nindividuals`), show.legend = TRUE, colour = bubble_colour) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_shape_manual (values = c(16), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_BA <- cbind (xy, `Basal\narea` = ind_agg %>% group_by (subplot) %>% summarize (BA = sum (BA)) %>% select (BA) %>% unlist)

g_BA <- ggplot (xy_BA, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = `Basal\narea`), show.legend = TRUE, colour = bubble_colour) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_shape_manual (values = c(16), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_maxBA <- cbind (xy, `Maximal\nbasal area` = ind_agg %>% group_by (subplot) %>% summarize (maxBA = max (BA)) %>% select (maxBA) %>% unlist)

g_maxBA <- ggplot (xy_maxBA, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = `Maximal\nbasal area`), show.legend = TRUE, colour = bubble_colour) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_shape_manual (values = c(16), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))


xy <- expand.grid (y = seq (5, 95, by = 10), x = seq (5, 95, by = 10))[,2:1]
xy_nobranch <- cbind (xy, `Number of\nbranches` = ind_agg %>% group_by (subplot) %>% summarize (branches = mean (no_branches)) %>% select (branches) %>% unlist)

g_nobranch <- ggplot (xy_nobranch, aes (x = x, y = y)) +
  geom_contour (data = IL_df, aes (x = I(x*10), y = I(y*10), z = z), linewidth = .5, colour = 'grey', linejoin = 'mitre', show.legend = FALSE) + 
  geom_point (aes (size = `Number of\nbranches`), show.legend = TRUE, colour = bubble_colour) + 
  labs (x = 'x [m]', y = 'y [m]') +
  theme_bw () + 
  scale_color_manual (values = c("black"), guide = 'none') + 
  scale_shape_manual (values = c(16), guide = 'none') + 
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1) + 
  scale_x_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand=c(0.01, 0.01), breaks = seq(0, 100, by = 10))

ggg <- g_specnum + g_noind + g_BA + g_maxBA + g_nobranch + plot_layout (ncol = 2) + plot_annotation(tag_levels = 'A')
ggsave (plot = ggg, filename = 'LFDP-plot-statistics.tiff', width = 15, height = 16, units = 'cm', dpi = 600, compression = 'lzw', scale = 1.3)

## Teabag decomposition vs other env. variables ----
env_decomp <- as.data.frame (scale (env_sub[,c('k_combust', 's_combust', 'Fe', 'evergr_BA', 'northeasterness', 'conif_BA')]))

gg_k_fe <- ggplot (tibble (resid_k = resid (lm (k_combust ~ evergr_BA, env_decomp)), resid_fe = resid (lm (Fe ~ evergr_BA, env_decomp))), aes (x = resid_fe, y = resid_k)) + 
  geom_point () + 
  geom_smooth (method = 'lm') + 
  xlab ('Fe (residual)') + 
  ylab ('decomposition rate (residual)') + 
  theme_bw ()

gg_k_ev <- ggplot (tibble (resid_k = resid (lm (k_combust ~ Fe, env_decomp)), resid_ev = resid (lm (evergr_BA ~ Fe, env_decomp))), aes (x = resid_ev, y = resid_k)) + 
  geom_point () + 
  geom_smooth (method = 'lm') + 
  xlab ('evergreen broadleaf BA (residual)') + 
  ylab ('decomposition rate (residual)') + 
  theme_bw ()

gg_s_ne <- ggplot (tibble (resid_s = resid (lm (s_combust ~ conif_BA, env_decomp)), resid_ne = resid (lm (northeasterness ~ conif_BA, env_decomp))), aes (x = resid_ne, y = resid_s)) + 
  geom_point () + 
  geom_smooth (method = 'lm') + 
  xlab ('northeasterness (residual)') + 
  ylab ('stabilization factor (residual)') + 
  theme_bw ()

gg_s_co <- ggplot (tibble (resid_s = resid (lm (s_combust ~ northeasterness, env_decomp)), resid_co = resid (lm (conif_BA ~ northeasterness, env_decomp))), aes (x = resid_co, y = resid_s)) +   geom_point () + 
  geom_smooth (method = 'lm') + 
  xlab ('coniferous BA (residual)') + 
  ylab ('stabilization factor (residual)') + 
  theme_bw ()

library (patchwork)
tiff ('LFDP-teabag-decomp.tiff', width = 16, height = 16, units = 'cm', res = 600, compression = 'lzw')
(gg_k_fe + gg_k_ev) / (gg_s_ne + gg_s_co) + plot_annotation(tag_levels = 'A')
dev.off ()

## Boxplots for vegetation types ----
# Physiognomic differences


