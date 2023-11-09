library(ggplot2)
library(plyr)


theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



p <- ggplot() + 
  geom_violin(data = resample.pdf,aes(x = period,y = borealization_index), adjust = 1.5)

p_build <- ggplot2::ggplot_build(p)$data

#This comes directly from the source of geom_violin
p_build <- transform(p_build,
                     xminv = x - violinwidth * (x - xmin),
                     xmaxv = x + violinwidth * (xmax - x))

p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y),
                 plyr::arrange(transform(p_build, x = xmaxv), -y))

#Add our fill variable
p_build$fill_group <- ifelse(p_build$y < -1,'Arctic',
                             ifelse(p_build$y > 2, 'Boreal', 'Middle'))
#This is necessary to ensure that instead of trying to draw
# 5 polygons, we're telling ggplot to draw 15 polygons

p_build$group1 <- with(p_build,interaction(factor(group),factor(fill_group)))

max(p_build$y[p_build$group1 == "5.Boreal"])
max(resample.pdf$borealization_index[resample.pdf$period == "1.5_to_2.0"])

unique(p_build$group1)

p_fill <- ggplot() + 
  geom_polygon(data = p_build,
               aes(x = x, y = y,group = group1,fill = fill_group, lty = NULL), alpha = 0.5) +
  coord_flip(ylim = c(-4.5, 4.5)) +
  scale_fill_manual(values = c(cb[3], "red", cb[1])) +
  scale_y_continuous(breaks = seq(-4,4, by = 2))

p_fill


###

p_fill <- ggplot() + 
  geom_polygon(data = p_build,
      aes(x = x,y = y,group = group1,fill = fill_group, lty = NULL), alpha = 0.5) +
  coord_flip(ylim = c(-4.5, 5)) +
  scale_fill_manual(values = c(cb[3], "red", cb[1])) +
  xlab("North Pacific warming") +
  ylab("Borealization index") +
  scale_x_continuous(labels = c("Preindustrial",
                                  "1950 to 0.5°",
                                  "0.5° to 1.0°",
                                  "1.0° to 1.5°",
                                  "1.5° to 2.0°"), breaks = c(1:5)) +
  scale_y_continuous(breaks = seq(-4,4, by = 2)) +
  geom_hline(yintercept = c(-1,2), color = c(cb[3], "red"), lty = 2, lwd = 0.8) +
  theme(legend.position = "none")


  

p_fill

