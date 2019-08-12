library(tidyverse) ##or ggplot2 and dplyr
load("flubasel_cutoff2_plotdata.Rda")

plotting_data_probs %>%
  ggplot(aes(x = group2, y = PointEst, group = group, colour = group)) +
  geom_point(position = position_dodge(width = 1), shape = 15, size = 3)  +
  geom_errorbar(aes(ymax = Upper, ymin = Lower), position = position_dodge(width = 1)) +
  theme_bw() +
  scale_color_OkabeIto() +
  facet_grid(.~group2, scales = "free") +
  guides(colour = guide_legend("From Group")) +
  theme(strip.text = element_blank(),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(colour = "black", size = rel(1.25)),
        axis.text.x = element_text(colour = "black", size = rel(1.25)),
        legend.box.margin = margin(0,0,0,0,unit = "pt"),
        legend.margin = margin(0,0,0,0),
        legend.position = "right" ,
        legend.direction = "vertical",
        legend.title.align = .5,
        legend.text = element_text(colour = "black", size = rel(1.25)),
        legend.title = element_text(colour = "black", size = rel(1.5))) +
  labs(x = "To Group", y = "Probability", title = "95% CI for multinomial probability", subtitle = "Probability of being clustered with 'to group' given you belong to 'from group'")

plotting_data_probs %>%
  ggplot(aes(x = 1, y = PointEst, group = paste(group,group2), colour = group)) +
  geom_point(position = position_dodge(width = 1), shape = 15, size = 3)  +
  geom_errorbar(aes(ymax = Upper, ymin = Lower), position = position_dodge(width = 1)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(group2~group, scales = "free", switch = "x") +
  guides(colour = "none") +
  theme(axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(colour = "black", size = rel(1.25)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.box.margin = margin(0,0,0,0,unit = "pt"),
        legend.margin = margin(0,0,0,0),
        legend.position = "right" ,
        legend.direction = "vertical",
        legend.title.align = .5,
        legend.text = element_text(colour = "black", size = rel(1.25)),
        legend.title = element_text(colour = "black", size = rel(1.5))) +
  labs(x = NULL, y = "Probability", title = "95% CI for multinomial probability", subtitle = "Probability of your cluster containing group Y given you belong to group X")
