library(tidyverse)
for (cutval in seq(1,5)){
  data_nicola_connection <- readr::read_csv(file = paste("./numberConnections_cuotff", cutval,".csv",sep="")) %>%
    dplyr::mutate(id = 1:n()) %>%
    tidyr::gather(group2, count, con.group1:con.group6) %>% group_by(id) %>% arrange(id)%>%
    dplyr::mutate(weight = 1/sum(count, na.rm = TRUE)) %>% filter(weight != Inf & count != 0) %>%
    # dplyr::mutate(weight = 1) %>% filter(weight != Inf & count != 0) %>%
    dplyr::mutate(group2 = as.numeric(stringr::str_extract(group2, "[0-9]"))) %>%
    dplyr::group_by(group, group2, id) %>%
    dplyr::group_by(age, group, group2, id, weight, count) %>%
    dplyr::summarize(count2 = purrr::map(count, ~rep(1, .x))) %>%
    tidyr::unnest() %>%
    dplyr::select(age, group, group2, id, weight, count)
  
  require(foreign)
  require(nnet)
  library(ggplot2)
  library(colorblindr)


  data_nicola_connection<- data_nicola_connection %>% group_by(group2) %>% dplyr::mutate(weight = n()) %>% group_by(group2, weight) %>% nest() %>% dplyr::mutate(weight = sum(weight)/weight) %>% dplyr::mutate(weight = weight /sum(weight)) %>% unnest() %>% mutate(group2 = as.character(group2), group = as.character(group)) %>% arrange(group, id)
  
  modelled <- multinom(data =  data_nicola_connection,
                       formula =  group2 ~ group - 1,
                       weights = weight,
                       maxit = 1000, Hess = TRUE)
  
  dictionary<- c("pre-school", "school", "adults unknown", "adults with kids", "adults no kids", "elderly")
  gdfs
  # save(plotting_data_probs,file = "flubasel_cutoff2_plotdata.Rda")
  plotting_data_probs <- right_join(data_nicola_connection %>%
                      group_by(group) %>%
                      dplyr::summarize(tot_nobs = n()) %>%
                      dplyr::mutate(group = factor(as.numeric(group), levels = 1:6, labels = dictionary, ordered = TRUE)),
                    cbind(data_frame(group = 1:6),
                          predict(modelled, newdata = data_frame(group = as.character(1:6)),
                                  type = "prob")) %>%
                      gather(group2, count,
                             as.character(1):as.character(6)) %>%
                      dplyr::mutate(group = factor(as.numeric(group),
                                                   levels = 1:6, labels = dictionary, ordered = TRUE),
                                    group2 = factor(as.numeric(group2), levels = 1:6,
                                                    labels = dictionary, ordered = TRUE))  %>%
                      group_by(group) %>%
                      dplyr::mutate(ods = (count/(1-count))) %>%
                      group_by(group2) %>%
                      dplyr::mutate(lods_rat = log(ods/ ods[group == group2]))) %>%
    dplyr::mutate(pred_count = tot_nobs * count) %>%
    dplyr::mutate(wilson = purrr::map2(pred_count, tot_nobs,
                                       ~ Hmisc::binconf(.x, .y, alpha = .05, method = "wilson") %>%
                                         as_data_frame())) %>%
    unnest()
  
  p1 = plotting_data_probs %>%
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
  ggsave(plot=p1, file=paste("./ConnectionsBetweenGroups_altfigure_cuotff", cutval,".pdf",sep=""))
  
  
  p2 = plotting_data_probs %>%
    ggplot(aes(x = 1, y = PointEst, group = paste(group,group2), colour = group)) +
    geom_point(position = position_dodge(width = 1), shape = 15, size = 3)  +
    geom_errorbar(aes(ymax = Upper, ymin = Lower), position = position_dodge(width = 1)) +
    # theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(group2~group, switch = "x") +
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
    labs(x = NULL, y = "Probability", title = paste("95% CI for multinomial probability", "cutoff", cutval), subtitle = "Probability of your cluster containing group Y given you belong to group X")
  # plot(plotting_data_probs)
  # require(grid)
  # grid.draw(plotting_data_probs)
  ggsave(plot=p2, file=paste("./ConnectionsBetweenGroups_cuotff", cutval,".pdf",sep=""))
  # ggsave(plot=p2, file=paste("./ConnectionsBetweenGroups_cuotff", cutval,".png",sep=""))
}
#
# cbind(data_frame(group = 1:6),
#       predict(modelled, newdata = data_frame(group = as.character(1:6)),
#               type = "prob", se = TRUE)) %>%
#   gather(group2, count,
#          as.character(1):as.character(6)) %>%
#   dplyr::mutate(group = factor(as.numeric(group), levels = 1:6, labels = dictionary, ordered = TRUE), group2 = factor(as.numeric(group2), levels = 1:6, labels = dictionary, ordered = TRUE)) %>%
#   ggplot(aes(x = group2, y = count, group = group, colour = group)) +
#   geom_point(position = position_dodge(width = 1), shape = 15, size = 3)  +
#   theme_bw() +
#   scale_color_OkabeIto() +
#   facet_grid(.~group2, scales = "free") +
#   theme(strip.text = element_text(colour = "black", size = rel(1)),
#         axis.title.y = element_text(size = rel(1.3), face = "bold"),
#         axis.text.y = element_text(colour = "black", size = rel(1.1)),
#         legend.text = element_text(colour = "black", size = rel(1)),
#         legend.box.margin = margin(0,0,0,0,unit = "pt"),
#         legend.margin = margin(0,0,0,0))
#
# cbind(data_frame(group = 1:6),
#       predict(modelled, newdata = data_frame(group = as.character(1:6)),
#               type = "prob")) %>%
#   gather(group2, count,
#          as.character(1):as.character(6)) %>%
#   dplyr::mutate(group = factor(as.numeric(group), levels = 1:6, labels = dictionary, ordered = TRUE), group2 = factor(as.numeric(group2), levels = 1:6, labels = dictionary, ordered = TRUE))  %>% group_by(group) %>% dplyr::mutate(ods = (count/(1-count))) %>% group_by(group2) %>% dplyr::mutate(lods_rat = log(ods/ ods[group == group2]))  %>%
#   ggplot(aes(x = group2, y = count, group = group, colour = group)) +
#   geom_point(position = position_dodge(width = 1), shape = 15, size = 3)  +
#   theme_bw() +
#   scale_color_OkabeIto() +
#   facet_grid(.~group2, scales = "free") +
#   theme(strip.text = element_text(colour = "black", size = rel(1)),
#         axis.title.y = element_text(size = rel(1.3), face = "bold"),
#         axis.text.y = element_text(colour = "black", size = rel(1.1)),
#         legend.text = element_text(colour = "black", size = rel(1)),
#         legend.box.margin = margin(0,0,0,0,unit = "pt"),
#         legend.margin = margin(0,0,0,0))
