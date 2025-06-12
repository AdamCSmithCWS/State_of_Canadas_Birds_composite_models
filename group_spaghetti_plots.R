


library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)



base_year <- 1970

# source("functions/specid_rename.R")
#

species_names <- readRDS("data/species_names.rds")


all_inds <- readRDS("data/all_socb_goal_indices.rds")


all_smoothed_indices <- readRDS("socb_smoothed_indices.rds") %>%
  left_join(.,species_names,
            by = "speciesID")

all_composites <- readRDS("output/annual_status_combine.rds")

species_groups <- readRDS("data/species_groups.rds")

groupIDs <- species_groups %>%
  select(groupName,groupID) %>%
  distinct()



all_composites_out <- read_csv("output/composite_indicators_all.csv")


# Example plots for methods -----------------------------------------------
all_groups <- species_groups %>%
  select(groupName,groupID,subgroupID) %>%
  distinct()

main_groups <- species_groups %>%
  filter(subgroupID == 0,
         groupName != "Galliformes: All") %>%
  select(groupName) %>%
  distinct() %>%
  #mutate(groupName = str_trim(gsub(": All","",x = groupName))) %>%
  unlist()

high_level_groups <- main_groups[c(1,2,3,4,5,6,7,8,9,13)]
main_composites <- all_composites %>%
  filter(groupName %in% high_level_groups)

#i <- which(all_groups$groupName == "Aerial Insectivores: Swifts, swallows, and nightjars")

for(i in 1:nrow(all_groups)){
grp <- all_groups[i,"groupName"]
grpID <- all_groups[i,"groupID"]

grp_labl <- species_groups %>%
  filter(groupName == grp) %>%
  select(groupName,groupNameFr) %>%
  distinct() %>%
  unlist()
grp_labl <- gsub("[[:punct:]]","",x = grp_labl)
grp_labl <- gsub("[[:blank:]]","_",x = grp_labl)
grp_labl <- paste(grp_labl,collapse = "-")


if(!file.exists(paste0("output/composite_fit_",grp_labl,".rds"))){next}

annual_status_difference <- readRDS(paste0("output/composite_fit_",grp_labl,".rds"))
inds_all <- readRDS(paste0("output/composite_data_",grp_labl,".rds"))
species_sel <- readRDS(paste0("output/composite_species_list_",grp_labl,".rds"))
inds_w_low_conf <- readRDS(paste0("output/composite_data_w_low_confid_",grp_labl,".rds"))


brks_pch <- c(-98,-95,-90,-75,-50,-33,0,50,100,300,500,1000,2000,5000)
brks_log <- log((brks_pch/100)+1)
brks_labs <- paste0(brks_pch,"%")

ylimu <- max(annual_status_difference$q97_5,1.1)
yliml <- min(annual_status_difference$q2_5,-1.1)




inds_all_plot <- all_smoothed_indices %>%
  filter(speciesID %in% species_sel$speciesID) %>%
  mutate(percent_dif = round((exp(scaled_status)-1)*100))


sp_y <- inds_all_plot %>%
  group_by(speciesID,english_name, french_name) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            mean_prec = 1/(mean(annual_diff_sd)^2),
            .groups = "drop") %>%
  mutate(prec_plot = scale(mean_prec, center = FALSE))



inds_label <- inds_all_plot %>%
  left_join(.,sp_y,
            by = c("speciesID","english_name","french_name")) %>%
  filter(year == last_year) %>%
  mutate(qual_dif = ifelse(percent_dif > 0,1,-1),
         lbl = paste0(english_name,":",percent_dif,"%"),
         lblf = paste0(french_name,":",percent_dif,"%"))

qual_difs <- inds_label %>%
  select(speciesID,qual_dif,prec_plot)

inds_all_plot <- inds_all_plot %>%
  left_join(qual_difs)


ylimu_spag <- max(inds_all_plot$scaled_status)
yliml_spag <- min(inds_all_plot$scaled_status)




tst2 <- ggplot(data = annual_status_difference,
               aes(x = year,y = mean))+
  geom_hline(yintercept = 0)+
  geom_line(data = inds_all_plot,
            aes(x = year,y = scaled_status,
                group = speciesID,
                colour = qual_dif,
                alpha = prec_plot),
            #alpha = 0.3,
            inherit.aes = FALSE)+
  geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
              alpha = 0.2)+
  geom_line()+
  ggrepel::geom_text_repel(data = inds_label,
                           aes(x = year,y = scaled_status,
                               label = lbl,
                               colour = qual_dif),
                           size = 1,
                           max.overlaps = 30,
                           min.segment.length = 0,
                           nudge_x = 6,
                           alpha = 1,
                           box.padding = 0.1,
                           segment.alpha = 0.3,
                           segment.size = 0.2,
                           hjust = "left")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs,
                     limits = c(yliml_spag,ylimu_spag))+
  scale_x_continuous(limits = c(1970,2045),
                     breaks = seq(1970,2020,by = 10),
                     expand = c(0,0))+
  colorspace::scale_color_continuous_diverging(rev = TRUE,
                                               mid = 0,
                                               n_interp = 3)+
  #labs(title = paste(grp_labl))+
  xlab("")+
  ylab("Percent change since first year")+
  theme_classic()+
  theme(legend.position = "none")

tst2

pdf(paste0("figures/",grp_labl,"example.pdf"),
    width = 8,
    height = 6)
print(tst2)

dev.off()

saveRDS(tst2,paste0("figures/saved_rds/",grp_labl,"example.rds"))

tst1 <- ggplot(data = annual_status_difference,
               aes(x = year,y = mean))+
  geom_hline(yintercept = 0)+
  geom_line(data = inds_all_plot,
            aes(x = year,y = scaled_status,
                group = speciesID,
                colour = qual_dif,
                alpha = prec_plot),
            inherit.aes = FALSE)+
  geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
              alpha = 0.2)+
  geom_line()+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs,
                     limits = c(yliml_spag,ylimu_spag))+
  scale_x_continuous(breaks = seq(1970,2020,by = 10),
                     expand = c(0,0))+
  colorspace::scale_color_continuous_diverging(rev = TRUE,
                                               mid = 0,
                                               n_interp = 3)+
  labs(caption = paste(grp_labl))+
  xlab("")+
  ylab("Percent change since first year")+
  theme_classic()+
  theme(legend.position = "none")

png(paste0("Figures/",grp_labl,"spag1.png"),
    res = 200,
    width = 5,height = 4,units = "in")
print(tst1)
dev.off()

saveRDS(tst1,paste0("figures/saved_rds/",grp_labl,"spag1.rds"))



#
# tst1 <- ggplot(data = annual_status_difference,
#                aes(x = year,y = mean))+
#   geom_hline(yintercept = 0)+
#   geom_line(data = inds_all_plot,
#             aes(x = year,y = scaled_status,
#                 group = speciesID,
#                 colour = qual_dif,
#                 alpha = prec_plot),
#             #alpha = 0.3,
#             inherit.aes = FALSE)+
#   geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
#               alpha = 0.2)+
#   geom_line()+
#   scale_y_continuous(breaks = brks_log,
#                      labels = brks_labs,
#                      limits = c(yliml_spag,ylimu_spag))+
#   scale_x_continuous(breaks = seq(1970,2020,by = 10),
#                      expand = c(0,0))+
#   colorspace::scale_color_continuous_diverging(rev = TRUE,
#                                                mid = 0,
#                                                n_interp = 3)+
#   #labs(title = paste(grp_labl))+
#   xlab("")+
#   ylab("Percent change since first year")+
#   theme_classic()+
#   theme(legend.position = "none")
#
# png(paste0("Figures/",grp_labl,"spag2.png"),
#     res = 200,
#     width = 5,height = 4,units = "in")
# print(tst1)
# dev.off()
#
# saveRDS(tst1,paste0("figures/saved_rds/",grp_labl,"spag2.rds"))
#


tst2 <- ggplot(data = annual_status_difference,
               aes(x = year,y = mean))+
  geom_hline(yintercept = 0)+
  geom_line(data = inds_all_plot,
            aes(x = year,y = scaled_status,
                group = speciesID,
                colour = qual_dif,
                alpha = prec_plot),
            #alpha = 0.3,
            inherit.aes = FALSE)+
  geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
              alpha = 0.2)+
  geom_line()+
  ggrepel::geom_text_repel(data = inds_label,
                           aes(x = year,y = scaled_status,
                               label = lblf,
                               colour = qual_dif),
                           size = 1,
                           max.overlaps = 30,
                           min.segment.length = 0,
                           nudge_x = 6,
                           alpha = 1,
                           box.padding = 0.1,
                           segment.alpha = 0.3,
                           segment.size = 0.2,
                           hjust = "left")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs,
                     limits = c(yliml_spag,ylimu_spag))+
  scale_x_continuous(limits = c(1970,2045),
                     breaks = seq(1970,2020,by = 10),
                     expand = c(0,0))+
  colorspace::scale_color_continuous_diverging(rev = TRUE,
                                               mid = 0,
                                               n_interp = 3)+
  #labs(title = paste(grp_labl))+
  xlab("")+
  ylab("Pourcentage de changement depuis la première année")+
  theme_classic()+
  theme(legend.position = "none")

tst2

pdf(paste0("figures/",grp_labl,"_example_Fr.pdf"),
    width = 8,
    height = 6)
print(tst2)

dev.off()
saveRDS(tst2,paste0("figures/saved_rds/",grp_labl,"_example_Fr.rds"))



}








