### combine composite indicators into single graph


library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)


base_year <- 1970


species_names <- readRDS("data/species_names.rds")


all_inds <- readRDS("data/all_socb_goal_indices.rds")


all_smoothed_indices <- readRDS("data/socb_smoothed_indices.rds") %>%
  left_join(.,species_names,
            by = "speciesID")

all_composites <- readRDS("output/annual_status_combine.rds")

published_groups <- readRDS("data/published_groups.rds")
species_groups <- readRDS("data/species_groups.rds")


groupIDs <- species_groups %>%
  filter(groupName %in% published_groups) %>%
  select(groupName,groupID) %>%
  distinct()


all_composites_out <- all_composites %>%
  inner_join(.,groupIDs) %>%
  rename(log_scale_indicator = mean,
         log_scale_indicator_sd = sd,
         log_scale_indicator_lci = q2_5,
         log_scale_indicator_uci = q97_5) %>%
  relocate(groupName,groupID,subgroupID,year,
           log_scale_indicator,log_scale_indicator_sd,log_scale_indicator_lci,log_scale_indicator_uci,
           percent_diff, percent_diff_lci, percent_diff_uci) %>%
  select(-c(ess_bulk,rhat,yearn2,model,variable)) %>%
  mutate(across(log_scale_indicator:percent_diff_uci,
                .fns = ~signif(.x,digits = 4)))

# Share composite_indicators_all.csv with Catherine for upload ------------
write_csv(all_composites_out,"composite_indicators_all.csv")



# Group-level models ------------------------------------------------------

high_level_groups <- groupIDs %>%
  filter(grepl(" All",groupName))


main_composites <- all_composites %>%
  filter(groupName %in% high_level_groups$groupName)


final_years <- main_composites %>%
  group_by(groupName) %>%
  summarise(last_year = max(year,na.rm = TRUE))

names_plot <- main_composites %>%
  inner_join(.,final_years,
             by = c("groupName",
                    "year" = "last_year")) %>%
  mutate(lbl = paste(groupName,round(percent_diff),"%")) %>%
  filter(row_number() < length(high_level_groups)+1)




brks_pch <- c(-98,-95,-90,-75,-50,-33,0,50,100,300,500,1000,2000,5000)
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis



overview <- ggplot(data = main_composites,
                   aes(x = year,y = mean, group = groupName))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = q2_5,ymax = q97_5,
                  fill = groupName),
              alpha = 0.2)+
  geom_line(aes(colour = groupName))+
  geom_text_repel(data = names_plot,
            aes(label = lbl),nudge_x = 10,
            size = 4,
            segment.alpha = 0.3)+
  coord_cartesian(xlim = c(1970,2040),
                  ylim = c(brks_log[3],brks_log[9]))+
  scale_color_viridis_d(aesthetics = c("fill","colour"))+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs)+
  theme(legend.position = "none")

overview


pdf("figures/Suggested_highlevel_composite_indicators.pdf",
    height = 8.5,
    width = 11)
print(overview)
dev.off()






