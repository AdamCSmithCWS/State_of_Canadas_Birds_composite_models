### combine composite indicators into single graph


library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)

source("functions/specid_rename.R")


base_year <- 1970



sp_tbl <- readRDS("data/SocbSpecies.rds") %>%
  rename_with(.,.fn = specid_rename)
group_tbl <- readRDS("data/Groups.rds") %>%
  rename_with(.,.fn = specid_rename)
trend_tbl <- readRDS("data/Trends.rds") %>%
  rename_with(.,.fn = specid_rename)
rank_tbl <- readRDS("data/SocbTrendRank.rds") %>%
  rename_with(.,.fn = specid_rename)
indices_tbl <- readRDS("data/TrendsIndices.rds") %>%
  rename_with(.,.fn = specid_rename)
goal_indices_tbl <- readRDS("data/TrendsIndicesGoals.rds") %>%
  rename_with(.,.fn = specid_rename)
species_names <- readRDS("data/species_names.rds") %>%
  rename_with(.,.fn = specid_rename)


all_inds <- readRDS("data/all_socb_goal_indices.rds")


all_smoothed_indices <- readRDS("socb_smoothed_indices.rds") %>%
  rename_with(.,
              .fn = specid_rename) %>%
  left_join(.,species_names,
            by = "speciesID")

all_composites <- readRDS("output/annual_status_combine.rds")

species_groups <- readRDS("data/species_groups.rds")

groupIDs <- species_groups %>%
  select(groupName,groupID) %>%
  distinct()


all_composites_out <- all_composites %>%
  left_join(.,groupIDs) %>%
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



species_groups <- readRDS("data/species_groups.rds")

main_groups <- species_groups %>%
  filter(subgroupID == 0,
         groupName != "Galliformes: All") %>%
  select(groupName) %>%
  distinct() %>%
  #mutate(groupName = str_trim(gsub(": All","",x = groupName))) %>%
  unlist()

# Group-level models ------------------------------------------------------


out_composites <- all_composites %>%
  select(groupName,year,mean,q2_5,q97_5,percent_diff,percent_diff_lci,percent_diff_uci) %>%
  mutate(across(mean:percent_diff_uci,
                .fns = ~signif(.x,digits = 4)))
# write_csv(out_composites,
#           "output/saved_draft_composite_trajectories.csv")
#

high_level_groups <- main_groups[c(1,2,3,4,5,6,7,8,9,13)]
main_composites <- all_composites %>%
  filter(groupName %in% high_level_groups)


final_years <- main_composites %>%
  group_by(groupName) %>%
  summarise(last_year = max(year))

names_plot <- main_composites %>%
  inner_join(.,final_years,
             by = c("groupName",
                    "year" = "last_year")) %>%
  mutate(lbl = paste(groupName,round(percent_diff),"%")) %>%
  filter(row_number() < length(high_level_groups)+1)




brks_pch <- c(-95,-90,-75,-50,-40,-25,0,25,50,100,200,300,500) # values of  percent change I’d like to show on the y-axis
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis



overview <- ggplot(data = main_composites,
                   aes(x = year,y = mean, group = groupName,
                       colour = groupName))+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_text_repel(data = names_plot,
            aes(label = lbl),nudge_x = 10,
            size = 4,
            segment.alpha = 0.3)+
  coord_cartesian(xlim = c(1970,2040),
                  ylim = c(brks_log[3],brks_log[9]))+
  scale_color_viridis_d()+
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






