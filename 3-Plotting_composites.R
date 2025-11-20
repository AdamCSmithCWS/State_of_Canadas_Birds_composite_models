### combine composite indicators into single graph


library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(naturecounts)
source("functions/specid_rename.R") #function to reconcile some alternate versions of the species_id column name in NatureCounts


#Data on species populations#
# status_tbl <- nc_query_table(username = "adam.smith",
#                table = "SocbStatus",
#                timeout = 120) %>%
#   rename_with(., .fn = specid_rename) %>%
#   mutate(speciesID = ifelse(speciesID == 12250,  # correcting the treatment of Western Flycatcher
#                             12231,speciesID))
#
#
# saveRDS(status_tbl,"data/SocbStatus.rds")
status_tbl <- readRDS("data/SocbStatus.rds")
status_join <- status_tbl %>%
  filter(popType == 1) %>%
  select(speciesID,popStatus) %>%
  mutate(qual_dif = ifelse(popStatus %in% c("LI","MI"),
                           "increase",
                           "stable"),
         qual_dif = ifelse(popStatus %in% c("LD","MD"),
                           "decrease",
                           qual_dif),
         qual_dif = ifelse(popStatus %in% c("DD"),
                           "data defficient",
                           qual_dif)) %>%
  distinct()


library(scico) #colour palettes for science that are vision-variation-friendly


base_year <- 1970


species_names <- readRDS("data/species_names.rds")


all_inds <- readRDS("data/all_socb_goal_indices.rds")


all_smoothed_indices <- readRDS("data/socb_smoothed_indices.rds") %>%
  left_join(.,species_names,
            by = "speciesID")


published_groups <- readRDS("data/published_groups.rds")


species_groups <- readRDS("data/species_groups.rds")# %>%
  # select(-populationID) %>%
  # distinct()


tst_sp_gr <- species_groups %>%
  group_by(groupName,speciesID,populationID) %>%
  summarise(n = n(),
            n_include = length(unique(include)))



all_composites <- readRDS("output/annual_status_combine.rds")

# shore <- all_composites %>% filter(groupName == "Shorebirds: All") %>%
#   group_by(year) %>%
#   sample_n(1)
#
# all_composites <- all_composites %>% filter(!groupName == "Shorebirds: All") %>%
#   bind_rows(shore)
#
# water <- all_composites %>% filter(groupName == "Waterfowl: All") %>%
#   group_by(year) %>%
#   sample_n(1)
#
# all_composites <- all_composites %>% filter(!groupName == "Waterfowl: All") %>%
#   bind_rows(water)
# saveRDS(all_composites,"output/annual_status_combine_2.rds")


## this csv file was manually modified to add more meaningful plotting names
# groupIDs <- species_groups %>%
#   filter(groupName %in% published_groups) %>%
#   select(groupName,groupID) %>%
#   distinct()
#
# write_csv(groupIDs,"data/group_plotting_names.csv")

groupIDs <- read.csv("data/group_plotting_names.csv")


high_level_groups <- groupIDs %>%
  filter(grepl(" All",groupName)) %>%
  mutate(facet = c(1,2,2,3,4,1,4,1,4,3))

sub_groups <- species_groups %>%
  select(subgroupID,groupName) %>%
  distinct()

groupIDs <- groupIDs %>%
  mutate(summary_plot = ifelse(groupName %in% high_level_groups$groupName,
                               "Yes","No")) %>%
  left_join(sub_groups,
            by = "groupName")




all_composites_out <- all_composites %>%
  inner_join(.,groupIDs) %>%
  rename(log_scale_indicator = mean,
         log_scale_indicator_sd = sd,
         log_scale_indicator_lci = q2_5,
         log_scale_indicator_uci = q97_5) %>%
  relocate(groupName,groupID,subgroupID,year,
           log_scale_indicator,log_scale_indicator_sd,
           log_scale_indicator_lci,log_scale_indicator_uci,
           percent_diff, percent_diff_lci, percent_diff_uci) %>%
  select(-c(ess_bulk,rhat,yearn2,model,variable)) %>%
  mutate(across(log_scale_indicator:percent_diff_uci,
                .fns = ~signif(.x,digits = 4))) %>%
  relocate(plotting_name,
           groupName,
           groupID,
           subgroupID,
           summary_plot,
           year,
           percent_diff,
           percent_diff_lci,
           percent_diff_uci,
           p_decrease,
           p_increase)

write_csv(all_composites_out,"composite_indicators_all.csv")


# Share composite_indicators_all.csv with Catherine for upload ------------
all_composites_socb_site <- all_composites_out %>%
  select(groupName,groupID,subgroupID,year,
         log_scale_indicator,log_scale_indicator_sd,log_scale_indicator_lci,log_scale_indicator_uci,
         percent_diff,percent_diff_lci,percent_diff_uci)
write_csv(all_composites_socb_site,"composite_indicators_all_socb_upload.csv")



start_years <- all_composites_out %>%
  group_by(plotting_name,
           groupName,
           groupID,
           summary_plot) %>%
  filter(year == min(year)) %>%
  rename(start_year = year) %>%
  ungroup() %>%
  select(groupName,start_year)


final_percent_change <- all_composites_out %>%
  group_by(plotting_name,
           groupName,
           groupID,
           summary_plot) %>%
  filter(year == max(year)) %>%
  select(plotting_name,
         groupName,
         groupID,
         summary_plot,
         year,
         percent_diff,
         percent_diff_lci,
         percent_diff_uci,
         p_decrease,
         p_increase) %>%
  mutate(across(percent_diff:p_increase,~signif(.x,3))) %>%
  left_join(start_years)

write_csv(final_percent_change,
          "final_percent_change_values_groups_subgroups.csv")

# Group-level models ------------------------------------------------------


high_level_groups
# groupName groupID
# 1               Waterfowl: All      33
# 2            Marine Birds: All      34
# 3              Shorebirds: All      35
# 4           Wetland Birds: All      36
# 5           Birds of Prey: All      38
# 6            Forest Birds: All      39
# 7         Grassland Birds: All      40
# 8     Aerial Insectivores: All      41
# 9            Arctic Birds: All      47
# 10 Long-Distance Migrants: All      52

main_composites <- all_composites %>%
  inner_join(high_level_groups,by = "groupName")


final_years <- main_composites %>%
  group_by(groupName) %>%
  summarise(last_year = max(year,na.rm = TRUE))

names_plot <- main_composites %>%
  inner_join(.,final_years,
             by = c("groupName",
                    "year" = "last_year")) %>%
  mutate(lbl = paste(gsub(pattern = ": All",
                          groupName,
                          replacement = "")))


high_level_groups_sorted <- names_plot %>%
  select(groupName,plotting_name,percent_diff) %>%
  mutate(plotting_name = factor(plotting_name),
         plotting_name = fct_reorder(.f = plotting_name,
                                     .x = percent_diff))










# end set up --------------------------------------------------------------


brks_pch <- c(-98,-95,-90,-75,-50,-33,0,50,100,300,500,1000,2000,5000)
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis



overview <- ggplot(data = main_composites,
                   aes(x = year,y = mean, group = groupName))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = q2_5,ymax = q97_5,
                  fill = groupName),
              alpha = 0.1)+
  geom_ribbon(aes(ymin = q10,ymax = q90,
                  fill = groupName),
              alpha = 0.15)+
  geom_ribbon(aes(ymin = q25,ymax = q75,
                  fill = groupName),
              alpha = 0.15)+
  geom_line(aes(colour = groupName))+
  geom_text_repel(data = names_plot,
            aes(label = lbl,
                colour = groupName),nudge_x = 10,
            nudge_y = 0.1,
            size = 3,
            segment.alpha = 0.6)+
  scale_x_continuous(breaks = seq(1970,2020, 10))+
  coord_cartesian(xlim = c(1970,2030))+
  scale_color_viridis_d(aesthetics = c("fill","colour"),
                        end = 0.9)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs)+
  facet_wrap(vars(facet), ncol = 2)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

overview


pdf("figures/Suggested_highlevel_composite_indicators.pdf",
    height = 8,
    width = 7)
print(overview)
dev.off()





# Standard overview plot --------------------------------------------------
brks_pch <- c(-98,-95,-90,-65,-50,-33,-20,0,25,50,100,300,500,1000,2000,5000)
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis

socb_pal <- read_csv("data/SOCB_palette.csv")
SOCB_palette <- socb_pal$hex
names(SOCB_palette) <- socb_pal$groupName


overview <- ggplot(data = main_composites,
                   aes(x = year,y = mean, group = groupName))+
  geom_hline(yintercept = 0)+
  # geom_ribbon(aes(ymin = q2_5,ymax = q97_5,
  #                 fill = groupName),
  #             alpha = 0.1)+
  # geom_ribbon(aes(ymin = q10,ymax = q90,
  #                 fill = groupName),
  #             alpha = 0.15)+
  # geom_ribbon(aes(ymin = q25,ymax = q75,
  #                 fill = groupName),
  #             alpha = 0.15)+
  geom_line(aes(colour = groupName))+
  # geom_text_repel(data = names_plot,
  #                 aes(label = lbl,
  #                     colour = groupName),nudge_x = 10,
  #                 nudge_y = 0.1,
  #                 size = 3,
  #                 segment.alpha = 0.6)+
  scale_x_continuous(breaks = seq(1970,2020, 10))+
  coord_cartesian(xlim = c(1970,2025),
                  ylim = c(log(0.25),log(2)))+
  # scale_color_viridis_d(aesthetics = c("fill","colour"),
  #                       end = 0.9)+
  scale_colour_manual(values = SOCB_palette)+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()#,
        #panel.grid.major.x = element_blank()
        )

overview


pdf("figures/Standard_highlevel_composite_indicators.pdf",
    height = 6,
    width = 7)
print(overview)
dev.off()





# Supplemental Spaghetti Plots --------------------------------------------


div_pal <- scico(palette = "romaO",n = 11)[c(2,5,9)]
names(div_pal) <- c("decrease","stable","increase")

pdf(paste0("figures/supplemental_figures.pdf"),
    width = 8.5,
    height = 11)


for(i in 1:nrow(groupIDs)){
  grp <- groupIDs[i,"groupName"]
  grpID <- groupIDs[i,"groupID"]
  grp_plot <- unlist(groupIDs[i,"plotting_name"])

  grp_labl1 <- species_groups %>%
    filter(groupName == grp) %>%
    select(groupName,groupNameFr) %>%
    distinct() %>%
    unlist() %>%
    unname()
  grp_labl <- gsub("[[:punct:]]","",x = grp_labl1)
  grp_labl <- gsub("[[:blank:]]","_",x = grp_labl)
  grp_labl <- paste(grp_labl,collapse = "-")


  if(!file.exists(paste0("output/composite_fit_",grp_labl,".rds"))){next}

  annual_status_difference <- readRDS(paste0("output/composite_fit_",grp_labl,".rds"))
  inds_all <- readRDS(paste0("output/composite_data_",grp_labl,".rds"))
  species_sel <- readRDS(paste0("output/composite_species_list_",grp_labl,".rds"))
  inds_w_low_conf <- readRDS(paste0("output/composite_data_w_low_confid_",grp_labl,".rds"))





  inds_all_plot <- all_smoothed_indices %>%
    filter(speciesID %in% species_sel$speciesID) %>%
    mutate(percent_dif = round((exp(scaled_status)-1)*100))


  sp_y <- inds_all_plot %>%
    group_by(speciesID,english_name, french_name) %>%
    summarise(first_year = min(year),
              last_year = max(year),
              mean_prec = 1/(mean(annual_diff_sd,na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(prec_plot = scale(mean_prec, center = FALSE))



  inds_label <- inds_all_plot %>%
    left_join(.,sp_y,
              by = c("speciesID","english_name","french_name")) %>%
    filter(year == last_year) %>%
    mutate(lbl = paste0(english_name),
           lblf = paste0(french_name)) %>%
    left_join(status_join, by = "speciesID")

  qual_difs <- inds_label %>%
    select(speciesID,qual_dif,prec_plot)


  n_status <- qual_difs %>%
    group_by(qual_dif) %>%
    summarise(n = n())

  n_dec <- n_status %>%
    filter(qual_dif == "decrease") %>%
    select(n) %>%
    unlist() %>%
    unname()

  n_stab <- n_status %>%
    filter(qual_dif == "stable") %>%
    select(n) %>%
    unlist() %>%
    unname()

  n_inc <- n_status %>%
    filter(qual_dif == "increase") %>%
    select(n) %>%
    unlist() %>%
    unname()

  stat_sum <- paste("This indicator includes",n_inc,"species with Canadian populations that have increased (moderate or large increase, in blue),",
                    n_stab,"species that have shown little change (in light green), and",
                    n_dec,"species that have decreased (moderate or large decrease, in red) over the long-term. Note: the colours representing each species population status categories are based on the species' status assessment in the published species account of The State of Canada's Birds and may differ from the change value at the end-point of the species' lines in the plot, which includes the precision-weighted smoothing in the first-stage of the analysis.")

  fig_cap <- paste(paste0("Figure S",i,"."),"Indicator of mean species percent population change for",grp_plot,"(black line with grey ribbon showing 95% CI of the mean), with coloured lines reflecting the approximate contribution each species makes to the indicator. These coloured lines are the cumulative annual precent change values of each species included in the indicator, and the transparency of each line shows the mean uncertainty of the estimates of the species' annual rates of population change (i.e., more transparent lines represent species with greater uncertainty in their annual status, and therefore lower weight in the estimation of the mean indicator line)",
                   stat_sum)

  fig_cap <- str_wrap(fig_cap, width = 130)

  inds_all_plot <- inds_all_plot %>%
    left_join(qual_difs,
              by = "speciesID")


  ylimu_spag <- max(inds_all_plot$scaled_status)
  yliml_spag <- min(inds_all_plot$scaled_status)

if(nrow(inds_label) > 75){
  sz = 1.75
}else{
  sz = 2.5
}


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
                             size = sz,
                             max.overlaps = 30,
                             min.segment.length = 0,
                             nudge_x = 10,
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
    scale_colour_manual(values = div_pal)+
    labs(caption = fig_cap)+
    xlab("")+
    ylab("Percent change since first year")+
    theme_bw()+
    theme(legend.position = "none",
          plot.caption = element_text(hjust = 0))

  print(tst2)


  saveRDS(tst2,paste0("figures/saved_rds/",grp_labl,"supplement.rds"))


}
dev.off()




# Spaghetti figure for publication ----------------------------------------

i = 20 # swallows swifts, nightjars

grp <- groupIDs[i,"groupName"]
grpID <- groupIDs[i,"groupID"]
grp_plot <- unlist(groupIDs[i,"plotting_name"])

grp_labl1 <- species_groups %>%
  filter(groupName == grp) %>%
  select(groupName,groupNameFr) %>%
  distinct() %>%
  unlist() %>%
  unname()
grp_labl <- gsub("[[:punct:]]","",x = grp_labl1)
grp_labl <- gsub("[[:blank:]]","_",x = grp_labl)
grp_labl <- paste(grp_labl,collapse = "-")


if(!file.exists(paste0("output/composite_fit_",grp_labl,".rds"))){next}

annual_status_difference <- readRDS(paste0("output/composite_fit_",grp_labl,".rds"))
inds_all <- readRDS(paste0("output/composite_data_",grp_labl,".rds"))
species_sel <- readRDS(paste0("output/composite_species_list_",grp_labl,".rds"))
inds_w_low_conf <- readRDS(paste0("output/composite_data_w_low_confid_",grp_labl,".rds"))





inds_all_plot <- all_smoothed_indices %>%
  filter(speciesID %in% species_sel$speciesID) %>%
  mutate(percent_dif = round((exp(scaled_status)-1)*100))


sp_y <- inds_all_plot %>%
  group_by(speciesID,english_name, french_name) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            mean_prec = 1/(mean(annual_diff_sd,na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(prec_plot = scale(mean_prec, center = FALSE))



inds_label <- inds_all_plot %>%
  left_join(.,sp_y,
            by = c("speciesID","english_name","french_name")) %>%
  filter(year == last_year) %>%
  mutate(lbl = paste0(english_name),
         lblf = paste0(french_name)) %>%
  left_join(status_join, by = "speciesID")

qual_difs <- inds_label %>%
  select(speciesID,qual_dif,prec_plot)


n_status <- qual_difs %>%
  group_by(qual_dif) %>%
  summarise(n = n())

n_dec <- n_status %>%
  filter(qual_dif == "decrease") %>%
  select(n) %>%
  unlist() %>%
  unname()

n_stab <- n_status %>%
  filter(qual_dif == "stable") %>%
  select(n) %>%
  unlist() %>%
  unname()

n_inc <- n_status %>%
  filter(qual_dif == "increase") %>%
  select(n) %>%
  unlist() %>%
  unname()

stat_sum <- paste("This indicator includes",n_inc,"species with Canadian populations that have increased (moderate or large increase, in blue),",
                  n_stab,"species that have shown little change (in light greeen), and",
                  n_dec,"species that have decreased (moderate or large decrease, in red) over the long-term. Note: the colours and the species population status categories are those from the published species page of the State of Canada's Birds and therefore may not agree with the %-change value at the the end-point of the species' line in the plot.")

fig_cap <- paste(paste0("Figure S",i,"."),"Indicator of mean species status for",grp_plot,"(black line with grey ribbon showing 95% CI of the mean), with coloured lines indicating the smoothed annual status of each specie, included in the indicator. The transparency of each species' line represents the mean uncertainty of the estimates of the species' annual rates of population change (more transparent lines represent species with greater uncertainty in their annual status, and therefore lower weight in the estimation of the mean indicator line)",
                 stat_sum)

fig_cap <- str_wrap(fig_cap, width = 130)

inds_all_plot <- inds_all_plot %>%
  left_join(qual_difs,
            by = "speciesID")


ylimu_spag <- max(inds_all_plot$scaled_status)
yliml_spag <- min(inds_all_plot$scaled_status)




tst3 <- ggplot(data = annual_status_difference,
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
                           size = 3,
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
  scale_x_continuous(limits = c(1970,2050),
                     breaks = seq(1970,2020,by = 10),
                     expand = c(0,0))+
  scale_colour_manual(values = div_pal)+
  #labs(caption = fig_cap)+
  xlab("")+
  ylab("Percent change since first year")+
  theme_bw()+
  theme(legend.position = "none")

print(tst3)

pdf("Figures/Spaghetti_example.pdf",
    width = 7,
    height = 5)
print(tst3)
dev.off()


saveRDS(tst3,paste0("figures/saved_rds/",grp_labl,"spaghetti_example_for)publication.rds"))



# group status bar graphs -------------------------------------------------

div_pal1 <- c(scico(palette = "romaO",n = 11)[c(2,5,9)],
              scico(palette = "tofino",n = 11)[c(8)],grey(0.5))
names(div_pal1) <- c("Decreased","Little Change","Increased","Status Assigned","Data Defficient")


status_join

status_full <- species_groups %>%
  left_join(status_join,"speciesID") %>%
  filter(include == "Y",
         groupName %in% published_groups)  %>%
  inner_join(high_level_groups_sorted,
             by = "groupName") %>%           #, relationship = "many-to-many"
  mutate(change_category = factor(qual_dif,
                           levels = c("decrease","stable","increase","data defficient","Assessed"),
                           ordered = TRUE,
                           labels = c("Decreased","Little Change","Increased","Data Defficient","Status Assigned")))

group_size <- status_full %>%
  group_by(plotting_name) %>%
  summarise(n_sp = n())%>%
  arrange(plotting_name)


group_size2 <- status_full %>%
  group_by(plotting_name,change_category) %>%
  summarise(n_sp = n()) %>%
  arrange(plotting_name)

status_only <- status_full %>%
  filter(change_category != "Data Defficient") %>%
  mutate(plot_column = "status")

status_data <- status_full %>%
  mutate(qual_dif = ifelse(qual_dif == "data defficient",qual_dif,"Assessed"),
         plot_column = "data",
         change_category = factor(qual_dif,
                                  levels = c("decrease","stable","increase","data defficient","Assessed"),
                                  ordered = TRUE,
                                  labels = c("Decreased","Little Change","Increased","Data Defficient","Status Assigned")))


plot_name_tmp <- paste0((group_size$plotting_name)," (",group_size$n_sp,")")

status_plot_data <- status_only %>%
  bind_rows(status_data) %>%
  mutate(Information = factor(plot_column,
                              levels = c("data","status"),
                              ordered = TRUE,
                              labels = c("All Species\n(percent of species\nin group)","Status assigned\n(percent of species\nwith assigned status)"))) %>%
  inner_join(group_size,by = "plotting_name") %>%
  arrange(plotting_name) %>%
  mutate(plotting_name2 = paste0(plotting_name," (",n_sp,")"),
         plotting_name2 = factor(plotting_name2,
                                ordered = TRUE,
                                levels = plot_name_tmp))

status_plot <- ggplot(data = status_plot_data)+
  geom_bar(aes(fill = change_category,
               y = plotting_name2),
           position = "fill")+
  xlab("Percent of species")+
  ylab("")+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = c("0","25","50","75","100"),
                     expand = expansion())+
  scale_y_discrete(expand = expansion())+
  scale_colour_manual(values = div_pal1,aesthetics = "fill",
                      name = "",
                      guide = guide_legend(reverse = TRUE))+
  facet_wrap(vars(Information),
             strip.position = c("top"),
             axes = "all_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 9),
        panel.spacing.x = unit(1,"cm"))


############## Figure out why there is one extra aerial insectivore and one extra forest bird
status_plot

pdf("Figures/Species_status_summary.pdf",
    width = 7,
    height = 5)
print(status_plot)
dev.off()






# Status and Confidence Summaries -----------------------------------------


all_inds <- readRDS("data/all_socb_goal_indices.rds") %>%
  select(-c(english_name,french_name,taxon_group,scientific_name,yearn))


species_confidence <- all_inds %>%
  select(speciesID,confidence) %>%
  distinct()


rank_tbl <- readRDS("data/SocbTrendRank.rds") %>%
  rename_with(.,.fn = specid_rename)

rank_tbl <- rank_tbl %>%
  select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
         resultsCode,popType,areaCode)%>%
  filter(goalTrend == "Y",
         popType == 1) #%>%
  distinct() %>%
  inner_join(species_confidence)



  table(species_confidence$confidence)





