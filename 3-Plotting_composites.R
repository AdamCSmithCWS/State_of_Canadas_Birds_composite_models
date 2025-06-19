### combine composite indicators into single graph


library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)

library(scico) #colour palettes for science that are vision-variation-friendly


base_year <- 1970


species_names <- readRDS("data/species_names.rds")


all_inds <- readRDS("data/all_socb_goal_indices.rds")


all_smoothed_indices <- readRDS("data/socb_smoothed_indices.rds") %>%
  left_join(.,species_names,
            by = "speciesID")


published_groups <- readRDS("data/published_groups.rds")
species_groups <- readRDS("data/species_groups.rds")


all_composites <- readRDS("output/annual_status_combine.rds")


## this csv file was manually modified to add more meaningful plotting names
# groupIDs <- species_groups %>%
#   filter(groupName %in% published_groups) %>%
#   select(groupName,groupID) %>%
#   distinct()
#
# write_csv(groupIDs,"data/group_plotting_names.csv")

groupIDs <- read.csv("data/group_plotting_names.csv")

all_composites_out <- all_composites %>%
  inner_join(.,groupIDs) %>%
  rename(log_scale_indicator = mean,
         log_scale_indicator_sd = sd,
         log_scale_indicator_lci = q2_5,
         log_scale_indicator_uci = q97_5) %>%
  relocate(groupName,groupID,year,
           log_scale_indicator,log_scale_indicator_sd,log_scale_indicator_lci,log_scale_indicator_uci,
           percent_diff, percent_diff_lci, percent_diff_uci) %>%
  select(-c(ess_bulk,rhat,yearn2,model,variable)) %>%
  mutate(across(log_scale_indicator:percent_diff_uci,
                .fns = ~signif(.x,digits = 4)))

# Share composite_indicators_all.csv with Catherine for upload ------------
write_csv(all_composites_out,"composite_indicators_all.csv")



# Group-level models ------------------------------------------------------

high_level_groups <- groupIDs %>%
  filter(grepl(" All",groupName)) %>%
  mutate(facet = c(1,2,2,3,4,1,4,1,4,3))

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





# Supplemental Spaghetti Plots --------------------------------------------


div_pal <- scico(palette = "romaO",n = 11)[c(3,5,10)]
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
    mutate(qual_dif = ifelse(percent_dif >= 33,"increase","stable"),
           qual_dif = ifelse(percent_dif <= -25,"decrease",qual_dif),
           lbl = paste0(english_name,":",percent_dif,"%"),
           lblf = paste0(french_name,":",percent_dif,"%"),
           lbl = ifelse(grepl("100%",lbl),
                        gsub("100%",">99.9%",lbl,fixed = TRUE),
                        lbl))

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

  stat_sum <- paste("This indicator has",n_inc,"species with Canadian populations that have increased by 33% or more (in blue),",n_stab,"species that have shown little change (in light greeen), and",
                    n_dec,"species that have decreased by 25% or more (in red) since the base year.")

  fig_cap <- paste(paste0("Figure S",i,"."),"Indicator of mean species status for",grp_plot,"(black line with grey ribbon showing 95% CI of the mean), with coloured lines indicating the smoothed annual status of each species that is included in the indicator. The transparency of each species' line represents the mean uncertainty of the estimates of the species' annual rates of population change (more transparent lines reprsent species with higher uncertainty in their annual status and therefore lower weight in the estimation of the mean indicator line)",
                   stat_sum)

  fig_cap <- str_wrap(fig_cap, width = 130)

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









