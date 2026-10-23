#comparing species lists

library(tidyverse)
library(readxl)
library(naturecounts)
library(patchwork)



source("functions/GAM_basis_function_mgcv.R") # function to set up the GAM smooth data for Stan model
source("functions/specid_rename.R") #function to reconcile some alternate versions of the species_id column name in NatureCounts

# Quantile functions to generate 95% CI instead of stan's default 90
q2_5 <- function(x)c(q2_5 = quantile(x,probs = c(0.025),
                                     names = FALSE))
q97_5 <- function(x)c(q97_5 = quantile(x,probs = c(0.975),
                                       names = FALSE))


re_download <- TRUE # set to TRUE if re-downloading the data files from NatureCounts database
## NOTE: if re_download is set to TRUE, this will overwrite the data files.
## IF the goal is to reproduce the estimates published in the 2024 State of Canada's Birds
## re_download must be set to FALSE
##


# extract NAtureCounts data -----------------------------------------

if(re_download){

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbTrendRank",
                 timeout = 120) -> rank_tbl
  saveRDS(rank_tbl,"data_2026/SocbTrendRank.rds")

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbSpecies",
                 timeout = 120) -> sp_tbl
  saveRDS(sp_tbl,"data_2026/SocbSpecies.rds")

  #Information on species' group#
  nc_query_table(username = "adam.smith",
                 table = "Groups",
                 timeout = 120) -> group_tbl
  saveRDS(group_tbl,"data_2026/Groups.rds")

  #Data on all annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndices",
                 timeout = 240) -> indices_tbl

  saveRDS(indices_tbl,"data_2026/TrendsIndices.rds")

  #Data on goal-based annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndicesGoals",
                 timeout = 120) -> goal_indices_tbl
  saveRDS(goal_indices_tbl,"data_2026/TrendsIndicesGoals.rds")

  #Data on species trends #
  nc_query_table(username = "adam.smith",
                 table = "Trends",
                 timeout = 120) -> trends_tbl
  saveRDS(trends_tbl,"data_2026/Trends.rds")


  #Species groups
  #
  species_groups <- nc_query_table(table = "SocbTrendGroups",
                                   username = "adam.smith",
                                   verbose = TRUE)%>%
    rename_with(., .fn = specid_rename) %>%
    mutate(groupName = gsub("\r\n",
                            "",
                            x = groupName,
                            fixed = TRUE),
           groupNameFr = gsub("\r\n",
                              "",
                              x = groupNameFr,
                              fixed = TRUE),
           include = toupper(include))%>%
    mutate(speciesID = ifelse(speciesID == 12250,  # correcting the treatment of Western Flycatcher
                              12231,speciesID)) %>%
    filter(!(groupName == "Wetland Birds: All" & subgroupID != 0))

  # wetland <- species_groups %>%
  #   filter(groupName == "Wetland Birds: All")

  saveRDS(species_groups,"data_2026/species_groups.rds")


  species_names <- naturecounts::search_species() %>%
    rename_with(.,.fn = specid_rename)

  saveRDS(species_names,"data_2026/species_names.rds")


  ### Data on species populations#
  status_tbl <- nc_query_table(username = "adam.smith",
                 table = "SocbStatus",
                 timeout = 120) %>%
    rename_with(., .fn = specid_rename) %>%
    mutate(speciesID = ifelse(speciesID == 12250,  # correcting the treatment of Western Flycatcher
                              12231,speciesID))


  saveRDS(status_tbl,"data_2026/SocbStatus.rds")


  status_tbl_ids <- unique(status_tbl$speciesID)

  if(any(!status_tbl_ids %in% species_groups$speciesID)){
    sp_miss <- status_tbl_ids[-which(status_tbl_ids %in% species_groups$speciesID)]
    sp_miss_names <- species_names[which(species_names$speciesID %in% sp_miss),"english_name"]
    warning(paste(length(sp_miss),"species with status are missing from species_groups",
                   paste(sp_miss_names$english_name,collapse = "; ")))
  }

  sp_group_ids <- unique(species_groups$speciesID)

  if(any(!sp_group_ids %in% status_tbl_ids)){
    sp_miss2 <- sp_group_ids[-which(sp_group_ids %in% status_tbl_ids)]
    sp_miss_names2 <- species_names[which(species_names$speciesID %in% sp_miss2),"english_name"]
    warning(paste(length(sp_miss2),"species in species_groups are missing from status table",
                  paste(sp_miss_names2$english_name,collapse = "; ")))
  }

}




species_groups <- readRDS("data_2026/species_groups.rds")
source("functions/specid_rename.R")# function that reconciles variations in nature counts column names


# removing one mistaken group - groups with the "All" suffix should be subgroupID 0

## manual re-setting of Western Flycatcher that was re-defined between the
## model runs and the publication of the report. This re-setting only affects
## the speciesID value, all data for the species are unaffected.
species_groups <- species_groups %>%
  filter(!(subgroupID != 0 & groupName == "Wetland Birds: All"))%>%
  mutate(speciesID = ifelse(speciesID == 12250,
                            12231,speciesID))


species_names <- readRDS("data_2026/species_names.rds")




sp_g <- species_groups %>%
  select(speciesID,groupName,include) %>%
  distinct() %>%
  pivot_wider(names_from = groupName,
              values_from = include) %>%
  left_join(species_names, by = "speciesID")

group_order <- c("Waterfowl",
            "Birds of Prey",
            "Wetland Birds",
            "Marine Birds",
            "Forest Birds",
            "Arctic Birds",
            "Long-Distance Migrants",
            "Shorebirds",
            "Aerial Insectivores",
            "Grassland Birds")

group_drop <- c("Edge/Early",
                "Urban",
                "Piscivores",
                "Galliformes",
                "Harvested",
                "Waterfowl: Ducks",
                "Dabblers",
                "Diving Ducks",
                "Treed Wetland Birds",
                "Marsh Birds",
                "Nesting in Canada",
                "Does not nest in Canada",
                "Forest Crop Birds",
                "Western",
                "Boreal Birds",
                "Carolinian Birds",
                "Arctic Tundra Birds",
                "Alpine Tundra Birds",
                "Long Distance Migrants",
                "Short Distance Migrants",
                "Arctic-breeders",
                "Boreal-breeders",
                "Coastal breeders",
                "Grassland breeders",
                "Sensitive to Linear Disturbance",
                "Shrub Steppe Birds")

is_y <- function(x){
  out <- x == "Y"
  out <- ifelse(any(out),
                TRUE,FALSE)
}
table_out <- sp_g %>%
  relocate(scientific_name,
           english_name,
           french_name,
           starts_with(group_order),
           speciesID) %>%
  select(-contains(group_drop)) %>%
  rename(NatureCounts_speciesID = speciesID)%>%
  rowwise() %>%
  mutate(included_in_composite_group = ifelse(is_y(c_across(4:31)),
                                              TRUE,FALSE))


survey_join <- data.frame(resultsCode = c("LINCOLN",
                                          "MIDWINTER",
                                          "CBC",
                                          "BBS",
                                          "WBPHS",
                                          "SCMP",
                                          "PRISM",
                                          "BCCWS",
                                          "No suitable data"),
                          Main_Survey = c("Lincoln Estimates of Population Size",
                                          "Midwinter Surveys",
                                          "Christmas Bird Count (CBC)",
                                          "North American Breeding Bird Survey (BBS)",
                                          "Waterfowl Breeding Population and Habitat Survey in Western Canada and the Northwestern United States",
                                          "Seabird Colony Monitoring Program",
                                          "Program for Regional and International Shorebird Monitoring (PRISM) - Fall Migration",
                                          "British Columbia Coastal Waterbird Survey",
                                          "No Suitable Data")
                          )
# manual re-setting of Red Knot PRISM data as full species ----------------
# manual step necessary because species account was modified between the
# model runs for the publication and the release of the report
# This manual re-set ensures that the species data are included and
# represent the published results.
rank_tbl <- readRDS("data_2026/SocbTrendRank.rds")%>%
  rename_with(.,.fn = specid_rename) %>%
  # mutate(popType = ifelse(speciesID == 4670 & subspeciesID == 47598 & resultsCode == "PRISM",
  #                         1,popType)) %>%
  select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
         resultsCode,popType,areaCode)%>%
  filter(goalTrend == "Y",
         popType == 1) %>%
  select(speciesID,
         resultsCode,areaCode)%>%
  distinct()

# rank_tbl <- rank_tbl %>%
#   select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
#          resultsCode,popType,areaCode)%>%
#   filter(goalTrend == "Y",
#          popType == 1) %>%
#   distinct()


  status_tbl <- readRDS("data_2026/SocbStatus.rds")
  status_goal_join <- status_tbl %>%
    filter(popType == 1) %>%
    select(speciesID,popStatusEn,objective,objectiveEn,
           goalEn,goalTrend) %>%
    distinct()




table_out_w_survey<- table_out %>%
  left_join(rank_tbl,
            by = c("NatureCounts_speciesID" = "speciesID"))%>%
  left_join(status_goal_join,
            by = c("NatureCounts_speciesID" = "speciesID"))%>%
  mutate(resultsCode = ifelse(is.na(resultsCode),
                              "No suitable data",
                              resultsCode)) %>%
  left_join(survey_join,
            by = "resultsCode") %>%
  select(-c(resultsCode,taxon_group)) %>%
  relocate(scientific_name,english_name,french_name,
           included_in_composite_group,
           Main_Survey,
           popStatusEn,
           objectiveEn,
           goalEn,
           goalTrend) %>%
  rename(Scientific_Name = scientific_name,
         English_Name = english_name,
         French_Name = french_name)



miss_sp <- table_out_w_survey %>% filter(is.na(Main_Survey))

write_excel_csv(table_out_w_survey,"Wide_format_SOCB_species_table2_2026.csv")



