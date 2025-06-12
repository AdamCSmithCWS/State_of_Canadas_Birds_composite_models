#comparing species lists

library(tidyverse)

species_groups <- readRDS("data/species_groups.rds")
source("functions/specid_rename.R")# function that reconciles variations in nature counts column names


# removing one mistaken group - groups with the "All" suffix should be subgroupID 0

## manual re-setting of Western Flycatcher that was re-defined between the
## model runs and the publication of the report. This re-setting only affects
## the speciesID value, all data for the species are unaffected.
species_groups <- species_groups %>%
  filter(!(subgroupID != 0 & groupName == "Wetland Birds: All"))%>%
  mutate(speciesID = ifelse(speciesID == 12250,
                            12231,speciesID))


species_names <- readRDS("data/species_names.rds") %>%
  mutate(speciesID = ifelse(speciesID == 12250,
                            12231,speciesID))




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

table_out <- sp_g %>%
  relocate(scientific_name,
           english_name,
           french_name,
           starts_with(group_order),
           speciesID) %>%
  select(-contains(group_drop)) %>%
  rename(NatureCounts_speciesID = speciesID)

# manual re-setting of Red Knot PRISM data as full species ----------------
# manual step necessary because species account was modified between the
# model runs for the publication and the release of the report
# This manual re-set ensures that the species data are included and
# represent the published results.
rank_tbl <- readRDS("data/SocbTrendRank.rds")%>%
  rename_with(.,.fn = specid_rename) %>%
  mutate(popType = ifelse(speciesID == 4670 & subspeciesID == 47598 & resultsCode == "PRISM",
                          1,popType)) %>%
  select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
         resultsCode,popType,areaCode)%>%
  filter(goalTrend == "Y",
         popType == 1) %>%
  select(speciesID,
         resultsCode,areaCode)%>%
  distinct()






table_out_w_survey<- table_out %>%
  left_join(rank_tbl,
            by = c("NatureCounts_speciesID" = "speciesID")) %>%
  relocate(scientific_name,english_name,french_name,
           resultsCode) %>%
  mutate(resultsCode = ifelse(is.na(resultsCode),
                              "No suitable data",
                              resultsCode)) %>%
  select(-taxon_group)



miss_sp <- table_out_w_survey %>% filter(is.na(resultsCode))

write_excel_csv(table_out_w_survey,"Wide_format_SOCB_species_table.csv")



