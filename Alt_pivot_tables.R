

# script to pivot from long version of regional information to wide form
# see email chain from Marcel on October 23, 2026

library(tidyverse)
library(readxl)

national <- read_excel("data/SOCBNat&RegionRanksV4.0_Oct2026.xlsx",
                       sheet = 1) |>
  group_by(english_name) |>
  filter(version == max(version)) |> # filter to the most recent version for a each species
  #filter(version == 2025) |>
  select(english_name,sar_rank,sar_designation,threat_score,final_score,rank_national) |>
  distinct()


regional <- read_excel("data/SOCBNat&RegionRanksV4.0_Oct2026.xlsx",
                       sheet = 2) |>
  filter(ebird_season %in% c("breeding","resident")) |>
  mutate(BCR_name = paste0(BCR,"_",descr)) |>
  #select(english_name,sar_rank,sar_designation,rank_regional,BCR_name) |>
  distinct() |>
  arrange(BCR_name)

region_sel <- regional |>
  pivot_wider(id_cols = c(english_name,sar_designation,sar_rank),
              names_from = BCR_name,
              values_from = rank_regional)

dups <- regional |>
  dplyr::summarise(n = dplyr::n(), .by = c(english_name, sar_designation, sar_rank, BCR_name)) |>
  dplyr::filter(n > 1L)



full <- national |> left_join(region_sel,
          by = c("english_name","sar_rank","sar_designation"))

write_csv(full,"draft_wide_format_scores.csv")
