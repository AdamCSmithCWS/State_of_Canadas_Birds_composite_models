### fitting the composite group indicator models for the State of Canada's Birds
###
library(tidyverse)
library(cmdstanr)
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


re_download <- FALSE # set to TRUE if re-downloading the data files from NatureCounts database
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
  saveRDS(rank_tbl,"data/SocbTrendRank.rds")

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbSpecies",
                 timeout = 120) -> sp_tbl
  saveRDS(sp_tbl,"data/SocbSpecies.rds")

  #Information on species' group#
  nc_query_table(username = "adam.smith",
                 table = "Groups",
                 timeout = 120) -> group_tbl
  saveRDS(group_tbl,"data/Groups.rds")

  #Data on all annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndices",
                 timeout = 240) -> indices_tbl

  saveRDS(indices_tbl,"data/TrendsIndices.rds")

  #Data on goal-based annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndicesGoals",
                 timeout = 120) -> goal_indices_tbl
  saveRDS(goal_indices_tbl,"data/TrendsIndicesGoals.rds")

  #Data on species trends #
  nc_query_table(username = "adam.smith",
                 table = "Trends",
                 timeout = 120) -> trends_tbl
  saveRDS(trends_tbl,"data/Trends.rds")


  #Species groups
  #
  species_groups <- nc_query_table(table = "SocbTrendGroups",
                                   username = "adam.smith",
                                   verbose = TRUE) %>%
    rename_with(., .fn = specid_rename) %>%
    mutate(groupName = gsub("\r\n",
                            "",
                            x = groupName,
                            fixed = TRUE),
           groupNameFr = gsub("\r\n",
                              "",
                              x = groupNameFr,
                              fixed = TRUE),
           include = toupper(include)) %>%
    filter(!(speciesID == 12231 & populationID == 8414)) # adjusting Western Flycatcher for recent change in taxonomy


  saveRDS(species_groups,"data/species_groups.rds")


  species_names <- naturecounts::search_species() %>%
    rename_with(.,.fn = specid_rename) %>%
    filter(speciesID %in% species_groups$speciesID)

  saveRDS(species_names,"data/species_names.rds")

}



# load previously downloaded datafiles ------------------------------------

#species population table including three generation length,
# canadian occurrence, population goals, etc.
sp_tbl <- readRDS("data/SocbSpecies.rds") %>%
  rename_with(.,.fn = specid_rename)

# surveys assocaited with species goals
rank_tbl <- readRDS("data/SocbTrendRank.rds") %>%
  rename_with(.,.fn = specid_rename)

# annual indices of abundance assocaited with species goal
goal_indices_tbl <- readRDS("data/TrendsIndicesGoals.rds") %>%
  rename_with(.,.fn = specid_rename)

# species inclusion table for composite groups
species_groups <- readRDS("data/species_groups.rds")%>%
  mutate(speciesID = ifelse(speciesID == 12250,  # correcting the treatment of Western Flycatcher
                            12231,speciesID))

# additional species names table
species_names <- readRDS("data/species_names.rds") %>%
  mutate(speciesID = ifelse(speciesID == 12250,  # correcting the treatment of Western Flycatcher
                            12231,speciesID))

# Prepare data for smooth models ------------------------------------------

base_year <- 1970
re_smooth <- FALSE

# id the best survey datasets for each species ----------------------------
rank_tbl <- rank_tbl %>%
  select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
         resultsCode,popType,areaCode)%>%
  filter(goalTrend == "Y",
         popType == 1) %>%
  distinct()

tst<- rank_tbl %>%
  anti_join(species_names)
if(nrow(tst) > 0){stop(paste(nrow(tst),"species are missing from at least one dataset"))}

sp_simple <- sp_tbl %>%
  select(speciesCode,speciesID,
         population,
         objective,
         popID,
         popType,
         confidence) %>% #,population,popType,popID) %>%
  filter(popType == 1) %>%
  left_join(.,species_names) %>%
  distinct()


traj_sel <- rank_tbl %>%
  inner_join(.,sp_simple)



goal_indices_tbl <- goal_indices_tbl %>%
  select(resultsCode,speciesID,
         year, index, indexUpperCI,indexLowerCI,
         areaCode)%>%
  distinct() %>%
  inner_join(.,traj_sel) %>%
  filter(year >= base_year) %>%
  distinct()


re_summarise_indices <- TRUE

if(re_summarise_indices){
all_inds <- goal_indices_tbl %>%
  inner_join(.,traj_sel) %>%
  filter(year >= base_year) %>%
  distinct()

# identify a few rows where the index of lower CI of the indices are 0 and
# therefore will throw NA values with a log transformation (Ross's Goose Lincoln estimates)
miss_inds <- all_inds %>%
  filter(index <= 0 |
           indexLowerCI <= 0 |
           is.infinite(indexUpperCI))

# fix the lower CI for these 0 values so they can be log-transformed
if(nrow(miss_inds) > 0){
  miss_sp <- unique(miss_inds$speciesCode)
  all_inds <- all_inds %>%
    rowwise() %>%
    mutate(indexLowerCI = ifelse(indexLowerCI == 0,max(0.1,index-indexUpperCI),indexLowerCI)) #fixes lower bound for a handful of
  warning(paste((paste(miss_sp,collapse = ", ")),"had missing index or CI information. Zeros have been replaced with arbitrarily small values"))
}



all_inds <- all_inds %>%
  filter(index > 0)



all_inds <- all_inds %>%
  mutate(ln_index = log(index),
         ln_lci = log(indexLowerCI),
         ln_uci = log(indexUpperCI),
         ln_index_sd = (ln_uci-ln_lci)/3.9) %>%
  group_by(speciesID) %>%
  filter(speciesCode != "norcro") %>%
  mutate(yearn = year-(min(year)-1)) %>% # sets a yearn value specific to each species
  arrange(speciesID,year)

saveRDS(all_inds,"data/all_socb_goal_indices.rds")

}else{

all_inds <- readRDS("data/all_socb_goal_indices.rds")

}


all_sp <- unique(all_inds$speciesID)

all_smoothed_indices <- NULL



# Generate species-level smooths ------------------------------------------


if(re_smooth){
  # Generate species-level smooths ------------------------------------------

  for(species in all_sp){

    inds <- all_inds %>%
      filter(speciesID == species) %>%
      select(speciesID,year,
             ln_index,ln_index_sd,yearn) %>%
      distinct()

    years <- c(min(inds$year):max(inds$year))
    n_years <- as.integer(length(years))
    n_indices <- as.integer(nrow(inds))
    n_knots <- as.integer(round(n_indices/3))
    start_year <- min(inds$year)

    gam_data <- gam_basis(years,
                          nknots = n_knots,
                          sm_name = "year")
    stan_data <- list(
      n_years = n_years,
      n_indices = n_indices,
      n_knots_year = gam_data$nknots_year,
      year = inds$yearn,
      ln_index = inds$ln_index,
      ln_index_sd = inds$ln_index_sd,
      year_basis = gam_data$year_basis
    )
    ## fit model with cmdstanr
    file <- "models/GAM_smooth_model.stan"
    mod <- cmdstan_model(file)

    fit <- mod$sample(data = stan_data,
                      parallel_chains = 4,
                      refresh = 0,
                      adapt_delta = 0.95,
                      show_messages = FALSE,
                      show_exceptions = FALSE)

    sum <- fit$summary(variables = NULL,
                       "mean",
                       "sd",
                       "ess_bulk",
                       "rhat",
                       q2_5 = q2_5,
                       q97_5 = q97_5)

    mx_rhat <- max(sum$rhat,na.rm = TRUE)
    if(mx_rhat > 1.05){
      fit <- mod$sample(data = stan_data,
                        parallel_chains = 4,
                        iter_warmup = 8000,
                        iter_sampling = 8000,
                        thin = 8,
                        refresh = 0,
                        show_exceptions = FALSE)
      sum <- fit$summary(variables = NULL,
                         "mean",
                         "sd",
                         "ess_bulk",
                         "rhat",
                         q2_5 = q2_5,
                         q97_5 = q97_5)
    }
    smooth_inds <- sum %>%
      filter(grepl("mu_smooth",variable)) %>%
      mutate(yearn = as.integer(str_extract_all(variable,
                                                "[[:digit:]]{1,}",
                                                simplify = TRUE)),
             speciesID = species,
             smooth_ind = mean,
             smooth_ind_sd = sd,
             smooth_ind_lci = q2_5,
             smooth_ind_uci = q97_5) %>%
      select(speciesID,yearn,smooth_ind,smooth_ind_sd,
             smooth_ind_lci,smooth_ind_uci)

    diffs <- sum %>%
      filter(grepl("annual_diffs",variable)) %>%
      mutate(yearn = as.integer(str_extract_all(variable,
                                                "[[:digit:]]{1,}",
                                                simplify = TRUE)),
             annual_diff = mean,
             annual_diff_sd = sd) %>%
      select(yearn,annual_diff,annual_diff_sd)

    scaled_status <- sum %>%
      filter(grepl("scaled_status",variable)) %>%
      mutate(yearn = as.integer(str_extract_all(variable,
                                                "[[:digit:]]{1,}",
                                                simplify = TRUE)),
             scaled_status = mean,
             scaled_status_sd = sd) %>%
      select(yearn,scaled_status,scaled_status_sd)

    scaled_log_status <- sum %>%
      filter(grepl("scaled_log_status",variable)) %>%
      mutate(yearn = as.integer(str_extract_all(variable,
                                                "[[:digit:]]{1,}",
                                                simplify = TRUE)),
             scaled_log_status = (mean),
             scaled_log_status_sd = sd) %>%
      select(yearn,scaled_log_status,scaled_log_status_sd)


    smooth_inds <- smooth_inds %>%
      inner_join(.,diffs,
                 by = c("yearn"))%>%
      inner_join(.,scaled_status,
                 by = c("yearn"))%>%
      inner_join(.,scaled_log_status,
                 by = c("yearn")) %>%
      mutate(year = yearn+(start_year-1))

    all_smoothed_indices <- bind_rows(all_smoothed_indices,smooth_inds)



    print(paste(species,"complete",round(which(all_sp == species)/length(all_sp),2)))

  }

  saveRDS(all_smoothed_indices,"data/socb_smoothed_indices.rds")

}

