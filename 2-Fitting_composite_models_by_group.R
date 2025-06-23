### fitting the composite group indicator models for the State of Canada's Birds
###
library(tidyverse)
library(cmdstanr)
library(readxl)
#library(naturecounts)
library(patchwork)

source("functions/GAM_basis_function_mgcv.R")

# Quantile functions to generate 95% CI instead of stan's default 90
q2_5 <- function(x)c(q2_5 = quantile(x,probs = c(0.025),
                                     names = FALSE))
q97_5 <- function(x)c(q97_5 = quantile(x,probs = c(0.975),
                                       names = FALSE))

q5 <- function(x)c(q5 = quantile(x,probs = c(0.05),
                                     names = FALSE))
q95 <- function(x)c(q95 = quantile(x,probs = c(0.95),
                                       names = FALSE))

q10 <- function(x)c(q10 = quantile(x,probs = c(0.1),
                                     names = FALSE))
q90 <- function(x)c(q90 = quantile(x,probs = c(0.9),
                                       names = FALSE))


q25 <- function(x)c(q25 = quantile(x,probs = c(0.25),
                                   names = FALSE))
q75 <- function(x)c(q75 = quantile(x,probs = c(0.75),
                                   names = FALSE))

q20 <- function(x)c(q20 = quantile(x,probs = c(0.20),
                                   names = FALSE))
q80 <- function(x)c(q80 = quantile(x,probs = c(0.80),
                                   names = FALSE))



base_year <- 1970

species_groups <- readRDS("data/species_groups.rds")

# group_names <- unique(species_groups$groupName)
#
# groups_published <- group_names[c(1,3,6,7,
#                                   23,24,25,
#                                   18,21,22,
#                                   8,
#                                   26,28,27,
#                                   41,45,46,44,
#                                   55,56,57,
#                                   11,
#                                   38,39,40,
#                                   33,34,35)]
#
# saveRDS(groups_published,"data/published_groups.rds")

groups_published <- readRDS("data/published_groups.rds")

species_names <- readRDS("data/species_names.rds")

all_smoothed_indices <- readRDS("data/socb_smoothed_indices.rds")


all_inds <- readRDS("data/all_socb_goal_indices.rds") %>%
  select(-c(english_name,french_name,taxon_group,scientific_name,yearn))%>%
  filter(!speciesID == 2840) #dropping black vulture which was removed from all indicators


species_confidence <- all_inds %>%
  select(speciesID,confidence) %>%
  distinct()


all_smoothed_indices <- all_smoothed_indices %>%
  left_join(.,species_names,
            by = "speciesID")


all_ind_compare <- all_smoothed_indices %>%
  full_join(., all_inds)



# Group loop --------------------------------------------------------------

drop_low_confidence <- TRUE
low_confid <- c("DD")
plot_low_confidence <- FALSE

annual_status_combine <- NULL

if(!drop_low_confidence){

  sp_drop_low_confidence <- species_confidence %>%
    filter(confidence %in% low_confid) %>%
    select(speciesID) %>%
    unlist

  all_smoothed_indices <- all_smoothed_indices %>%
  filter(!speciesID %in% sp_drop_low_confidence)
}


re_fit_all <- TRUE
groups_to_refit <- groups_published[c(12,23,25)] #groups_published

for(grp in groups_published){



    grp_labl <- species_groups %>%
      filter(groupName == grp) %>%
      select(groupName,groupNameFr) %>%
      distinct() %>%
      unlist()
    grp_labl <- gsub("[[:punct:]]","",x = grp_labl)
    grp_labl <- gsub("[[:blank:]]","_",x = grp_labl)
    grp_labl <- paste(grp_labl,collapse = "-")



    species_sel <- species_groups %>%
    filter(groupName == grp,
           include == "Y") %>%
    distinct()



  # Dropping low confidence species -----------------------------------------


  if(drop_low_confidence){
sp_drop_low_confidence <- species_confidence %>%
  filter(confidence %in% low_confid) %>%
  select(speciesID) %>%
  unlist


  }else{
    sp_drop_low_confidence <- NULL
  }

  if(nrow(species_sel) == 0){next}
  if(nrow(species_sel) != length(unique(species_sel$speciesID))){
    stop(paste("Problem with species list for",grp))
  }
  n_sp_w_data <- length(which(species_sel$speciesID %in% all_smoothed_indices$speciesID))


  species_sel <- species_sel %>%
    mutate(sufficient_data = ifelse(speciesID %in% all_smoothed_indices$speciesID,
                                    TRUE,
                                    FALSE),
           sufficient_confidence = ifelse(speciesID %in% sp_drop_low_confidence,
                                   FALSE,
                                   TRUE)) %>%
    left_join(.,species_names,
              by = "speciesID") %>%
    select(-taxon_group)

  n_sp_w_low_conf <- length(which(species_sel$speciesID %in% sp_drop_low_confidence))

  print(paste("There are data for",(n_sp_w_data-n_sp_w_low_conf),"of",nrow(species_sel),
              "in the",grp,"group"))

if((n_sp_w_data-n_sp_w_low_conf)/nrow(species_sel) < 0.2){
  print(paste("Skipping",grp,"because only",round((n_sp_w_data-n_sp_w_low_conf)/nrow(species_sel),2)*100,"% of species have data"))
    next}

  ### Drop the base-year values and other species
inds_all <- all_smoothed_indices %>%
  filter(speciesID %in% species_sel$speciesID,
         !speciesID %in% sp_drop_low_confidence) %>%
  mutate(species_ind = as.integer(factor(speciesID)))#

years_w_GT_50 <- inds_all %>%
  group_by(year) %>%
  summarise(n_sp = n(),.groups = "drop") %>%
  mutate(p_sp = n_sp/max(n_sp,na.rm = TRUE)) %>%
  filter(p_sp > 0.5)
min_yr <- range(years_w_GT_50$year)[1]
max_yr <- range(years_w_GT_50$year)[2]


base_yr <- max(base_year,min_yr)

inds <- inds_all %>%
  group_by(speciesID,species_ind) %>%
  mutate(yearn2 = year-base_yr) %>%
  filter(year <= max_yr,
         year > base_yr,
         scaled_log_status_sd > 0,
         scaled_status_sd > 0,
         annual_diff_sd > 0) %>%
  arrange(species_ind,yearn2)


## track the start and end years for each species
sp_y <- inds %>%
  group_by(species_ind,speciesID) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            first_yearn2 = min(yearn2),
            last_yearn2 = max(yearn2),
            .groups = "drop")

inds_w_low_conf <- all_smoothed_indices %>%
  filter(speciesID %in% species_sel$speciesID) %>%
  group_by(speciesID) %>%
  mutate(yearn2 = year-(base_yr-1)) %>%
  filter(year <= max_yr,
         year >= base_yr,
         scaled_log_status_sd > 0,
         scaled_status_sd > 0,
         annual_diff_sd > 0) %>%
  arrange(speciesID,yearn2)

sp_y_low_conf <- inds_w_low_conf %>%
  group_by(speciesID) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            .groups = "drop")

# number of years and species
n_years <- max(inds$yearn2,na.rm = TRUE)
n_species <- max(inds$species_ind,na.rm = TRUE)


species <- matrix(data = NA,
                  nrow = n_years,
                  ncol = n_species)
ln_index <- matrix(data = 0,
                   nrow = n_years,
                   ncol = n_species)
ln_index_sd <- matrix(data = 0,
                      nrow = n_years,
                      ncol = n_species)
n_species_year <- vector("integer",length = n_years)

re_fit <- ifelse(re_fit_all | grp %in% groups_to_refit,
                 TRUE,FALSE)

if(re_fit){

for(y in 1:n_years){
  tmp <- inds %>%
    filter(yearn2 == y) %>%
    arrange(species_ind)
  n_sp <- nrow(tmp)
  sp_inc <- as.integer(tmp$species_ind)
  species[y,1:n_sp] <- sp_inc
  if(n_sp < n_species){
    sp_miss <- c(1:n_species)[-sp_inc]
    species[y,c((n_sp+1):n_species)] <- sp_miss
  }
  n_species_year[y] <- n_sp
  for(s in sp_inc){
    ln_index[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"annual_diff"])
    ln_index_sd[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"annual_diff_sd"])
  }
}

stan_data2 <- list(n_years = n_years,
                   n_species = n_species,
                   n_species_year = n_species_year,
                   species = species,
                   ln_index = ln_index,
                   ln_index_sd = ln_index_sd)


file2 <- "models/State_of_Birds_Model_differences.stan"
mod2 <- cmdstan_model(file2)

fit2 <- mod2$sample(data = stan_data2,
                    parallel_chains = 4,
                    refresh = 0,
                    adapt_delta = 0.9,
                    iter_warmup = 2000,
                    iter_sampling = 4000,
                    show_exceptions = FALSE,
                    seed = 2024)

fit2$save_object(paste0("output/stan_fit_",grp_labl,".rds"))



sum2 <- fit2$summary(variables = NULL,
                     "mean",
                     "sd",
                     "ess_bulk",
                     "rhat",
                     q2_5 = q2_5,
                     q5 = q5,
                     q10 = q10,
                     q20 = q20,
                     q25 = q25,
                     q75 = q75,
                     q80 = q80,
                     q90 = q90,
                     q95 = q95,
                     q97_5 = q97_5)

mx_rhat2 <- max(sum2$rhat,na.rm = TRUE)



if(mx_rhat2 > 1.02){
  fit2 <- mod2$sample(data = stan_data2,
                     parallel_chains = 4,
                     iter_warmup = 4000,
                     iter_sampling = 4000,
                     thin = 2,
                     refresh = 0,
                     adapt_delta = 0.9,
                     show_exceptions = FALSE)
  sum2 <- fit2$summary(variables = NULL,
                       "mean",
                       "sd",
                       "ess_bulk",
                       "rhat",
                       q2_5 = q2_5,
                       q5 = q5,
                       q10 = q10,
                       q20 = q20,
                       q25 = q25,
                       q75 = q75,
                       q80 = q80,
                       q90 = q90,
                       q95 = q95,
                       q97_5 = q97_5)
  mx_rhat2 <- max(sum2$rhat,na.rm = TRUE)
}

annual_status_difference <- sum2 %>%
  filter(grepl("annual_status",variable,fixed = TRUE)) %>%
  mutate(yearn2 = as.integer(str_extract_all(variable,
                                             "[[:digit:]]{1,}",
                                             simplify = TRUE)),
         model = "difference",
         year = yearn2+(base_yr-1),
         percent_diff = (exp(mean)-1)*100,
         percent_diff_lci = (exp(q2_5)-1)*100,
         percent_diff_uci = (exp(q97_5)-1)*100)


p_gt <- function(x,threshold = 0){
  p <- length(which(x > threshold))/length(x)
  return(p)
}

p_lt <- function(x,threshold = 0){
  p <- length(which(x < threshold))/length(x)
  return(p)
}

p_decrease <- fit2$draws(variables = "annual_status",
                         format = "df") %>%
  pivot_longer(cols = starts_with("annual_status"),
               values_to = "annual_status",
               names_prefix = "annual_status\\[",
               names_to = "yearn2") %>%
  mutate(yearn2 = as.integer(gsub("\\]","",yearn2))) %>%
  group_by(yearn2) %>%
  summarise(p_decrease = p_lt(annual_status),
            p_increase = p_gt(annual_status))

annual_status_difference <- annual_status_difference %>%
  left_join(p_decrease, by = "yearn2")

sigma2 <- sum2 %>%
  filter(grepl("sigma",variable))

annual_status_difference <- annual_status_difference %>%
  mutate(groupName = grp,
         model = "difference")
inds_all <- inds_all %>%
  mutate(groupName = grp)
species_sel <- species_sel %>%
  mutate(groupName = grp)

  saveRDS(annual_status_difference,paste0("output/composite_fit_",grp_labl,".rds"))
  saveRDS(inds_all,paste0("output/composite_data_",grp_labl,".rds"))
  saveRDS(inds_w_low_conf,paste0("output/composite_data_w_low_confid_",grp_labl,".rds"))
  saveRDS(species_sel,paste0("output/composite_species_list_",grp_labl,".rds"))

    }else{#end if re-fit
      annual_status_difference <- readRDS(paste0("output/composite_fit_",grp_labl,".rds"))
      inds_all <- readRDS(paste0("output/composite_data_",grp_labl,".rds"))
      species_sel <- readRDS(paste0("output/composite_species_list_",grp_labl,".rds"))
      inds_w_low_conf <- readRDS(paste0("output/composite_data_w_low_confid_",grp_labl,".rds"))


    }
annual_status_combine <- bind_rows(annual_status_combine,annual_status_difference)

saveRDS(annual_status_combine,"output/annual_status_combine.rds")

}


