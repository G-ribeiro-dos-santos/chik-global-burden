library(dplyr)
library(reactable)

rm(list=ls())

chik_cfr = read.csv("0_Data/cfr_measures.csv") #cfr per 1000
cfr = data.frame(age_min=seq(0,80,by=10), age_max=c(seq(9,79,by=10),110),cfr = chik_cfr |> group_by(Age) |> mutate(cfr=mean(CFR)/1000) |> ungroup() |> select(cfr) |> unique())
load("2_Burden/df_params.Rdata")
source("2_Burden/util_functions.R")

list_files = list.files("3_Vaccine_impact/simulation_by_country/", pattern="*.Rdata")
load("0_Data/demog_data_df.Rdata")

list_files = data.frame(file=list_files,
                        country=stringr::str_match(list_files, "_\\s*(.*?)\\s*.Rdata")[,2])

#Get continents
list_files$continent = countrycode::countrycode(sourcevar = list_files[, "country"],
                                                origin = "country.name",
                                                destination = "continent")

list_files$continent[list_files$country=="Kosovo"] = "Europe"
list_files$continent[list_files$country=="United_States_of_America"] = "Americas"
list_files$continent[list_files$country=="Saint_Lucia"] = "Americas"
list_files$continent[list_files$country=="Micronesia"] = "Oceania"

#Get CHIK status classification
list_files$classification = all_countries$status[sapply(list_files$country, function(x) which(all_countries$country==x))]
list_files$who_region = all_countries$who_region[sapply(list_files$country, function(x) which(all_countries$country==x))]

df_sims = lapply(list_files$file, function(file){

  classification = list_files$classification[list_files$file==file]
  country = list_files$country[list_files$file==file]
  continent = list_files$continent[list_files$file==file]
  who_region = list_files$who_region[list_files$file==file]
  print(country)
  pop_size = sum(all_countries[all_countries$country==country,c(as.character(1:99),"100+")])
  
  load(paste0("3_Vaccine_impact/simulation_by_country/",file))
  
  if(classification=="Endemic"){
    
    df_sims = df %>%
      mutate(bs_ID=rep(1:length(sims), nrow(df)/length(sims)))
    
  }else if(classification=="Epidemic"){
    
    df_sims =df %>%
      mutate(bs_ID=rep(1:length(sims),nrow(df_params)))
    
  }else{#No transmission
    
    df_sims = df %>%
      mutate(bs_ID=rep(1:length(sims),nrow(df_params)))
    
    pop_size=0
  }
  
  df_sims %>% mutate(country=country, who_region=who_region, continent=continent, classification=classification, pop_size=pop_size)
  
}) %>% do.call(what = rbind)


df_sims_light = df_sims %>% filter(grepl("scenario", Scenario) | Scenario %in% c("BC_cepi", "BC_gavi"))

save(df_sims, file="3_Vaccine_impact/simulations_sensitivity.Rdata")
save(df_sims_light, file="3_Vaccine_impact/simulations_basecase.Rdata")

