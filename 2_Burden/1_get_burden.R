library(dplyr)

rm(list=ls())

# This script was written to run on a multi-core CPU from an HPC

source("2_Burden/util_functions.R")
load("0_Data/demog_data_df.Rdata")
load("2_Burden/df_params.Rdata")
chik_cfr = read.csv("0_Data/cfr_measures.csv") #cfr per 1000
cfr = data.frame(age_min=seq(0,80,by=10), age_max=c(seq(9,79,by=10),110),cfr = chik_cfr |> group_by(Age) |> mutate(cfr=mean(CFR)/1000) |> ungroup() |> select(cfr) |> unique())
rm(chik_cfr)

#Overall parameters
district_size = 10e6
bootstrap = 1e4
num_cores = parallel::detectCores()

time_start = Sys.time()

df_burden = parallel::mclapply(which(all_countries$status!="NoTransmission"), mc.cores=num_cores, function(task_id){

  #Get country params
  country = all_countries$country[task_id]
  pop_size_eff = (all_countries$pop_size*all_countries$prop_eff_pop)[task_id]
  n_district = max(round(pop_size_eff/district_size), 1)
  s0 = all_countries$s0[task_id]
  demog_vect_eff = all_countries[task_id, c(as.character(1:99),"100+")]*all_countries$prop_eff_pop[task_id]
  demog_vect_district = demog_vect_eff/sum(demog_vect_eff)*district_size
  life_expectancy =all_countries$life_expectancy[task_id]
  classification = all_countries$status[task_id]
  who_region = all_countries$who_region[task_id]
  reliable = all_countries$reliable_healthcare[task_id]
  
  print(country)
  
  #Get transmission dynamics parameters
  if(classification=="Epidemic"){
    Var = df_params$Var_outbreak[1]
    Lambda = df_params$Lambda_outbreak[1]
    Mu = df_params$Mu_outbreak[1]
  }
  
  if(classification=="Endemic"){
    Var = df_params$Var_outbreak_endemic[1]
    Lambda = df_params$Lambda_outbreak_endemic[1]
    Mu = 1
  }
  
  if(classification=="NoTransmission"){
    Var = df_params$Var_outbreak[1]
    Lambda = df_params$Lambda_outbreak[1]
    Mu = -1
  }
  
  #Get fois
  alpha = ((1-Lambda)/Var-(1/Lambda))*Lambda^2
  beta = alpha*(1/Lambda-1)
  national_foi = rbeta(bootstrap, alpha, beta)
  
  if(reliable){
    Var = (all_countries$case_foi_sd[task_id])^2
    Lambda = all_countries$case_foi[task_id]
    Mu = df_params$Mu_outbreak[1]
    national_foi = rnorm(n = bootstrap, mean = Lambda, sd = sqrt(Var))
    national_foi[national_foi<0] = 0
  }
  
  #Get burden estimates
  df_burden = lapply(1:bootstrap, function(bs){

    n_outbreak = sum(runif(n_district,0,1)<=Mu)
    
    infections = national_foi[bs]*demog_vect_district*s0
    cases = infections*df_params$prob_symptoms[1]
    deaths = unlist(sapply(1:length(cases), function(a) cases[a]*cfr$cfr[which(cfr$age_min<=a & cfr$age_max>=a)]))
    ylls = get_yll(age_vect = (1:length(cases))-0.5,
                   n_infections_vect = infections,
                   parms = df_params[1,],
                   regular_life_expectancy = life_expectancy,
                   cfr_df = cfr, infected_to_case = df_params$prob_symptoms[1])
    ylds = get_yld(age_vect = (1:length(cases))-0.5,
                   n_infections_vect = infections,
                   parms = df_params[1,],
                   regular_life_expectancy = life_expectancy,
                   cfr_df = cfr, infected_to_case = df_params$prob_symptoms[1])
    dalys = ylls+ylds
    national_burden = data.frame(infections = sum(infections)*n_outbreak,
                                 cases = sum(cases)*n_outbreak,
                                 deaths= sum(deaths)*n_outbreak,
                                 ylls = sum(ylls)*n_outbreak,
                                 ylds = sum(ylds)*n_outbreak,
                                 dalys = sum(dalys)*n_outbreak,
                                 infections_u5 = sum(infections[1:5])*n_outbreak,
                                 cases_u5 = sum(cases[1:5])*n_outbreak,
                                 deaths_u5 = sum(deaths[1:5])*n_outbreak,
                                 ylls_u5 = sum(ylls[1:5])*n_outbreak,
                                 ylds_u5 = sum(ylds[1:5])*n_outbreak,
                                 dalys_u5 = sum(dalys[1:5])*n_outbreak)
    
    national_burden = national_burden*pop_size_eff/(n_district*district_size)
    return(national_burden)
    
  }) %>% do.call(what=rbind)
  
  df_burden = as.data.frame(df_burden) %>% 
    mutate(country=country, bs_ID=1:bootstrap) %>%
    relocate(country, bs_ID)
  
}) %>% do.call(what=rbind)
print(Sys.time()-time_start)

#No transmission countries
burden_notrans = data.frame(country=all_countries$country[which(all_countries$status=="NoTransmission")],
                            bs_ID = 1,
                            infections = 0,
                            cases = 0,
                            deaths= 0,
                            ylls = 0,
                            ylds = 0,
                            dalys = 0,
                            infections_u5 = 0,
                            cases_u5 = 0,
                            deaths_u5 = 0,
                            ylls_u5 = 0,
                            ylds_u5 = 0,
                            dalys_u5 = 0)

#Collate and save
df_burden = rbind(df_burden, burden_notrans)

write.csv(df_burden, "2_Burden/burden_estimates.csv", row.names = F)





