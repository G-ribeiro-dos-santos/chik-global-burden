#Doc on seed setting for parallel simulations ; https://irudnyts.github.io/setting-a-seed-in-r-when-using-parallel-simulation/

# This script was written to run on a multi-core CPU from an HPC - task_id is an argument passed on from the bash script being the array ID from the job array (from 1 to 190, one run by country)

library(parallel)
library(deSolve)
library(ggplot2)
library(dplyr)

rm(list=ls())

time_start = Sys.time()

RNGkind("L'Ecuyer-CMRG")
n_cores = detectCores()

#Arguments passed
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

load("0_Data/demog_data_df.Rdata")
chik_cfr = read.csv("0_Data/cfr_measures.csv") #cfr per 1000
cfr = data.frame(age_min=seq(0,80,by=10), age_max=c(seq(9,79,by=10),110),cfr = chik_cfr |> group_by(Age) |> mutate(cfr=mean(CFR)/1000) |> ungroup() |> select(cfr) |> unique())

source("2_Burden/util_functions.R")

rm(df_design, chik_cfr)

bootstrap = 1000
Time_vect = 1:20 #years

#Age structure, life expectancy and S0
demog_vect = all_countries[task_id, c(as.character(1:99),"100+")]
life_expectancy =all_countries$life_expectancy[task_id]
prop_no_vacc_vect = unlist((c(rep(0,14),all_countries[task_id, paste0("Births_",15:49)],rep(0,51)))*demog_vect + #Pregnant
                             c(rep(1.6,39),rep(2.3,10),rep(4.4,10),rep(3.9,10),rep(3.1,10),rep(2.5,21))/100*demog_vect)/sum(unlist(demog_vect)) #Immunocompromised

S_0_init_who = all_countries$s0[task_id]

Age_structure = unlist((demog_vect/sum(demog_vect)))

flag_infection = rep(F, bootstrap)
time_infection = rep(-1e3, bootstrap)

foi_dist = "beta"

#Get parameters table
load("2_Burden/df_params.Rdata")

#Get long term abacus
load("3_Vaccine_impact/abacus_long_term.Rdata")
abacus_long_term = df_abacus %>% filter(country==all_countries$country[task_id])

#Run simulations
#Epidemic----------------
if(all_countries$status[task_id]=="Epidemic"){
  
  if(all_countries$reliable_healthcare[task_id]){
    df_params$Lambda_outbreak = all_countries$case_foi[task_id]
    df_params$Var_outbreak = (all_countries$case_foi_sd[task_id])^2
    foi_dist="normal"
  }
  
  get_beta_longterm_epi = stats::approxfun(abacus_long_term$foi, y = abacus_long_term$beta0, method="linear", rule = 2)
  beta_vect = NULL
  
  #Vaccine impact
  sims=c()
  for(row_vec in list(which(grepl("cepi",rownames(df_params)) & !grepl("dyn",rownames(df_params))),
                      which(grepl("cepi",rownames(df_params)) & grepl("dyn",rownames(df_params)) & grepl("01",rownames(df_params))),
                      which(grepl("cepi",rownames(df_params)) & grepl("dyn",rownames(df_params)) & grepl("02",rownames(df_params))),
                      which(grepl("cepi",rownames(df_params)) & grepl("dyn",rownames(df_params)) & grepl("03",rownames(df_params))),
                      which(grepl("cepi",rownames(df_params)) & grepl("dyn",rownames(df_params)) & grepl("04",rownames(df_params))))){
    
    print(paste0("r",row_vec[1]))
    
    sims=c(sims,
           run_sims_epidemic(parameters_df = df_params[row_vec,],
                             time_vect = Time_vect,
                             n_bootstrap = bootstrap ,
                             num_cores = n_cores,
                             status="Epidemic",
                             lambda_pattern = NULL,
                             pre_run=T,
                             foi_distribution=foi_dist))
  }
  
  sims_temp = list()
  for(bs in 1:bootstrap){
    sims_temp[[bs]]=1
    for(id in 1:10){sims_temp[[bs]]=c(sims_temp[[bs]],sims[[(1:length(sims))[((id-1)*bootstrap)+bs]]])}
    sims_temp[[bs]] = sims_temp[[bs]][-1]
  }
  
  sims=sims_temp
  
  bc_pop_size = df_params$pop_size==1e7
  
  df = rbind(get_output_df(list_bs=sims,
                           parameters_df=df_params,
                           pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                           selec_scenarios=bc_pop_size,
                           district_size =10e6,
                           status = "Epidemic"),
             
             get_output_df(list_bs=sims,
                           parameters_df=df_params,
                           pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                           selec_scenarios=!bc_pop_size,
                           district_size =1e6,
                           status = "Epidemic"))
  
  #Endemic-----------------------------
}else if(all_countries$status[task_id]=="Endemic"){
  
  df_params$Lambda_outbreak[grepl("cepi" ,rownames(df_params))] = df_params$Lambda_outbreak_endemic[grepl("cepi" ,rownames(df_params))]
  df_params$Var_outbreak[grepl("cepi" ,rownames(df_params))] = df_params$Var_outbreak_endemic[grepl("cepi" ,rownames(df_params))]
  
  abacus_long_term_end = abacus_long_term %>% filter(classification=="Endemic")
  get_beta_longterm_end = stats::approxfun(abacus_long_term_end$foi, y = abacus_long_term_end$beta0, method="linear", rule = 2)
  abacus_long_term_epi = abacus_long_term %>% filter(classification=="Epidemic")
  get_beta_longterm_epi = stats::approxfun(abacus_long_term_epi$foi, y = abacus_long_term_epi$beta0, method="linear", rule = 2)
  
  beta_vect = NULL

  #Vaccine impact
  #Remove sensitivity analysis in endemic countries for computational purposes
  select_vect = which(!grepl("sens",rownames(df_params)))
  df_params = df_params[select_vect,]
  
  sims = run_sims_epidemic(parameters_df = df_params[grepl("cepi", rownames(df_params)),],
                           time_vect =Time_vect,
                           n_bootstrap = bootstrap,
                           status="Endemic",
                           lambda_pattern = NULL,
                           num_cores = n_cores,
                           pre_run = T,
                           foi_distribution=foi_dist)
  
  sims_gavi = run_sims_epidemic(parameters_df = df_params[grepl("gavi", rownames(df_params)),],
                                time_vect =Time_vect,
                                n_bootstrap = bootstrap,
                                status="Epidemic",
                                lambda_pattern = NULL,
                                num_cores = n_cores,
                                pre_run = T,
                                foi_distribution=foi_dist)
  
  bc_pop_size = df_params$pop_size[grepl("cepi", rownames(df_params))]==1e7
  df = get_output_df(list_bs=sims,
                     parameters_df=df_params[grepl("cepi", rownames(df_params)),],
                     pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                     selec_scenarios=bc_pop_size,
                     district_size =10e6,
                     status = "Endemic")
  
  bc_pop_size = df_params$pop_size[grepl("gavi", rownames(df_params))]==1e7
  df_gavi = get_output_df(list_bs=sims,
                          parameters_df=df_params[grepl("gavi", rownames(df_params)),],
                          pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                          selec_scenarios=bc_pop_size,
                          district_size =10e6,
                          status = "Epidemic")
  
  sims_temp = list()
  for(bs in 1:bootstrap){
    sims_temp[[bs]]=c(sims[[bs]],sims_gavi[[bs]])
  }
  
  sims = sims_temp
  df = rbind(df, df_gavi)
  
  
  #NoTransmission---------------------
}else {
  beta_vect = NULL

  sims = run_sims_epidemic(parameters_df = df_params,
                           time_vect = Time_vect,
                           num_cores = n_cores,
                           n_bootstrap = bootstrap,
                           status = "NoTransmission",
                           pre_run = T,
                           foi_distribution=foi_dist)
  
  bc_pop_size = df_params$pop_size==1e7
  
  df = rbind(get_output_df(list_bs=sims,
                           parameters_df=df_params,
                           pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                           selec_scenarios=bc_pop_size,
                           district_size =10e6,
                           status = "NoTransmission"),
             
             get_output_df(list_bs=sims,
                           parameters_df=df_params,
                           pop_size=sum(demog_vect)*all_countries$prop_eff_pop[task_id],
                           selec_scenarios=!bc_pop_size,
                           district_size =1e6,
                           status = "NoTransmission"))
}

#Save data----
save(df_params, df, sims, file=paste0("3_Vaccine_impact/output_by_country/",task_id,"_",all_countries$country[task_id],".Rdata"))

print(Sys.time()-time_start)

