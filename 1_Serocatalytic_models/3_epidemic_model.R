library(iterators)
library(tidyr)
library(dplyr)
library(foreach)
library(pomp)
library(spatPomp)
library(doParallel)
library(doRNG)
library(ggplot2)

rm(list=ls())

time_start = Sys.time()

# This script was written to run on a multi-core CPU from an HPC - task_id is an argument passed on from the bash script being the array ID from the job array (from 1 to 198)

#Arguments passed
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#print(paste0("Available cores : ", availableCores()))

#Getting nb of cores used for parallelising
print(paste0("detectCores : ", detectCores(), " cores"))
registerDoParallel(cores=detectCores())

#Making sure we get diff seeds for each parallelised task
registerDoRNG(13456123)

#Parameters--------------------------------------
N_p= 1e3 #Number of particles
N_rep = 20 #Nb. of filter replicates to get likelihood estimate and uncertainty
ncores = detectCores()

range_mu = exp(seq(log(0.02), log(1), length.out = 50))
range_lambda = seq(0.05, 0.5,length.out=50)

fixedParNames = c("mu","lambda")
expandedParNames = NULL

list_subsets = list("Africa","America",c("Asia","Oceania"),c("Africa","America"),c("Africa","Asia","Oceania"),c("America","Asia","Oceania"))
excluded_continent = list_subsets[[((task_id-1)%%6)+1]]

#Import the data---------------------------------------------------
load("1_Serocatalytic_models/stan_input_global.Rdata")

tapply(sub_seroprev$Region, list(sub_seroprev$continent, sub_seroprev$classification), function(x) length(unique(x)))

data_sub =  sub_seroprev[sub_seroprev$classification=="Epidemic" & !sub_seroprev$continent %in% excluded_continent,]

colnames(data_sub)[colnames(data_sub)=="N_seropos"] = "N_pos"
data_sub$Age_mean = floor(0.5*(data_sub$Age_max+data_sub$Age_min))
data_sub = data_sub[data_sub$Age_mean<=30,]

data_sub$sub_ID = as.numeric(as.factor(paste0(data_sub$Region,data_sub$ID)))
data_sub = data_sub[!is.na(data_sub$N_tot),]
data_sub = data_sub[data_sub$N_tot!=0,]

R1 = length(unique(data_sub$sub_ID))

data_sub$N_pos = as.integer(data_sub$N_pos) ; data_sub$N_tot = as.integer(data_sub$N_tot)
data_sub = pivot_wider(data_sub, id_cols=Age_mean, names_from=sub_ID, values_from=c(N_pos,N_tot))
data_sub = data_sub[order(data_sub$Age_mean),c("Age_mean",paste0("N_pos_",1:R1),paste0("N_tot_",1:R1))]

#Filtering and Prediction functions--------------------------
rinit_c = spatPomp_Csnippet (unit_statenames = c("foi","cum_foi"),
                             code = "for (int u=0; u < U; u++){
                       foi[u] = 0;
                       cum_foi[u] = 0;}
                   ")

rprocess_c = discrete_time(spatPomp_Csnippet (method = "rprocess", unit_statenames = c("foi","cum_foi"), unit_paramnames=c("mu","lambda"),code = "
                                     for (int u=0; u < U; u++){
                                     if(runif(0,1)<mu[u*lambda_unit]){
                                     foi[u] = rbeta(10,(10/lambda[u*lambda_unit])-10);
                                     }else{
                                     foi[u]=0;
                                     }
                                     cum_foi[u] = cum_foi[u] + foi[u];
                                     }
                                     "),
                           delta.t = 1)

rmeasure_c = spatPomp_Csnippet (method = "rmeasure", unit_statenames =  c("cum_foi"),
                                unit_obsnames =  c("N_pos","N_tot"), code = "
                      for (int u=0; u < U; u++){
                      N_tot[u] = 100;
                      N_pos[u] = rbinom(N_tot[u], 1-exp(-cum_foi[u]));}
                      ")


dmeasure_c = spatPomp_Csnippet(method = "dmeasure", unit_statenames =  c("cum_foi"),  unit_obsnames =  c("N_pos","N_tot"), 
                               code="
       lik = 0;
      for (int u=0; u < U; u++){
        if(ISNA(N_pos[u])){
          lik = lik;
        }else{
          lik = lik + dbinom(N_pos[u], N_tot[u], 1-exp(-cum_foi[u]-1e-6), 1);
        }
      }
        if(!give_log) lik = exp(lik);
                      ")

dunit_measure_c = spatPomp_Csnippet (code="
       lik = 0;
      if(ISNA(N_pos)){
          lik = lik;
        }else{
          lik = lik + dbinom(N_pos, N_tot, 1-exp(-cum_foi-1e-6), 1);
       }
          if(!give_log) lik = exp(lik);
                      ")

#Building the spatPOMP model-----------------------------------------------------
data_test = data_sub %>% tidyr::pivot_longer(cols = -Age_mean, names_to = "variable", values_to = "values") %>%
  mutate(Age_mean=as.integer(Age_mean),
         loc = paste0("loc_", stringr::str_split_i(variable,"_",-1)),
         variable = paste0("N_",stringr::str_split_i(variable,"_",2))) %>%
  tidyr::pivot_wider(id_cols = c(Age_mean, loc), names_from = variable, values_from = values) %>% as.data.frame()

paramnames <- c(if (length(fixedParNames) > 0) {
  paste0(fixedParNames, "1")
}, if (length(expandedParNames) > 0) {
  paste0(rep(expandedParNames, each = R1), 1:R1)
})

set_expanded <- Csnippet(paste0("const int ", expandedParNames, 
                                "_unit = 1;\n", collapse = " "))
set_fixed <- Csnippet(paste0("const int ", fixedParNames, 
                             "_unit = 0;\n", collapse = " "))
globals1 <- Csnippet(paste(set_expanded, set_fixed, 
                           sep = "\n"))

mod1 =  spatPomp(data_test, times="Age_mean", t0=0, units="loc", #Data, time and spatial units
                 
                 rinit=rinit_c, 
                 
                 rprocess=rprocess_c,
                 
                 rmeasure=rmeasure_c,
                 
                 dunit_measure=dunit_measure_c,
                 
                 unit_statenames=c("foi","cum_foi"),
                 
                 paramnames=paramnames,
                 
                 globals  = globals1,
                 
                 partrans=parameter_trans(logit=paramnames)) #Variable declaration and transformation

#MLE search
rw_char = lapply(expandedParNames, function(par) paste0(par,1:R1,"=0.2", collapse = ","))
rw_char = paste0("rw_sd(", do.call(what=paste0, args=list(rw_char, collapse=",")), ")")
rw_vect = eval(parse(text=rw_char))

ranges = as.matrix(data.frame(mu=c(0.05, 0.25), lambda=c(0.05, 0.25)))

starts = expand.grid(mu=range_mu, lambda=range_lambda)

n_complement = ncores - (nrow(starts) %% ncores)
starts = rbind(starts, data.frame(mu=rep(0.1, n_complement), lambda=0.4))

start_pars = lapply(1:nrow(starts),  function(x) eval(parse(text=paste0("c(",paste0("lambda",1:R1,"=",starts$lambda[x],collapse=","),",",paste0("mu",1:R1,"=",starts$mu[x],collapse=","),")")))) %>% do.call(what=rbind)

idx_task = (1:ncores) + ncores * ((task_id-1)%/%6)


estimates = foreach (pars=idx_task,
                     .combine=rbind, .packages=c("spatPomp"),
                     .errorhandling="remove", .inorder=FALSE) %dopar% {
                       ll = logmeanexp(replicate(logLik(bpfilter(mod1,
                                                                 params=start_pars[pars,],
                                                                 block_size=1,
                                                                 Np=N_p)),n=N_rep),se=T)
                       data.frame(as.list(start_pars[pars,c("lambda1","mu1")]),loglik=ll[1],loglik.se=ll[2],excluded_cont=((task_id-1)%%6)+1)
                     }

save(start_pars, list_subsets, estimates, file=paste0("1_Serocatalytic_models/epidemic_outputs/output_bpf_map_subset_",task_id,".Rdata"))

print(Sys.time()-time_start)
