library(readxl)
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(stringr)

rm(list=ls())

check_cmdstan_toolchain(fix=T)

load("1_Serocatalytic_models/stan_input_global.Rdata")

#Get who regions
load("0_Data/demog_data_df.Rdata")
sub_seroprev= sub_seroprev %>% mutate(matched_name = all_countries$country[stringdist::amatch(Country, all_countries$country, method='jw', maxDist = Inf)])
sub_seroprev %>% filter(Country!=matched_name) %>% dplyr::select(Country, matched_name) %>% unique()
sub_seroprev = sub_seroprev %>% mutate(who_region = all_countries$who_region[sapply(matched_name, function(x) which(all_countries$country==x))])
sub_seroprev %>% dplyr::select(Country, who_region) %>% unique() %>% dplyr::select(who_region) %>% table()

#Select locations and restrict age
removed_areas = c("Odisha","West Bengal","Bihar","Tripura","Meghalaya","Assam","Punjab", #Remove places in India with no cirulation
                  "Kitui", "Central Nyanza", "Kainji Lakes") #Remove old serostudies


list_subsets = list("Africa","America",c("Asia","Oceania"),c("Africa","America"),c("Africa","Asia","Oceania"),c("America","Asia","Oceania"),c("Africa","America","Asia","Oceania"))

#Run MCMC
for(subset in list_subsets){
  print(subset)
  select_vect =  data$classification=="Epidemic" & !colnames(data$N_pos) %in% removed_areas & data$continent %in% subset
  age_max = 30
  
  #replace NA with zeroes
  data$N_tot[is.na(data$N_tot)] = 0
  data$N_pos[is.na(data$N_pos)] = 0

  #Run model
  data_sub = list(A = data$A,
                  A2 = data$A2,
                  R1 = sum(select_vect),
                  N_tot = data$N_tot[,select_vect],
                  N_pos = data$N_pos[,select_vect],
                  age_min = data$age_min[,select_vect],
                  age_max = data$age_max[,select_vect],
                  max_group_loc=data$max_group_loc[select_vect],
                  #YearsFixed=data$YearsFixed,
                  rho_fixed=data$rho_fixed)

  #Restrict to age_max
  data_sub$max_group_loc = apply(data_sub$age_max,2,function(x) max(which(x<age_max &x>0)))
  
  #Remove that dont have young age groups
  is_finite = !is.infinite(data_sub$max_group_loc)
  
  data_sub = list(A = data_sub$A,
                  A2 = data_sub$A2,
                  R1 = sum(is_finite),
                  N_tot = data_sub$N_tot[,is_finite],
                  N_pos = data_sub$N_pos[,is_finite],
                  age_min = data_sub$age_min[,is_finite],
                  age_max = data_sub$age_max[,is_finite],
                  max_group_loc=data_sub$max_group_loc[is_finite],
                  rho_fixed=data_sub$rho_fixed)

  mod <- cmdstan_model( "1_Serocatalytic_models/foi_endemic_global_hyperprior.stan", pedantic=F)

  list_init = list(lambda_hyper=0.02, var_hyper=0.01, lambda1 = rep(c(0.02),data_sub$R1))
  output <- mod$sample(data=data_sub, chains=4, parallel_chains=4, iter_sampling=2000, refresh=100, iter_warmup=500,
                       init=rep(list(list_init),4))
  output <- rstan::read_stan_csv(output$output_files())
  chains<-rstan::extract(output, inc_warmup=FALSE)
  save(chains, output, file=paste0("1_Serocatalytic_models/chains_epidemic_hyperprior_",paste0(subset,collapse=""),".Rdata"))
}


# Get endemic FOI
chains = lapply(list_subsets, function(subset){
  load(paste0("5_Stan/chains_epidemic_hyperprior_",paste0(subset,collapse=""),".Rdata"))
  chains})

names(chains) = lapply(list_subsets, function(subset) paste0(subset, collapse=""))

foi_fit = lapply(list_subsets, function(subset){
  select_vect =  data$classification=="Epiemic" & !colnames(data$N_pos) %in% removed_areas & data$continent %in% subset
  data_sub = list(A = data$A,
                  A2 = data$A2,
                  R1 = sum(select_vect),
                  N_tot = data$N_tot[,select_vect],
                  N_pos = data$N_pos[,select_vect],
                  age_min = data$age_min[,select_vect],
                  age_max = data$age_max[,select_vect],
                  max_group_loc=data$max_group_loc[select_vect],
                  #YearsFixed=data$YearsFixed,
                  rho_fixed=data$rho_fixed)

  subset = paste0(subset, collapse="")
  print(subset)
  foi = 1-exp(-chains[[subset]]$lambda1)
  foi = apply(foi,2,function(x) c(median(x), quantile(x,probs=c(0.025,0.975),na.rm=T)))
  lapply(1:ncol(foi), function(r1) data.frame(age=0:100) %>% mutate(Region=colnames(data_sub$N_tot)[r1],
                                                                    serop_mean=1-exp(-foi[1,r1]*age),
                                                                    serop_lower=1-exp(-foi[2,r1]*age),
                                                                    serop_upper=1-exp(-foi[3,r1]*age))) %>%
    do.call(what=rbind) %>% mutate(subset=subset)
}) %>% do.call(what=rbind)


#Hyperprior
lambda = sapply(list_subsets, function(subset){
  subset = paste0(subset, collapse="")
  quantile(chains[[subset]]$lambda_hyper,c(0.5,0.025,0.975))
}) %>% t() %>% as.data.frame() %>% rename(med = '50%', lo='2.5%', hi='97.5%') %>% mutate(subset = sapply(list_subsets, function(x) paste0(x,collapse="\n")))

lambda$subset = gsub(replacement = "",pattern = "Oceania", x=lambda$subset)


par_paper = lambda[c(1,2,3,5),]
par_paper$subset = c("Africa", "America","Asia", "Asia &\nAfrica")


p_paper = par_paper %>%
  mutate(subset = factor(subset, levels=subset)) %>%
  ggplot(aes(x=subset, y=med, ymin=lo, ymax=hi))+
  geom_pointrange(size=1.5, linewidth=1.5) +
  theme_minimal()+
  scale_y_log10() +
  ylab("Mean\nepidemic FOI") +
  xlab("Subset")+
  theme(plot.title = element_text(size=25)) +
  theme(axis.title.x = element_text(size=20), axis.text.x = element_text(size=15, angle = 45,hjust=1,vjust=1)) +
  theme(axis.title.y = element_text(size=20), axis.text.y =  element_text(size=15)) +
  expand_limits(y=0) +
  #theme(legend.position="none") +
  theme(legend.title = element_text(color = "black", size = 40), legend.text = element_text(color = "black", size = 30)) +
  theme(panel.grid.minor = element_blank())

p_paper

ggsave(p_paper, file="1_Serocatalytic_models/sensitivity_continent_estimates_mcmc_epidemic.pdf")
write.csv(lambda, file="1_Serocatalytic_models/estimates_endemic_in_epidemic.csv", row.names = F)
