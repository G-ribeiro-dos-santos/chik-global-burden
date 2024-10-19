library(dplyr)
library(ggplot2)

rm(list=ls())

#Subsets
list_files = list.files(path = "1_Serocatalytic_models/epidemic outputs/", full.names = T)

load(list_files[[1]])
list_subsets
list_subsets = c("No_Africa","No_America","No_Asia","Asia","America","Africa")

df_estimates = lapply(list_files, function(x){load(x);estimates}) %>% do.call(what=rbind) %>% filter(mu1!=0.1) %>% mutate(subset=list_subsets[excluded_cont]) #REMOVING COMPLEMENT RUNS

#Normalise likelihood
df_estimates = df_estimates %>% group_by(subset) %>%
  mutate(loglik_norm = loglik-max(loglik))

#Best pairs
best_pair = df_estimates %>% group_by(excluded_cont) %>% filter(loglik==max(loglik))

df_estimates %>%

  ggplot(aes(x=1/mu1, y=lambda1, fill=loglik_norm, z=loglik_norm))+
  geom_raster()+
  scale_x_log10()+
  facet_wrap(~subset, nrow=2) +
  
  geom_point(data=best_pair %>% 
               filter(subset %in% c("Global", "Asia", "America", "Africa")) %>%
               mutate(subset = factor(subset, levels=c("Global","Africa","Asia","America"))), aes(x=1/mu1, y=lambda1), col="green") +
  
  geom_vline(data=best_pair %>% 
               filter(subset %in% c("Global", "Asia", "America", "Africa")) %>%
               mutate(subset = factor(subset, levels=c("Global","Africa","Asia","America"))), aes(xintercept=1/mu1), col="green", linetype="dashed") +
  
  geom_hline(data=best_pair %>% 
               filter(subset %in% c("Global", "Asia", "America", "Africa")) %>%
               mutate(subset = factor(subset, levels=c("Global","Africa","Asia","America"))), aes(yintercept=lambda1), col="green", linetype="dashed") +
  
  scale_fill_gradient2(name="Distance to\nmax. LL",midpoint = -100, low = "black", mid = "purple", high = "orange", limits=c(-200,0), na.value = 1)+
  theme_minimal() +
  xlab("Outbreak freq.") +
  ylab("Outbreak FOI") +
  theme(plot.title = element_text(size=25)) +
  theme(axis.title.x = element_text(size=20), axis.text.x = element_text(angle=0, size=15)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(angle=0, size=15)) +
  theme(legend.title = element_text(color = "black", size = 20), legend.text = element_text(color = "black", size = 15)) +
  #theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(0.3, 'inches')) + 
  theme(strip.text = element_text(size=15))

df_estimates_crop = df_estimates
df_estimates_crop = df_estimates %>% filter(max(loglik)-loglik<200)

#Likelihood profiles
sp=0.5
mu_vect = seq(0.03,1,by=0.0001)
log_mu = lapply(unique(df_estimates$subset), function(sbset) predict(loess(loglik~mu1,  df_estimates_crop %>% filter(subset==sbset) %>% group_by(mu1) %>% filter(rank(-loglik)<=2) %>% ungroup() , span=sp), newdata=mu_vect))
mu_range = sapply(log_mu, function(log_mu) c(mu_vect[which.max(log_mu)], range(mu_vect[which(max(log_mu, na.rm=T)-log_mu<1.92)]))) %>% 
  t() %>% as.data.frame() %>% rename(max=V1,lower=V2,upper=V3) %>% mutate(subset=unique(df_estimates$subset))

lambda_vect = seq(0.01,1,by=0.0001)
log_lambda = lapply(unique(df_estimates$subset), function(sbset) predict(loess(loglik~lambda1,  df_estimates_crop %>% filter(subset==sbset) %>% group_by(lambda1) %>% filter(rank(-loglik)<=2) %>% ungroup() , span=sp), newdata=lambda_vect))
lambda_range = sapply(log_lambda, function(log_lambda) c(lambda_vect[which.max(log_lambda)], range(lambda_vect[which(max(log_lambda, na.rm=T)-log_lambda<1.92)]))) %>%
  t() %>% as.data.frame() %>% rename(max=V1,lower=V2,upper=V3) %>% mutate(subset=unique(df_estimates$subset))

p1 = df_estimates_crop %>%
  group_by(mu1, subset) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() %>%
  ggplot()+
  facet_wrap(~subset, scales="free_y") +
  geom_point(aes(x=mu1,y=loglik))+
  geom_smooth(aes(x=mu1,y=loglik), method="loess",span=sp)+
  scale_x_log10()+
  geom_vline(data=mu_range, aes(xintercept =max))+
  geom_vline(data=mu_range, aes(xintercept =lower), linetype="dashed")+
  geom_vline(data=mu_range, aes(xintercept =upper), linetype="dashed")+
  theme_minimal() +
  scale_x_continuous(name="Outbreak freq.") +
  ylab("Log. Lik") +
  theme(plot.title = element_text(size=25)) +
  theme(axis.title.x = element_text(size=20), axis.text.x = element_text(angle=0, size=15)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(angle=0, size=15)) +
  theme(legend.title = element_text(color = "black", size = 20), legend.text = element_text(color = "black", size = 15)) +
  theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(0.3, 'inches')) + 
  theme(strip.text = element_text(size=15))

p2 = df_estimates_crop %>%
  group_by(lambda1, subset) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() %>%
  ggplot()+
  facet_wrap(~subset, scales="free_y") +
  geom_point(aes(x=lambda1,y=loglik))+
  geom_smooth(aes(x=lambda1,y=loglik), method="loess",span=sp)+
  geom_vline(data=lambda_range, aes(xintercept =max))+
  geom_vline(data=lambda_range, aes(xintercept =lower), linetype="dashed")+
  geom_vline(data=lambda_range, aes(xintercept =upper), linetype="dashed")+
  theme_minimal() +
  scale_x_continuous(name="Outbreak FOI") +
  ylab("Log. Lik") +
  theme(plot.title = element_text(size=25)) +
  theme(axis.title.x = element_text(size=20), axis.text.x = element_text(angle=0, size=15)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(angle=0, size=15)) +
  theme(legend.title = element_text(color = "black", size = 20), legend.text = element_text(color = "black", size = 15)) +
  theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(0.3, 'inches')) + 
  theme(strip.text = element_text(size=15))

cowplot::plot_grid(p1, p2)

par_estimates = rbind(mu_range%>%mutate(var="mu"), lambda_range%>%mutate(var="lambda")) %>% tidyr::pivot_wider(id_cols = subset, names_from = var, values_from = c(max, lower, upper))

par_clean = par_estimates %>% mutate(lambda = paste0(max_lambda, " [", lower_lambda, "-", upper_lambda,"]"),
                                     mu = paste0(max_mu, " [", lower_mu, "-", upper_mu,"]")) %>%
  select(subset, mu , lambda)

par_paper = par_estimates[c(2,3,4,5),]
par_paper$subset = c("Asia", "America", "Africa", "Asia &\nAfrica")

p_paper = par_paper %>% 
  ggplot(aes(x=max_mu, y=max_lambda, label=subset, col=subset=="Global")) + 
  geom_point() + 
  geom_errorbar(aes(xmin=lower_mu, xmax=upper_mu)) +
  geom_errorbar(aes(ymin=lower_lambda, ymax=upper_lambda)) +
  ggrepel::geom_text_repel(fontface="bold")+
  scale_y_continuous(name="Epidemic FOI")+
  expand_limits(y=c(0,0.55),x=c(0.09,0.75)) +
  scale_x_log10(name="Annual\nprobability") +
  scale_color_manual(values=c(1,2))+
  theme_minimal() +
  theme(legend.position="none") +
  theme(plot.title = element_text(size=25)) +
  theme(axis.title.x = element_text(size=20), axis.text.x = element_text(angle=0, size=15)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(angle=0, size=15)) +
  theme(legend.title = element_text(color = "black", size = 20), legend.text = element_text(color = "black", size = 15)) +
  theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(0.3, 'inches')) + 
  theme(strip.text = element_text(size=10))

p_paper

write.csv(par_clean, file ="1_Serocatalytic_models/estimates_epidemic_formatted.csv")
ggsave(p_paper, file="3_Plots/plots_report/sensitivity_continent_estimates_smc.pdf")
