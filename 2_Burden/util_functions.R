get_yll = function(age_vect, n_infections_vect, parms, regular_life_expectancy, cfr_df, infected_to_case){
  
  prop_infected_that_die = sapply(floor(age_vect), function(age) infected_to_case * cfr_df$cfr[which(cfr_df$age_min<=age & cfr_df$age_max>=age)] * parms$mortality)
  
  years_lost_death =  sapply(age_vect, function(age) max(regular_life_expectancy-age,0))
  
  yll = prop_infected_that_die * years_lost_death * n_infections_vect
  
  return(yll)
  
}

get_yld = function(age_vect, n_infections_vect, parms, regular_life_expectancy, cfr_df, infected_to_case){
  
  prop = c(Asymptomatic=1-infected_to_case,mild=infected_to_case*0.5,moderate=infected_to_case*0,severe=infected_to_case*0.5,chronic=infected_to_case*parms$case_to_chronic)
  
  disability_weight = c(Asymptomatic=0, mild=0.006, moderate=0.051, severe=0.133, chronic=0.233)
  
  duration = c(Asymptomatic=0,mild=6/365,moderate=6/365,severe=6/365,chronic=1)
  
  years_lost = disability_weight*duration*prop
  
  yld = years_lost*n_infections_vect
  
  return(yld)
}

get_deaths = function(age_vect, n_infections_vect, parms, cfr_df, infected_to_case){
  
  prop_infected_that_die = sapply(floor(age_vect), function(age) infected_to_case * cfr_df$cfr[which(cfr_df$age_min<=age & cfr_df$age_max>=age)] * parms$mortality)
  
  n_death = prop_infected_that_die * n_infections_vect
  
  return(n_death)
}

get_beta0 = function(foi, S0, abacus_foi2){
  #Looking at abacus doing super basic linear interpolation
  lower_row_S0 = findInterval(S0, as.numeric(rownames(abacus_foi2)))
  if(lower_row_S0==1){
    vect_foi_along_beta0 = abacus_foi2[lower_row_S0+1,]
  }
  else if (lower_row_S0==nrow(abacus_foi2)){
    vect_foi_along_beta0 = abacus_foi2[lower_row_S0,]
  }else{
    weight_lower_S0 = 1 - (S0-as.numeric(rownames(abacus_foi2))[lower_row_S0]) / diff(as.numeric(rownames(abacus_foi2))[c(lower_row_S0,lower_row_S0+1)])
    vect_foi_along_beta0 = weight_lower_S0* abacus_foi2[lower_row_S0,] + (1-weight_lower_S0) * abacus_foi2[lower_row_S0+1,]
  }
  beta0 = approx(y=as.numeric(colnames(abacus_foi2)), x=vect_foi_along_beta0, xout=foi)$y
  beta0 = ifelse(is.na(beta0), as.numeric(-log(1-foi)/(foi*S0)), beta0)
  return (beta0)
}

sir_model_age_epidemic <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #Time-varying beta
    T_outbreak = 1 #year (for abacus)
    beta =  ifelse(time<T_outbreak,25+2*beta0*(1-time/T_outbreak), 25)

    #Vectorise compartments
    S <- state[grep("^S_",names(state))]
    I <- state[grep("^I_[0-9]",names(state))]
    R <- state[grep("^R_",names(state))]
    V <- state[grep("^V_",names(state))]
    IV <- state[grep("^IV_[0-9]",names(state))]
    RV <- state[grep("^RV_",names(state))]
    Icum = state[grep("^Icum",names(state))]
    IVcum = state[grep("^IVcum",names(state))]
    
    #Vectorise parameters
    target_age = unlist(parameters[grep("target_age",names(parameters))])
    
    #Vaccination
    if(status=="Epidemic"){
      if(sum(Icum)>(vacc_threshold*pop_size/prob_symptoms)){flag_infection[bs]<<-T} #CAREFUL ABOUT SYNTAX HERE <<- different from <-
      if(flag_infection[bs] & time_infection[bs]==-1e3){time_infection[bs] <<- time}
      if((time-time_infection[bs])>v_delay & time <= (time_infection[bs] + v_delay + d_uptake) & #If within vacc campaign window
         sum(I)>0.1*(vacc_threshold*pop_size/prob_symptoms)){ #and if outbreak is still going on
        vaccinate = daily_vaccines}else{vaccinate=0}
    }else if(status=="Endemic"){vaccinate = ifelse(time>v_delay & time <= (v_delay + d_uptake),daily_vaccines,0)} #If within vacc campaign window
    
    
    #ODE
    # Total population
    N <- sum(S) + sum(I) + sum(R) + sum(V) + sum(IV) + sum(RV)
    
    # Force of infection
    lambda <- beta * (sum(I)+sum(IV)) / N
    
    #Vaccination
    S_to_V = S*0
    S_to_V[target_age] = S_to_V[target_age]+vaccinate*S[target_age]/(sum(S[target_age])+sum(R[target_age]))
    R_to_RV = R*0
    R_to_RV[target_age] = vaccinate*R[target_age]/(sum(S[target_age])+sum(R[target_age]))
    
    #Infection
    S_to_I = lambda*S
    V_to_IV = lambda*V*(1-v_eff_inf)
    
    #Recovery
    I_to_R = sigma*I
    IV_to_RV = sigma*IV
    
    dS = - S_to_I - S_to_V
    dI = S_to_I - I_to_R
    dIV = V_to_IV - IV_to_RV
    dR = I_to_R - R_to_RV
    dV = S_to_V - V_to_IV
    dRV = R_to_RV + IV_to_RV
    
    ddoses = vaccinate
    dIcum = S_to_I
    dIVcum = V_to_IV
    
    return(list(c(dS,dI,dIV,dR,dV,dRV,dIcum,dIVcum,ddoses)))
  })
}

sir_model_age_epidemic_twodose <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #Time-varying beta
    T_outbreak = 1 #year (for abacus)
    beta = ifelse(time<T_outbreak,25+2*beta0*(1-time/T_outbreak), 25)

    
    #Vectorise compartments
    S <- state[grep("^S_",names(state))]
    I <- state[grep("^I_[0-9]",names(state))]
    R <- state[grep("^R_",names(state))]
    V <- state[grep("^V_",names(state))]
    IV <- state[grep("^IV_[0-9]",names(state))]
    RV <- state[grep("^RV_",names(state))]
    V1 = state[grep("^V1_",names(state))]
    V2 = state[grep("^V2_",names(state))]
    Icum = state[grep("^Icum",names(state))]
    IVcum = state[grep("^IVcum",names(state))]
    
    #Vectorise parameters
    target_age = unlist(parameters[grep("target_age",names(parameters))])
    
    #Vaccination
    if(status=="Epidemic"){
      if(sum(Icum)>(vacc_threshold*pop_size/prob_symptoms)){flag_infection[bs]<<-T} #CAREFUL ABOUT SYNTAX HERE <<- different from <-
      if(flag_infection[bs] & time_infection[bs]==-1e3){time_infection[bs] <<- time}
      if((time-time_infection[bs])>v_delay & time <= (time_infection[bs] + v_delay + d_uptake) & #If within vacc campaign window
         sum(I)>0.1*(vacc_threshold*pop_size/prob_symptoms)){ #and if outbreak is still going on
        vaccinate = daily_vaccines}else{vaccinate=0}
    }else if(status=="Endemic"){vaccinate = ifelse(time>v_delay & time <= (v_delay + d_uptake),daily_vaccines,0)} #If within vacc campaign window
    
    
    #ODE
    # Total population
    N <- sum(S) + sum(I) + sum(R) + sum(V) + sum(IV) + sum(RV)
    
    # Force of infection
    lambda <- beta * (sum(I)+sum(IV)) / N
    
    #1st dose
    S_to_V1 = S*0
    S_to_V1[target_age] = S_to_V1[target_age]+vaccinate*S[target_age]/(sum(S[target_age])+sum(R[target_age]))
    R_to_RV = R*0
    R_to_RV[target_age] = vaccinate*R[target_age]/(sum(S[target_age])+sum(R[target_age]))
    
    #2nd dose
    V1_to_V2 = V1/((28+v_delay)/365) #28 days til 2nd dose
    
    #Infection
    S_to_I = lambda*S
    V1_to_IV = lambda*V1*(1-0.5*v_eff_inf)
    V2_to_IV = lambda*V2*(1-v_eff_inf)
    
    #Recovery
    I_to_R = sigma*I
    IV_to_RV = sigma*IV
    
    dS = - S_to_I - S_to_V1
    dI = S_to_I - I_to_R
    dIV = V1_to_IV + V2_to_IV - IV_to_RV
    dR = I_to_R - R_to_RV
    
    dV1 = S_to_V1 - V1_to_IV - V1_to_V2
    dV2 = V1_to_V2 - V2_to_IV
    dV = dV1 + dV2
    
    dRV = R_to_RV + IV_to_RV
    
    ddoses = sum(S_to_V1 + V1_to_V2 + 2*R_to_RV)
    dIcum = S_to_I
    dIVcum = V1_to_IV + V2_to_IV

    return(list(c(dS,dI,dIV,dR,dV,dRV,dIcum,dIVcum,ddoses,dV1,dV2)))
  })
}

run_sims_epidemic = function(parameters_df, num_cores, time_vect, n_bootstrap, lambda_pattern=NULL, status, pre_run, pre_run_status=NULL, foi_distribution="beta"){

  n_scenarios = nrow(parameters_df)
  
  #Transformed parameters - Beta distribution of outbreaks mean Lambda and variance Var
  alpha = ((1-parameters_df$Lambda_outbreak[1])/parameters_df$Var_outbreak[1]-(1/parameters_df$Lambda_outbreak[1]))*parameters_df$Lambda_outbreak[1]^2
  beta = alpha*(1/parameters_df$Lambda_outbreak[1]-1)
  
  list_bs = mclapply(X=1:n_bootstrap, mc.cores=num_cores, mc.set.seed=T, FUN = function(bs) {
    
    tryCatch({
      lambda_fixed = 1
      if(foi_distribution=="normal"){ #Has the country reliable case data
        while(lambda_fixed>0.4 | lambda_fixed<0){lambda_fixed = rnorm(1, parameters_df$Lambda_outbreak[1] , sqrt(parameters_df$Var_outbreak[1]))} #Avoid pb with long tail distributions (especially when base foi is really low)
      }else{
        while(lambda_fixed>0.4){lambda_fixed = rbeta(1, alpha ,beta)} #Avoid pb with long tail distributions (especially when base foi is really low)
      }
      
      #Get a fixed beta in some scenarios (whether it's for abacus creation or longer term runs for vaccine impact)
      if(status=="Endemic"){
        lambda_vect_pre = rep(1,30)
        if(pre_run){ #vaccine impact
          beta_fixed_pre = get_beta_longterm_end(lambda_fixed)
          beta_fixed = get_beta_longterm_end(lambda_fixed)
          beta0_vect =rep(beta_fixed,n_scenarios)
        }else if(!is.null(beta_vect)){ #abacus creation
          beta_fixed_pre = beta_vect[1]
          beta0_vect = beta_vect
        }
      }
      if(status=="Epidemic"){
        if(is.null(pre_run_status)){lambda_vect_pre = rep(c(1,0,0,0,0,0,0,0),15)[1:30]
        }else if(pre_run_status=="Endemic"){lambda_vect_pre = rep(1,30)
        }else{lambda_vect_pre = rep(c(1,0,0,0,0,0,0,0),15)[1:30]
        }
        if(pre_run){
          beta_fixed_pre = get_beta_longterm_epi(lambda_fixed)
          beta_fixed = get_beta_longterm_epi(lambda_fixed)
          beta0_vect =rep(beta_fixed,n_scenarios)
        }else if(!is.null(beta_vect)){
          beta_fixed_pre = beta_vect[1]
          beta0_vect = beta_vect
        }
      }
      if(status=="NoTransmission"){
        lambda_vect_pre = rep(0,30)
        beta_fixed_pre = 0
        beta_fixed = 0
        beta0_vect =0
      }
      
      S_0_init = ifelse(status=="NoTransmission",1,S_0_init_who)
      
      if(!is.null(beta_vect) | pre_run){S_0_init = 1} #If creating abacus or if vaccine impact estimate
      
      age_groups_df = data.frame(min=c(0,5,10,12,18,20,30,40,50,60,70,80)+1,
                                 max=c(5,10,12,18,20,30,40,50,60,70,80,100))
      
      names_vect = as.vector(sapply(c("S","I","IV","R","V","RV"), function(x) paste0(x, "_", stringr::str_pad(1:nrow(age_groups_df),2,pad="0"))))
      
      
      Population_df =  data.frame(Susceptible_Population = parameters_df$pop_size[1] *  Age_structure * S_0_init,
                                  Infected_Population = 0 * Age_structure,
                                  Infected_Vaccinated_Population = 0 * Age_structure,
                                  Recovered_Population = Age_structure * (1-S_0_init) * parameters_df$pop_size[1], 
                                  Vaccinated_Population = 0 * Age_structure,
                                  Recov_Vaccinated_Population = 0 * Age_structure)
      
      Birth_cohort = data.frame(Susceptible_Population = round(parameters_df$pop_size[1] *  Age_structure[1]),
                                Infected_Population = 0,
                                Infected_Vaccinated_Population = 0,
                                Recovered_Population = 0,
                                Vaccinated_Population=0,
                                Recov_Vaccinated_Population=0)
      
      #Start recording
      Record_pop = lapply(1:n_scenarios,function(x){
        tmp = c(colSums(Population_df),rep(0,27))
        cats = c("infections","cases","deaths","dalys","ylls","ylds")
        names(tmp) = c(colnames(Population_df), cats, paste0(cats,"_averted"), paste0(cats,"_u5"), paste0(cats,"_u5_averted"), "doses", "foi","beta0")
        return(tmp)})
      
      if(pre_run){
        print(paste0("pre_run_",bs))
        #Run endemically without vaccination for 30 years to get to equilibrium point
        for(t in 1:30){ #time_vect
          if(lambda_vect_pre[t]>0){
            #Introduce infections (1 infection per 1e6 ppl - to be consistent with abacus)
            intro_inf = parameters_df$pop_size[1]/1e6
            Population_df[1,c("Susceptible_Population","Infected_Population")] = Population_df[1,c("Susceptible_Population","Infected_Population")] + c(-intro_inf, intro_inf) #S to I
            
            init_state = c()
            for(a in 1:nrow(age_groups_df)){
              out = colSums(Population_df[age_groups_df$min[a]:age_groups_df$max[a],])
              names(out)=paste0(c("S","I","IV","R","V","RV"),"_",stringr::str_pad(a,2,pad="0"))
              init_state=c(init_state,out)
            }
            init_state = c(init_state[names_vect],rep(0,2*nrow(age_groups_df)+1))
            names(init_state)=c(names_vect,
                                paste0("Icum_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0")),
                                paste0("IVcum_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0")),
                                "doses")
            
            #Run ODE system
            flag_infection[bs] <- F
            time_infection[bs] <- -1e3
            
            output = as.data.frame(ode(y=init_state, 
                                       func = sir_model_age_epidemic,
                                       c(parameters_df[1,],
                                         status=status,
                                         bs=bs,
                                         beta0=beta_fixed_pre,
                                         target_age= rep(T, nrow(age_groups_df)),
                                         sigma=parameters_df$Sigma[1],
                                         daily_vaccines=0),
                                       times = seq(0,2,by=1/365)))
            
            #empty the I compartment
            output[,grepl("^Icum_", colnames(output))] = output[,grepl("^I_", colnames(output))] + output[,grepl("^Icum_", colnames(output))]
            output[,grepl("^R_", colnames(output))] = output[,grepl("^I_", colnames(output))] + output[,grepl("^R_", colnames(output))]
            output[,grepl("^I_", colnames(output))] = 0
            #empty the IV compartment
            output[,grepl("^IVcum_", colnames(output))] = output[,grepl("^IV_", colnames(output))] + output[,grepl("^IVcum_", colnames(output))]
            output[,grepl("^RV_", colnames(output))] = output[,grepl("^IV_", colnames(output))] + output[,grepl("^RV_", colnames(output))]
            output[,grepl("^IV_", colnames(output))] = 0
            
            #Update Population
            final_state = output[nrow(output),colnames(output)!="time"][names(init_state)]
            infected_by_age_group_by_status = final_state[grepl("cum",names(final_state))]
            final_state = final_state[!grepl("cum",names(final_state))]
            
            updated_pop = Birth_cohort #Just to initialise, will remove first row after the loop 
            for(a in 1:nrow(age_groups_df)){
              pop_temp = (as.matrix(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]/sum(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]))) %*%
                (as.matrix(final_state[grepl(stringr::str_pad(a,2,pad="0"),names(final_state))]))
              colnames(pop_temp)=names(Birth_cohort)
              updated_pop = rbind(updated_pop, pop_temp)
            }
            updated_pop = updated_pop[-1,]
            Population_df = updated_pop
            
            foi = sum(infected_by_age_group_by_status)/sum(output[1,grepl("^V_",colnames(output))] + output[1,grepl("^S_",colnames(output))])
            Record_pop = lapply(1:n_scenarios,function(x){rbind(Record_pop[[x]], c(colSums(Population_df),rep(0,25),foi,beta_fixed_pre))})
          }else{
            Record_pop = lapply(1:n_scenarios,function(x){rbind(Record_pop[[x]], c(colSums(Population_df),rep(0,27)))})
          }
          #Deaths and birth of population
          new_pop = rbind(Birth_cohort, Population_df[-nrow(Population_df),])
          age_structure = rowSums(new_pop)/sum(new_pop)
          adjustment = Age_structure/age_structure
          adjustment[is.nan(adjustment)] = 0 #Solve countries like Eritrea that have 0 pop in 99 bin
          new_pop = apply(new_pop, 2, FUN =  function(y) y*adjustment*parameters_df$pop_size[1]/sum(new_pop))
          Population_df = as.data.frame(new_pop)
        }
      }
      
      #Create one population per scenario
      Population_df = rep(list(Population_df),n_scenarios)
      Birth_cohort = rep(list(Birth_cohort), n_scenarios)
      #Scale to population size of each parameters set
      Population_df = lapply(1:length(Population_df), function(x) Population_df[[x]]/sum(Population_df[[x]])*parameters_df$pop_size[x])
      Birth_cohort = lapply(1:length(Birth_cohort), function(x) Birth_cohort[[x]]/sum(Birth_cohort[[x]])*parameters_df$pop_size[x]*Age_structure[1])
      
      #Running the simulation
      if(status=="NoTransmission"){lambda_pattern=c(rep(0,length(time_vect)))}
      if(status=="Endemic"){lambda_pattern=rep(1,length(time_vect))}
      if(is.null(lambda_pattern)){lambda_pattern = (runif(length(time_vect),0,1)<parameters_df$Mu_outbreak[1])} #Draw outbreaks beforehand
      
      print(paste0("run_",bs))
      for(t in time_vect){ #time_vect
        #print(t)
        #Loss of vaccination coverage
        Population_df = lapply(1:length(Population_df), function(x) {
          y = Population_df[[x]]
          y[,c("Susceptible_Population",
               "Vaccinated_Population",
               "Recovered_Population",
               "Recov_Vaccinated_Population")] = 
            
            y[,c("Susceptible_Population",
                 "Vaccinated_Population",
                 "Recovered_Population",
                 "Recov_Vaccinated_Population")] +
            
            c((Population_df[[x]]$Vaccinated_Population*(1-exp(-1/parameters_df$v_dur[x]))), #V to S
              - (Population_df[[x]]$Vaccinated_Population*(1-exp(-1/parameters_df$v_dur[x]))), #V to S
              (Population_df[[x]]$Recov_Vaccinated_Population*(1-exp(-1/parameters_df$v_dur[x]))), #RV to R
              -(Population_df[[x]]$Recov_Vaccinated_Population*(1-exp(-1/parameters_df$v_dur[x])))) #RV to R
          
          return(y)})
        
        outbreak_flag = F
        if(lambda_pattern[t]>0){
          
          outbreak_flag = T
          
          #Introduce infections (1 infection per 1e6 ppl - to be consistent with abacus)
          intro_inf = sapply(1:n_scenarios, function(x) parameters_df$pop_size[x]/1e6)
          
          Population_df = lapply(1:n_scenarios, function(x) {
            y = Population_df[[x]]
            y[1,c("Susceptible_Population","Infected_Population")] = 
              y[1,c("Susceptible_Population","Infected_Population")] +
              c(-intro_inf[x], intro_inf[x]) #S to I
            return(y)})
          
          init_state = lapply(Population_df, function(x){
            ret=c()
            for(a in 1:nrow(age_groups_df)){
              out = colSums(x[age_groups_df$min[a]:age_groups_df$max[a],])
              names(out)=paste0(c("S","I","IV","R","V","RV"),"_",stringr::str_pad(a,2,pad="0"))
              ret=c(ret,out)
            }
            ret = c(ret[names_vect],rep(0,2*nrow(age_groups_df)+1))
            names(ret)=c(names_vect,
                         paste0("Icum_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0")),
                         paste0("IVcum_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0")),
                         "doses")
            ret
          })
          
          #Get target pop for vaccines
          targ_age = lapply(1:n_scenarios, function(x){if(parameters_df$target_group[x]=="6m"){
            target_age = age_groups_df$min>0
          }else if(parameters_df$target_group[x]=="12yo"){
            target_age = age_groups_df$min>12
          }else if(parameters_df$target_group[x]=="18yo"){
            target_age = age_groups_df$min>18
          }else if(parameters_df$target_group[x]=="12_60_end"){
            if(status=="Endemic"){target_age = age_groups_df$min>12 & age_groups_df$max<=60}else{target_age = age_groups_df$min>12}
          }else if(parameters_df$target_group[x]=="60_end"){
            if(status=="Endemic"){target_age = age_groups_df$min>60}else{target_age = age_groups_df$min>12}
          }
            return(target_age)})
          
          #Get the number of daily vaccines during vaccine campaign period (depends on target group)
          target_pop = lapply(1:n_scenarios, function(x){
            target_age =  (sapply(1:sum(targ_age[[x]]), function(a) age_groups_df[targ_age[[x]],]$min[a]:age_groups_df[targ_age[[x]],]$max[a])) %>% do.call(what=c)
            
            vulnerable_pop = ifelse(parameters_df$vuln_pop[x]=="preg_and_imm", 0, sum(parameters_df$pop_size[x]*prop_no_vacc_vect[target_age]))
            #Blind to infection?
            if(parameters_df$infec_status_at_vac[x]=="All"){
              target_pop = max(sum(rowSums(Population_df[[x]][target_age,c("Susceptible_Population","Recovered_Population")]))-vulnerable_pop,0)
            }
            else if (parameters_df$infec_status_at_vac[x]=="Uninfected"){
              target_pop = max(sum((Population_df[[x]][target_age,c("Susceptible_Population")]))-vulnerable_pop,0)
            }
            #No vaccination
            target_pop = ifelse(is.infinite(parameters_df$vacc_threshold[x]),0,target_pop)
            return(target_pop)
          })
          
          #Total population targeted (including already vaccinated)
          target_pop_tot = lapply(1:n_scenarios, function(x){
            target_age =  (sapply(1:sum(targ_age[[x]]), function(a) age_groups_df[targ_age[[x]],]$min[a]:age_groups_df[targ_age[[x]],]$max[a])) %>% do.call(what=c)
            
            vulnerable_pop = ifelse(parameters_df$vuln_pop[x]=="preg_and_imm", 0, sum(parameters_df$pop_size[x]*prop_no_vacc_vect[target_age]))
            #Blind to infection?
            if(parameters_df$infec_status_at_vac[x]=="All"){
              target_pop = max(sum(rowSums(Population_df[[x]][target_age,]))-vulnerable_pop,0)
            }
            else if (parameters_df$infec_status_at_vac[x]=="Uninfected"){
              target_pop = max(sum((Population_df[[x]][target_age,]))-vulnerable_pop,0)
            }
            #No vaccination
            target_pop = ifelse(is.infinite(parameters_df$vacc_threshold[x]),0,target_pop)
            return(target_pop)
          })
          
          daily_vaccines = lapply(1:n_scenarios, function(x){
            #Linear coverage
            parameters_df$vacc_coverage[x]*target_pop[[x]]/(parameters_df$d_uptake[x])
          })
          
          #Is it a vaccination year in endemic scenarios? (frequency equal to duration of vaccination)
          if(status=="Endemic"){daily_vaccines = lapply(1:n_scenarios, function(x) ifelse((((t-1)%%parameters_df$vaccination_frequency[x])==0), daily_vaccines[[x]], 0))}
          
          #Get beta_0 according to lambda_fixed (for burden estimates)
          if(!pre_run & is.null(beta_vect)){
            beta0_vect = sapply(1:n_scenarios, function(x) get_beta0(foi=lambda_fixed,
                                                                     S0 = sum(init_state[[x]][grepl("^S_",names(init_state[[x]]))] +
                                                                                init_state[[x]][grepl("^V_",names(init_state[[x]]))])/
                                                                       sum(init_state[[x]]),
                                                                     abacus_foi2 = abacus_foi_12months))
          }
          
          #Run ODE system
          output = list()
          for(x in 1:nrow(parameters_df)){
            #print(paste0(t,"_",x))
            
            # Debug
            flag_infection[bs] <- F
            time_infection[bs] <- -1e3
            # 
            #flag_infection[bs] <<- F
            #time_infection[bs] <<- -1e3
            # 
            if(parameters_df$n_dose[x]==1){
              mod = as.data.frame(ode(y=init_state[[x]], 
                                      func = sir_model_age_epidemic,
                                      c(parameters_df[x,],
                                        status=status,
                                        bs=bs,
                                        beta0=beta0_vect[x],
                                        target_age= targ_age[[x]],
                                        sigma=parameters_df$Sigma[x],
                                        daily_vaccines=as.numeric(daily_vaccines[[x]])),
                                      times = seq(0,2,by=1/365)))
              
              
            }else if(parameters_df$n_dose[x]==2){
              
              init_state[[x]][paste0("V1_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0"))] = 0
              init_state[[x]][paste0("V2_",stringr::str_pad(1:nrow(age_groups_df),2,pad="0"))] = init_state[[x]][grep("^V_",names(init_state[[x]]))]
              
              mod = as.data.frame(ode(y=init_state[[x]], 
                                      func = sir_model_age_epidemic_twodose,
                                      c(parameters_df[x,],
                                        status=status,
                                        bs=bs,
                                        beta0=beta0_vect[x],
                                        target_age= targ_age[[x]],
                                        sigma=parameters_df$Sigma[x],
                                        daily_vaccines=as.numeric(daily_vaccines[[x]])),
                                      times = seq(0,2,by=1/365)))
              
              mod = mod[,!grepl("^V[1-2]",colnames(mod))]
              init_state[[x]] = init_state[[x]][!grepl("^V[1-2]",names(init_state[[x]]))]
            }
            
            #empty the I compartment
            mod[,grepl("^Icum_", colnames(mod))] = mod[,grepl("^I_", colnames(mod))] + mod[,grepl("^Icum_", colnames(mod))]
            mod[,grepl("^R_", colnames(mod))] = mod[,grepl("^I_", colnames(mod))] + mod[,grepl("^R_", colnames(mod))]
            mod[,grepl("^I_", colnames(mod))] = 0
            #empty the IV compartment
            mod[,grepl("^IVcum_", colnames(mod))] = mod[,grepl("^IV_", colnames(mod))] + mod[,grepl("^IVcum_", colnames(mod))]
            mod[,grepl("^RV_", colnames(mod))] = mod[,grepl("^IV_", colnames(mod))] + mod[,grepl("^RV_", colnames(mod))]
            mod[,grepl("^IV_", colnames(mod))] = 0
            
            output[[x]] = mod
          }
          
          #Update Population
          final_state = lapply(1:n_scenarios, function(x) output[[x]][nrow(output[[x]]),colnames(output[[x]])!="time"][names(init_state[[x]])])
          infected_by_age_group_by_status = lapply(final_state, function(x) x[grepl("cum",names(x))])
          final_state = lapply(final_state, function(x) x[!grepl("cum",names(x))])
          
          Population_df = lapply(1:n_scenarios, function(x){
            updated_pop = Birth_cohort[[x]] #Just to initialise, will remove first row after the loop 
            for(a in 1:nrow(age_groups_df)){
              pop_temp = (as.matrix(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]/sum(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]))) %*%
                (as.matrix(final_state[[x]][grepl(stringr::str_pad(a,2,pad="0"),names(final_state[[x]]))]))
              colnames(pop_temp)=names(Birth_cohort[[x]])
              updated_pop = rbind(updated_pop, pop_temp)
            }
            updated_pop = updated_pop[-1,]
            return(updated_pop)
          })
          
          #Measure stuff
          foi = lapply(1:n_scenarios, function(x) c(foi=sum(infected_by_age_group_by_status[[x]])/sum(output[[x]][1,grepl("^V_",colnames(output[[x]]))] + output[[x]][1,grepl("^S_",colnames(output[[x]]))])))
          
          
          infections_by_status = lapply(1:n_scenarios, function(x){
            updated_pop = rep(0,2)
            for(a in 1:nrow(age_groups_df)){
              pop_temp = (as.matrix(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]/sum(Age_structure[age_groups_df$min[a]:age_groups_df$max[a]]))) %*%
                (as.matrix(infected_by_age_group_by_status[[x]][grepl(stringr::str_pad(a,2,pad="0"),names(infected_by_age_group_by_status[[x]]))]))
              updated_pop = rbind(updated_pop, pop_temp)
            }
            return(updated_pop[-1,])
          })
          
          infections = lapply(infections_by_status, rowSums)
          
          ylls = lapply(1:n_scenarios, function(x){
            
            yll_unvaccinated  = get_yll(age_vect = 1:length(Age_structure)-0.5,
                                        n_infections_vect = infections_by_status[[x]][,1],
                                        parms = parameters_df[x,],
                                        regular_life_expectancy = life_expectancy,
                                        cfr_df = cfr,
                                        infected_to_case = parameters_df$prob_symptoms[x])
            
            yll_vaccinated =  get_yll(age_vect = 1:length(Age_structure)-0.5,
                                      n_infections_vect = infections_by_status[[x]][,2],
                                      parms = parameters_df[x,],
                                      regular_life_expectancy = life_expectancy,
                                      cfr_df = cfr,
                                      infected_to_case = parameters_df$prob_symptoms[x] * (1-parameters_df$v_eff_dis[x])/(1-parameters_df$v_eff_inf[x]))
            
            return(yll_vaccinated+yll_unvaccinated)})
          
          ylds = lapply(1:n_scenarios, function(x){
            
            yld_unvaccinated  = get_yld(age_vect = 1:length(Age_structure)-0.5,
                                        n_infections_vect = infections_by_status[[x]][,1],
                                        parms = parameters_df[x,],
                                        regular_life_expectancy = life_expectancy,
                                        cfr_df = cfr,
                                        infected_to_case = parameters_df$prob_symptoms[x])
            
            yld_vaccinated =  get_yld(age_vect = 1:length(Age_structure)-0.5,
                                      n_infections_vect = infections_by_status[[x]][,2],
                                      parms = parameters_df[x,],
                                      regular_life_expectancy = life_expectancy,
                                      cfr_df = cfr,
                                      infected_to_case = parameters_df$prob_symptoms[x] * (1-parameters_df$v_eff_dis[x])/(1-parameters_df$v_eff_inf[x]))
            
            return(yld_vaccinated+yld_unvaccinated)})
          
          dalys = lapply(1:n_scenarios, function(x){ylls[[x]]+ylds[[x]]})
          
          cases = lapply(1:n_scenarios, function(x){
            cases_unvaccinated = infections_by_status[[x]][,1]*parameters_df$prob_symptoms[x]
            cases_vaccinated =  infections_by_status[[x]][,2]*parameters_df$prob_symptoms[x] * (1-parameters_df$v_eff_dis[x])/(1-parameters_df$v_eff_inf[x])
            return(cases_unvaccinated+cases_vaccinated)})
          
          deaths = lapply(1:n_scenarios, function(x){
            
            deaths_unvaccinated = get_deaths(age_vect = 1:length(Age_structure)-0.5,
                                             n_infections_vect = infections_by_status[[x]][,1],
                                             parms = parameters_df[x,],
                                             cfr_df = cfr,
                                             infected_to_case = parameters_df$prob_symptoms[x])
            
            deaths_vaccinated = get_deaths(age_vect = 1:length(Age_structure)-0.5,
                                           n_infections_vect = infections_by_status[[x]][,2],
                                           parms = parameters_df[x,],
                                           cfr_df = cfr,
                                           infected_to_case = parameters_df$prob_symptoms[x]*(1-parameters_df$v_eff_dis[x])/(1-parameters_df$v_eff_inf[x]))
            
            return(deaths_vaccinated+deaths_unvaccinated)})
          
          #Summarise stuff
          stats_df = lapply(1:n_scenarios, function(x){
            cbind(infections=infections[[x]], cases=cases[[x]], deaths=deaths[[x]], dalys=dalys[[x]], ylls=ylls[[x]], ylds=ylds[[x]])
          })
          
          sum_stats = lapply(1:n_scenarios, function(x) colSums(stats_df[[x]]))
          
          sum_stats_under5 = lapply(1:n_scenarios, function(x) {
            tmp = colSums(stats_df[[x]][1:5,])
            names(tmp) = paste0(names(tmp),"_u5")
            tmp})
          
          sum_stats_averted = lapply(1:n_scenarios, function(x){
            if(parameters_df$pop_size[x]!=parameters_df$pop_size[1]){
              tmp = sum_stats[[1]]*parameters_df$pop_size[x]/parameters_df$pop_size[1]-sum_stats[[x]] #Adjust to pop_size
            }else{
              tmp = sum_stats[[1]]-sum_stats[[x]]
            }
            names(tmp) = paste0(names(tmp),"_averted")
            tmp}) #Assuming first scenarios is the control one
          
          sum_stats_under5_averted = lapply(1:n_scenarios, function(x){
            
            if(parameters_df$pop_size[x]!=parameters_df$pop_size[1]){
              tmp = sum_stats_under5[[1]]*parameters_df$pop_size[x]/parameters_df$pop_size[1]-sum_stats_under5[[x]]
            }else{
              tmp = sum_stats_under5[[1]]-sum_stats_under5[[x]]
            }
            
            names(tmp) = paste0(names(tmp),"_averted")
            tmp}) #Assuming first scenarios is the control one
          
          #Doses if we vaccinate everyone blindly
          doses = lapply(1:n_scenarios, function(x){
            
            doses_eff = c("doses"=final_state[[x]]$doses)
            
            doses_eff_target = parameters_df$vacc_coverage[x]*target_pop[[x]]*parameters_df$n_dose[x]
            
            doses_tot = (doses_eff/doses_eff_target)* #adjust in case campaign stops prematurely in epidemic scenarios
              doses_eff*(target_pop_tot[[x]]/target_pop[[x]])
            
            return(doses_tot)
          })
          
          #Record stuff
          Record_pop = lapply(1:n_scenarios,function(x){
            tmp = rbind(Record_pop[[x]],
                        c(colSums(Population_df[[x]]),
                          sum_stats[[x]],
                          sum_stats_averted[[x]],
                          sum_stats_under5[[x]],
                          sum_stats_under5_averted[[x]],
                          doses[[x]],
                          foi[[x]],
                          beta0_vect[x]))
            
            return(tmp)})
        }else{
          #Record no outbreak years in epi scenario
          Record_pop = lapply(1:n_scenarios,function(x){rbind(Record_pop[[x]], c(colSums(Population_df[[x]]),rep(0,27)))})
        }
        
        #Deaths and birth of population
        Population_df = lapply(1:n_scenarios , function(x) {
          new_pop = rbind(Birth_cohort[[x]], Population_df[[x]][-nrow(Population_df[[x]]),])
          age_structure = rowSums(new_pop)/sum(new_pop)
          adjustment = Age_structure/age_structure
          adjustment[is.nan(adjustment)] = 0 #Solve countries like Eritrea that have 0 pop in 99 bin
          new_pop = apply(new_pop, 2, FUN =  function(y) y*adjustment*parameters_df$pop_size[x]/sum(new_pop))
          return(as.data.frame(new_pop))
        })
        
      }
      Record_pop = lapply(Record_pop, function(x){
        tmp=as.data.frame(x)
        if(pre_run){year_vect = c(-30:0,time_vect)}else{year_vect=c(0,time_vect)}
        tmp = cbind(data.frame(Year=year_vect,tmp))})
      return(Record_pop)
    },error=function(e){e})
  })
  
  return(list_bs)}

get_output_df=function(list_bs, parameters_df, pop_size, district_size, status, selec_scenarios){

  n_districts = max(round(pop_size/district_size), 1)
  
  #If endemic get foi of subdistricts centered around a national foi else randomly sample
  if(status=="Endemic"){
    alpha = ((1-parameters_df$Lambda_outbreak[1])/parameters_df$Var_outbreak[1]-(1/parameters_df$Lambda_outbreak[1]))*parameters_df$Lambda_outbreak[1]^2
    beta = alpha*(1/parameters_df$Lambda_outbreak[1]-1)
    national_foi = rbeta(length(list_bs), alpha, beta)
    fois = sapply(list_bs, function(x)x[[1]]$foi[x[[1]]$Year==1])
    #Sample subdistricts based on proximity to national foi
    bs_index = as.vector((sapply(national_foi, function(nat_foi) sample(x = 1:length(fois), size = n_districts, replace = T, prob = dnorm(x =(nat_foi-fois)^2/parameters_df$Var_outbreak[1], mean = 0,sd = 1)))))
  }else{
    bs_index = sample(x = rep(1:length(list_bs), each=n_districts), replace = T, size = n_districts*length(list_bs))
  }
  list_bs = list_bs[bs_index]
  y_tmp = list_bs[[1]][[1]]
  sum_stats = array(unlist(lapply(list_bs, function(x) sapply(x[selec_scenarios], function(y) c(colSums(y[y$Year>0,colnames(y)!="Year"&!grepl("Population",colnames(y))]),
                                                                                                n_outbreaks=sum(y$foi>0))))),
                    
                    dim=c(1+sum(colnames(y_tmp)!="Year"&!grepl("Population",colnames(y_tmp))), sum(selec_scenarios), length(list_bs)))
  
  dimnames(sum_stats)[[1]] =  names(c(colSums(y_tmp[,colnames(y_tmp)!="Year"&!grepl("Population",colnames(y_tmp))]), n_outbreaks=sum(y_tmp$foi>0)))
  dimnames(sum_stats)[[2]] = rownames(parameters_df[selec_scenarios,])
  sum_stats = aperm(sum_stats, c(3,1,2))
  rm(y_tmp)
  n_years = 20
  df_ = lapply(dimnames(sum_stats)[[3]], function(vac_status){
    out = sum_stats[ , , vac_status] %>%
      as_tibble %>%
      
      #Group sum by n_districts iterations to simulate outbreaks in subdistricts
      mutate(grp=rep(c(1:(length(list_bs)/n_districts)),each=n_districts)) %>%
      select(-foi) %>%
      group_by(grp) %>%
      summarise(across(-n_outbreaks,~sum(.)*pop_size/(n_districts*district_size)),
                across(n_outbreaks,sum)) %>% ungroup() %>% select(-grp)
    
    Cols = grep('_averted$', colnames(out), value = T)
    
    #Proportion averted
    prop = out[ , Cols] / (out[ , Cols] + out[ , gsub('_averted$', '', Cols)])
    colnames(prop) = paste0("prop_",Cols)
    
    #Averted per dose
    per_dose = out[ , Cols] / out$doses
    colnames(per_dose) = paste0(Cols,"_per_dose")
    
    out = cbind(out,prop,per_dose) %>% mutate(Scenario = vac_status, .before = 1)
    
    return(out)
    
  }) %>%
    do.call(what = rbind)
  return(df_)
}