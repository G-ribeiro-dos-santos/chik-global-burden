//Spatial correlation based on https://aheblog.com/2016/12/07/geostatistical-modelling-with-r-and-stan/
//https://betanalpha.github.io/assets/case_studies/gaussian_processes.html
//https://natelemoine.com/fast-gaussian-process-models-in-stan/
//http://gaussianprocess.org/gpml/chapters/RW.pdf (chapter 2 - focusing on 2.7!)

data {
  int<lower=0> A; //the number of years
  int<lower=0> A2; //the number of maximum foi groups
  int<lower=0> R1; //the number of locations with data
  array[A2, R1] int<lower=0> N_tot; //the number of individuals
  array[A2, R1] int<lower=0> N_pos; //the number of seropositives
  array[A2, R1] int<lower=0> age_min; // Age min per group per location
  array[A2, R1] int<lower=0> age_max; // Age max per group per location
  array[R1] int<lower=0> max_group_loc; // Number of foi groups per location
  // int <lower=1, upper=A2> ind_by_age[A]; //
  // int<lower=0> YearsFixed; // Number of years with fixed lambda
  // real rho_fixed;
  
  //int <lower=0> R2; //the number of new points(for prediction)
  // matrix[R1+R2,R1+R2] dist; //distances between points
  //row_vector[2] coords[R1]; //coordinates of observed points
  //row_vector[2] coords_pred[R2]; //coordinates of points to predict
}
transformed data {
  //int<lower=1> R;
  //R = R1 + R2;
}
parameters {
  real<lower=0, upper=1> lambda_hyper;
  real<lower=0, upper=1> var_hyper;
  array[R1] real<lower=0, upper=1> lambda1; //
}
transformed parameters {
  array[A, R1] real cum_foi; // cumulative foi by age
  array[A2, R1] real cum_foi_by_group; //cumulative foi by age group
  
  for (r1 in 1 : R1) {
    cum_foi[1, r1] = lambda1[r1];
    for (j in 2 : A) {
      cum_foi[j, r1] = cum_foi[j - 1, r1] + lambda1[r1];
    }
    
    for (a2 in 1 : max_group_loc[r1]) {
      cum_foi_by_group[a2, r1] = mean(cum_foi[(1 + age_min[a2, r1]) : (
                                      1 + age_max[a2, r1]), r1]);
    }
  }
}
model {
  target += beta_lpdf(lambda_hyper | 1, 1);
  target += beta_lpdf(var_hyper | 1, 1);
  
  for (r1 in 1 : R1) {
    target += beta_lpdf(lambda1[r1] | ((1 - lambda_hyper) / var_hyper
                                       - (1 / lambda_hyper))
                                      * lambda_hyper ^ 2, (((1 - lambda_hyper)
                                                            / var_hyper
                                                            - (1
                                                               / lambda_hyper))
                                                           * lambda_hyper ^ 2)
                                                          * (1 / lambda_hyper
                                                             - 1));
    
    for (a2 in 1 : max_group_loc[r1]) {
      N_pos[a2, r1] ~ binomial(N_tot[a2, r1],
                               1 - exp(-cum_foi_by_group[a2, r1]));
    } //Serocatalytic model
  }
}
