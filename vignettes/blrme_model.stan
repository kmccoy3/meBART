data {
    int<lower=0> N;       // number of cases
    vector[N] y;          // outcome (variate)
    vector[N] x_meas;   // measurement of x
    real<lower=0> tau;     // measurement noise

    int<lower=0> N_test;
    vector[N_test] x_test;

    // int<lower=0> N_func;
    // vector[N_func] x_func;
    // 
    // int<lower=0> N_nonoise;
    // vector[N_nonoise] x_nonoise;
}
parameters {
    vector[N] x;    // unknown true value
    real mu_x;          // prior location
    real sigma_x;       // prior scale
    real alpha;           // intercept
    real beta;            // slope
    real<lower=0> sigma;  // outcome noise
}
model {
    x ~ normal(mu_x, sigma_x);  // prior
    x_meas ~ normal(x, tau);    // measurement model
    y ~ normal(alpha + beta * x, sigma);
    alpha ~ normal(0, 10);
    beta ~ normal(0, 10);
    sigma ~ cauchy(0, 5);
}
generated quantities {
   vector[N_test] y_test;
   for(i in 1:N_test){
    y_test[i] = alpha+beta*x_test[i];
   }
   
   //    for(i in 1:N_test){
   //  y_test[i] = normal_rng(alpha+beta*x_test[i], sigma);
   // }

   //  vector[N_func] y_func;
   //  for(j in 1:N_func){
   //      y_func[j] = normal_rng(alpha+beta*x_func[j], sigma);
   // }
   // 
   //  vector[N_nonoise] y_nonoise;
   //  for(k in 1:N_nonoise){
   //      y_nonoise[k] = normal_rng(alpha+beta*x_nonoise[k], sigma);
   // }
}
