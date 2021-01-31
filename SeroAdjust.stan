data {
  int<lower=0> N;  
  int<lower=0,upper=1> y[N]; 
}
parameters {
    // # sensitivity
    real<lower=0,upper=1> kse;
    
    // # specificity
      real<lower=0,upper=1> ksp;
      
      // probability of seropositivity 
      real<lower=0, upper=1> p;
}

model {
    kse ~ beta(11,1); //#mean=0.938;95 quantitle=0.995 this is expected to match the center and 95 quantitle of the sens range
    ksp ~ beta(15,1); //#mean=0.954;95 quantitle=0.996 this is expected to match the center and 95 quantitle of the spec range
    p ~ beta(1,1);
    y ~ bernoulli(kse*p+(1-ksp)*(1-p));
}