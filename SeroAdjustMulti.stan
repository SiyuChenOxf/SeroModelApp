data {
  int n_days;
  int n_days2;
  int daily_death[n_days];
  int sero[n_days2];
  int t[n_days];
  int t2[n_days2];
}

parameters {
   real <lower= 0,upper=1> beta;  
   real <lower= 0,upper=1> gamma;
}

model {  
  real store_value;
  int  index;
  
  beta~ beta(1,1);
  gamma~ beta(1,1); 

  
  store_value = 0;
  index = 1;
  for (i in 1:n_days2) {
    for (n in index:t2[i]){
      store_value = store_value + exp(beta*t[n])*(1-gamma)/gamma*daily_death[n];
    }
    sero[i] ~ neg_binomial_2(store_value/exp(beta*t2[i]),100);
    index = t2[i]+1;
  };

}




