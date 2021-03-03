// NCP and stan model block learning
// Started by Cat 26 January 2021
// Level: Species on INTERCEPTS and SLOPES

data {
  
	int<lower=1> N;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	vector[N] y; 		// response
	vector[N] canadian; 	// canadian predictor
	vector[N] herbivore; 	// herbivore predictor
		
	}
	
transformed data {
  vector[N] inter_canadianherb;                 

  inter_canadianherb = canadian .* herbivore; 
}

parameters {
  
  real mu_a_sp;   
  real mu_b_canadian_sp;     
  real mu_b_herb_sp;
  real mu_b_ch_sp; // slope of canadian x herbivore effect
  
  real<lower=0> sigma_b_canadian_sp;
  real<lower=0> sigma_b_herb_sp;
  real<lower=0> sigma_b_ch_sp;
  real<lower=0> sigma_a_sp;
  real<lower=0> sigma_y; 
  
  real a_sp[n_sp]; // intercept for species
  
  vector[n_sp] b_canadian_raw; // slope of canadian effect 
  vector[n_sp] b_herb_raw; // slope of herb length effect 
  vector[n_sp] b_ch_raw;
  
	}

transformed parameters {
  vector[N] yhat;

  vector[n_sp] b_canadian = mu_b_canadian_sp + sigma_b_canadian_sp*b_canadian_raw; 
  vector[n_sp] b_herb = mu_b_herb_sp + sigma_b_herb_sp*b_herb_raw; 
  vector[n_sp] b_ch = mu_b_ch_sp + sigma_b_ch_sp*b_ch_raw;
  
  for(i in 1:N){    
    yhat[i] = a_sp[sp[i]] + // indexed with species
		          b_canadian[sp[i]] * canadian[i] + 
		          b_herb[sp[i]] * herbivore[i] +
		          b_ch[sp[i]] *  inter_canadianherb[i];
	      }
	      
}

model {
	b_canadian_raw ~ normal(0, 1);
	b_herb_raw ~ normal(0, 1);
	b_ch_raw ~ normal(0, 1);
	
	a_sp ~ normal(mu_a_sp, sigma_a_sp); 
	
	/// NOTE ON NCP: since we used NCP, the b_canadian etc are assumed to be partially pooled, if we did 
	 // b_canadian(mu_b_canadian_sp, sigma_b_canadian_sp) it would confuse the model because it would be redundant
	  
  target += normal_lpdf(to_vector(b_canadian) | 0, 50); // just another way to write normal()
	target += normal_lpdf(to_vector(b_herb) | 0, 30);
	target += normal_lpdf(to_vector(b_ch) | 0, 10);
	      
        mu_a_sp ~ normal(300, 20);
        sigma_a_sp ~ normal(0, 20);
        
        sigma_y ~ normal(0, 10);
	      
	y ~ normal(yhat, sigma_y);

}
