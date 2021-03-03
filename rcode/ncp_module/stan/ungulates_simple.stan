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
  
  vector[n_sp] b_canadian; // slope of canadian effect 
  vector[n_sp] b_herb; // slope of herb length effect 
  vector[n_sp] b_ch;
  
	}

transformed parameters {
  vector[N] yhat;
  
  for(i in 1:N){    
    yhat[i] = a_sp[sp[i]] + // indexed with species
		          b_canadian[sp[i]] * canadian[i] + 
		          b_herb[sp[i]] * herbivore[i] +
		          b_ch[sp[i]] *  inter_canadianherb[i];
	      }
	      
}

model {
	a_sp ~ normal(mu_a_sp, sigma_a_sp); 
	
	target += normal_lpdf(to_vector(b_canadian) | mu_b_canadian_sp, sigma_b_canadian_sp); // just another way to write normal()
	target += normal_lpdf(to_vector(b_herb) | mu_b_herb_sp, sigma_b_herb_sp);
	target += normal_lpdf(to_vector(b_ch) |  mu_b_ch_sp, sigma_b_ch_sp);
	
        mu_a_sp ~ normal(400, 75);
        sigma_a_sp ~ normal(0, 50);
        
        sigma_y ~ normal(0, 100);
	      
	y ~ normal(yhat, sigma_y);

}

