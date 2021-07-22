### Sample simulation code for understanding NCPs and stan model blocks
## Let's pretend to look at migration patterns in animals
# During peak migration season, how far do animals travel per day?
## Predictor: 0/1 if American or Canadian. Canadian's are 1s

set.seed(12321)

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### For faster analyses and better fit simulation data, let's follow Geoff's code ##### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

#  1) Let's make the observations much higher than the actual data to build a good model.
nsp = 10 # number of species
ntot = 100 # numbers of obs per species. 

sample_a <- list(can.env = rbinom(1000, 1, 0.5),
                 herb.env = rbinom(1000, 1, 0.5))

model.parameters <- list(intercept = 300,
                         can.coef = 50,
                         herb.coef = -20,
                         canxherb = 5)

#  2) Now, we will make varying intercepts
env.samples <- sapply(sample_a, FUN = function(x){
  sample(x, size = nsp * ntot, replace = TRUE)})

# Determine which environmental variables interact
intrxnname <- names(model.parameters)[4] # interaction terms
names.temp <- gsub("x", "|", intrxnname) # remove text to align with colnames
env.pairs <- sapply(1:length(names.temp), FUN = function(X){
  grep(pattern = names.temp[X], x = colnames(env.samples))
})
# Add these interactions (product) to env.samples        
env.interactions <- sapply(1:ncol(env.pairs), FUN = function(X){
  apply(env.samples[, env.pairs[, X]], MARGIN = 1, FUN = prod)
})
env.samples2 <- cbind(env.samples, env.interactions)
# Create model matrix
mm <- model.matrix(~env.samples2)

#  4) We need to make a random intercept model for each species
parameters.temp <- matrix(unlist(model.parameters), ncol = length(model.parameters), nrow = nsp * ntot, byrow = TRUE)

# Which parameters are random?
random.regex <- grep(pattern = paste(c("intercept", "can.coef", "herb.coef", "canxherb"), collapse = "|"), x = names(model.parameters))

# Generate random parameters (by species)
parameters.temp[, 1] <- sapply(1:nsp, FUN = function(x){
  rep(rnorm(n = 1, mean = model.parameters[[random.regex[1]]], sd = 10), ntot)})
parameters.temp[, 2] <- sapply(1:nsp, FUN = function(x){
  rep(rnorm(n = 1, mean = model.parameters[[random.regex[2]]], sd = 10), ntot)})
parameters.temp[, 3] <- sapply(1:nsp, FUN = function(x){
  rep(rnorm(n = 1, mean = model.parameters[[random.regex[3]]], sd = 5), ntot)})
parameters.temp[, 4] <- sapply(1:nsp, FUN = function(x){
  rep(rnorm(n = 1, mean = model.parameters[[random.regex[4]]], sd = 2), ntot)})
# Calculate response
response <- sapply(1:nrow(env.samples), FUN = function(x){
  rnorm(n = 1, mean = mm[x, ] %*% parameters.temp[x, ], sd = 2)})

distall <- cbind(data.frame(species = as.vector(sapply(1:nsp, FUN = function(x) rep(x, ntot))),
                                               distance = response, canadian = env.samples[,1], herbivore = env.samples[,2]))

#write.csv(testdata_provmethod_intrxn, file="output/testdata_provmethod_intrxn.csv", row.names = FALSE)

#  7) Let's do a quick lmer model to test the fake data
modtest <- lmer(distance ~ canadian + herbivore + canadian*herbivore + (canadian + herbivore + canadian*herbivore|species), data=distall) ## Quick look looks good!


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
if(FALSE){
  # Step 1: Set up years, days per year, temperatures, sampling frequency, required GDD (fstar)
  nspps <- 10
  ninds <- 100 
  nobs <- nspps*ninds
  
  ### These are our response variables
  y_dist <- 300  ### mu_a_sp in model output
  y_dist_speciessd <- 50 ### sigma_a_sp in model output
  
  ## Sigma_y to be added at the end
  sigma_y <- 2 
  
  #### Incorporate Predictors 
  can_effect <- 20  ## Canadian effect, this is saying that if sites are from 1 degree north, they travel 20 fewer kms
  can_sd <- 5 ## Canadian effect sd
  
  distspp <- round(rnorm(nspps, y_dist, y_dist_speciessd), digits=0)
  df.dist <- as.data.frame(cbind(species=rep(1:nspps, each=ninds), ind=rep(1:ninds), 
                                 distspp=rep(distspp, each=ninds)))
  
  df.dist$canadian <- rbinom(nobs, 1, 0.5)
  
  df.dist$dist.can <- df.dist$canadian * rep(rnorm(n=nspps, mean=can_effect, sd=can_sd), each=ninds)
  
  df.dist$distance <- df.dist$distspp + df.dist$dist.can + rnorm(n=nobs, mean=0, sd=sigma_y)
  
  ##### Clean up the data frame to prepare for analyses
  distall <- df.dist[!duplicated(df.dist),]
}

