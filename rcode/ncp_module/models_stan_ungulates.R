### Running through different models to better understand NCPs and also how to use the Stan model blocks
# Started by Cat on 26 January 2021

### For more information from Stan users manual follow this link: https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html

### Simulation Question:
# During peak migration season, how far do ungulates travel per day?
# Predictor: 0/1 if American or Canadian. Canadian's are 1s


# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())


# Load Libraries
library(RColorBrewer)
library(viridis)
library(lme4)
library(ggplot2)
library(gridExtra)
library(rstan)
library(shiny)

# Set Working directory
setwd("~/Documents/git/la/rcode/ncp_learning/")

source("simdata.R")

#### This sim code has...
# Intercept = 300
# mu_b_canadian_sp = 50
# mu_b_herb_sp = -20
# mu_b_ch_sp = 5
# sigma_b_canadian_sp = 10
# sigma_b_canadian_sp = 5
# sigma_b_ch_sp = 5
# sigma_y = 2

### Prep for all muplots
cols <- adjustcolor("indianred3", alpha.f = 0.3) 
my.pal <-rep(viridis_pal(option="viridis")(9),2)
my.pch <- rep(15:18, each=10)
alphahere = 0.4

datalist.simple <- with(distall, 
                     list(y = distance, 
                          canadian = canadian, ### for simple: 
                          herbivore = herbivore,
                          sp = as.numeric(as.factor(species)),
                          N = nrow(distall),
                          n_sp = length(unique(distall$species))
                     )
)


##################################################################
########## MODEL 1: simple, no NCP with interaction ##############
##################################################################
### okay so for model 1, this is simple, no interaction and no NCP
simple = stan('stan/ungulates_simple.stan', data = datalist.simple,
                      iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(simple)$summary

modelhere <- simple
source("muplot.R")


##################################################################
################ MODEL 1: NCP, with interaction ##################
##################################################################
#NOTE: This model is slower... interesting! 
ncpsimple = stan('stan/ungulates_ncp_simple.stan', data = datalist.simple,
              iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(ncpsimple)$summary

modelhere <- ncpsimple
source("muplot.R")

### Uh-oh! mu plot is a mess! Those species level estimates are terribly wrong. 
# Okay let's see if we can fix it...
noncps <- modoutput[!grepl("_raw", rownames(modoutput)),]
source("muplot_ncp.R")

## Phew! All set now. But let's see if we can get more accurate estimates (fix the ESS issue) by adding in more priors...
ncptweak = stan('stan/ungulates_ncp_tweak.stan', data = datalist.simple,
                 iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(ncptweak)$summary
noncps <- modoutput[!grepl("_raw", rownames(modoutput)),]
source("muplot_ncp.R")


###### Cool! Well that works. But isn't kind of annoying scrolling through ALL of the yhats to see your model results?
## Let's start playing around with the blocks
# Sandbox 1: move yhats to model block rather than transformed parameters. A little faster!
ncp_rmyhat = stan('stan/ungulates_ncp_yhat.stan', data = datalist.simple,
                iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(ncpyhat)$summary
noncps <- modoutput[!grepl("_raw", rownames(modoutput)),]
source("muplot_ncp.R")

# But Wait! We really should look at those yhats for checking out our posterior predictive checks
# Sandbox 2: Let's add a generated quantities block and remove the ncp just to keep playing
simple_yhatmod_genquant = stan('stan/ungulates_yhat_genquant.stan', data = datalist.simple,
                  iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(simple_yhatmod_genquant)$summary
source("muplot.R")

# Sandbox 3: Finally, let's just test the timing between using transformed parameter block vs generated quantities
simple_yhatmod_transpara = stan('stan/ungulates_yhat_transpara.stan', data = datalist.simple,
                               iter = 1000, warmup=500, chains=2) ### , control=list(adapt_delta=0.99, max_treedepth=15)

modoutput <- summary(simple_yhatmod_transpara)$summary
source("muplot.R")

##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### And the timing is...##### ##### ##### ##### ###
##### ##### ##### ##### ##### ##### ##### ##### ##### 


get_elapsed_time(simple)
##warmup  sample
#chain:1 8.21052 2.41631
#chain:2 8.14878 1.95720

get_elapsed_time(ncpsimple)
##warmup  sample
#chain:1 65.6045 18.7121
#chain:2 41.0756 17.8176

get_elapsed_time(ncptweak)
##warmup  sample
#chain:1 34.7202 21.5377
#chain:2 24.0256 18.3064

get_elapsed_time(ncp_rmyhat)
##warmup  sample
#chain:1 36.6496 20.7792
#chain:2 35.7168 23.0576

get_elapsed_time(simple_yhatmod_genquant)
#warmup  sample
#chain:1 8.75932 1.88479
#chain:2 6.58924 1.99405

get_elapsed_time(simple_yhatmod_transpara)
## warmup  sample
#chain:1 7.94259 2.04537
#chain:2 7.60744 1.99668

