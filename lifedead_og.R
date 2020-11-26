# Illustrates an issue w/ local maxima in the likelihood when 
# fitting model combining live recapture and dead recoveries in Jags 
# using SSM formulation of multistate models

# If initial values are generated from the ld.init() function, we failed at recovering the values
# used to simulate the data; 

# everything goes well if we use the true latent states (data are simulated) as initial values 

# we carry out a classical freq analysis with several sets of inits and 
# demonstrate the presence of local minima in the likelihood

# Bayes with MCMC is not immune to the issue, see e.g. page 456-… in the gentle introduction to Mark.

# It remains to modify the ld.init() function to have more randomness in the inits generation.

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # To avoid error message if nothing to replace
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Multinomial trials for state transitions
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }




# Define mean survival, transitions, recapture, 
# as well as number of occasions, states, observations and released individuals 
ss <- 0.8
ff <- 0.6
rr <- 0.1
pp <- 0.5
n.occasions <- 10  
n.states <- 4
n.obs <- 3
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Releases in study area

set.seed(2020)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      ss*ff, ss*(1-ff), 1-ss, 0,
      0,   ss,       1-ss, 0,
      0,   0,       0,   1,
      0,   0,       0,   1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pp, 0, 1-pp,
      0, 0, 1,
      0, rr, 1-rr,
      0, 0, 1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH
zsimul <- sim$CH.TRUE

# Compute date of first capture
get.first <- function(x) min(which(x!=0))
first <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3


# Specify model in BUGS language
sink("lifedead.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# ss: true survival probability
# ff: fidelity probability
# rr: recovery probability
# pp: recapture/resighting probability
# -------------------------------------------------
# States (S):
# 1 alive in study area
# 2 alive outside study area
# 3 recently dead and recovered
# 4 recently dead, but not recovered, or dead (absorbing)
# Observations (O):
# 1 seen alive
# 2 recovered dead
# 3 neither seen nor recovered
# -------------------------------------------------

# Priors
ss ~ dunif(0, 1)     # Prior for mean survival
ff ~ dunif(0, 1)     # Prior for mean fidelity
rr ~ dunif(0, 1)     # Prior for mean recovery
pp ~ dunif(0, 1)     # Prior for mean recapture

# Define state-transition and observation matrices 	
   # Define probabilities of state S(t+1) given S(t)
      ps[1,1] <- ss * ff
      ps[1,2] <- ss * (1 - ff)
      ps[1,3] <- 1 - ss
      ps[1,4] <- 0
      ps[2,1] <- 0
      ps[2,2] <- ss
      ps[2,3] <- 1 - ss
      ps[2,4] <- 0
      ps[3,1] <- 0
      ps[3,2] <- 0
      ps[3,3] <- 0
      ps[3,4] <- 1
      ps[4,1] <- 0
      ps[4,2] <- 0
      ps[4,3] <- 0
      ps[4,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,1] <- pp
      po[1,2] <- 0
      po[1,3] <- 1 - pp
      po[2,1] <- 0
      po[2,2] <- 0
      po[2,3] <- 1
      po[3,1] <- 0
      po[3,2] <- rr
      po[3,3] <- 1 - rr
      po[4,1] <- 0
      po[4,2] <- 0
      po[4,3] <- 1

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,first[i]] <- y[i,first[i]]
   for (t in (first[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], 1:4])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], 1:3])
      } #t
   } #i
}
",fill = TRUE)
sink()


# Bundle data
jags.data <- list(y = rCH, 
                  first = first, 
                  n.occasions = dim(rCH)[2], 
                  nind = dim(rCH)[1])

# Initial values
ld.init <- function(ch, f){
   ch[ch==3] <- NA
   v2 <- which(ch==2, arr.ind = T)
   ch[v2] <- 3
   for (i in 1:nrow(v2)){
      ifelse(v2[i,2]!=ncol(ch), ch[v2[i,1], (v2[i,2]+1):ncol(ch)] <- 4, next)}
   for (i in 1:nrow(ch)){
      m <- max(which(ch[i,]==1))
      ch[i,f[i]:m] <- 1
      }
   for (i in 1:nrow(v2)){
      u1 <- min(which(ch[v2[i,1],]==1))
      ch[v2[i],u1:(v2[i,2]-1)] <- 1
      }
   for (i in 1:nrow(ch)){
      for (j in f[i]:ncol(ch)){
         if(is.na(ch[i,j])==1) ch[i,j] <- 1
         }
      ch[i,f[i]] <- NA
      }
   return(ch)
   }

inits <- function(){list(ss = runif(1, 0, 1), 
                         ff = runif(1, 0, 1), 
                         pp = runif(1, 0, 1),
                         rr = runif(1, 0, 1), 
                         z = ld.init(rCH, first))}  

# Parameters monitored
parameters <- c("ss", "ff", "rr", "pp")

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 2

# Call JAGS from R
library(jagsUI)
lifedead <- jags(data = jags.data, 
                 inits = inits, 
                 parameters.to.save = parameters, 
                 model.file = "lifedead.jags", 
                 n.chains = nc, 
                 n.thin = nt, 
                 n.iter = ni, 
                 n.burnin = nb, 
                 parallel = TRUE)

print(lifedead, digit = 3)

#traceplot(lifedead)


# JAGS output for model 'lifedead.jags', generated by jagsUI.
# Estimates based on 2 chains of 2000 iterations,
# adaptation = 100 iterations (sufficient),
# burn-in = 1000 iterations and thin rate = 1,
# yielding 2000 total samples from the joint posterior. 
# MCMC ran in parallel for 1.315 minutes at time 2020-11-26 18:52:20.
# 
# mean     sd    2.5%     50%    97.5% overlap0 f  Rhat n.eff
# ss         0.989  0.002   0.985   0.989    0.992    FALSE 1 1.001  2000
# ff         0.466  0.019   0.431   0.466    0.503    FALSE 1 1.015   138
# rr         0.974  0.024   0.912   0.981    0.999    FALSE 1 1.003   583
# pp         0.483  0.030   0.425   0.483    0.542    FALSE 1 1.019   111
# deviance 961.411 44.868 871.082 962.546 1046.672    FALSE 1 1.023    95
# 
# Successful convergence based on Rhat values (all < 1.1). 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 996.3 and DIC = 1957.745 
# DIC is an estimate of expected predictive error (lower is better).






# instead of using the ld.init() function to generate inits for the latent states
# we use the true latent states

zinit <- zsimul
for (i in 1:nrow(zinit)){
   zinit[i, first[i]] <- NA
}

inits <- function(){list(ss = runif(1, 0, 1), 
                         ff = runif(1, 0, 1), 
                         pp = runif(1, 0, 1),
                         rr = runif(1, 0, 1), 
                         z = zinit)}  

# Parameters monitored
parameters <- c("ss", "ff", "rr", "pp")

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 2

# Call JAGS from R
library(jagsUI)
lifedead <- jags(data = jags.data, 
                 inits = inits, 
                 parameters.to.save = parameters, 
                 model.file = "lifedead.jags", 
                 n.chains = nc, 
                 n.thin = nt, 
                 n.iter = ni, 
                 n.burnin = nb, 
                 parallel = TRUE)

print(lifedead, digit = 3)


# 
# 
# JAGS output for model 'lifedead.jags', generated by jagsUI.
# Estimates based on 2 chains of 2000 iterations,
# adaptation = 100 iterations (sufficient),
# burn-in = 1000 iterations and thin rate = 1,
# yielding 2000 total samples from the joint posterior. 
# MCMC ran in parallel for 1.347 minutes at time 2020-11-26 18:54:31.
# 
# mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# ss          0.809  0.008    0.793    0.809    0.825    FALSE 1 1.000  2000
# ff          0.571  0.023    0.525    0.570    0.615    FALSE 1 1.010  2000
# rr          0.089  0.012    0.067    0.089    0.115    FALSE 1 1.006   251
# pp          0.475  0.029    0.417    0.475    0.532    FALSE 1 1.007  2000
# deviance 1291.918 45.178 1206.354 1290.108 1382.250    FALSE 1 1.019  2000
# 
# Successful convergence based on Rhat values (all < 1.1). 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 1020.8 and DIC = 2312.743 
# DIC is an estimate of expected predictive error (lower is better).




# Let's do the frequentist analysis, 
# and play around with the inits to douche check whether we have local minima in the deviance


# Fit multistate model 
# Maximum-likelihood approach
# see Pradel (2005), Gimenez et al. (2012)

# -log(lik) 
devMULTIEVENT <- function(b, data, eff, e, garb, nh, km1){
   
   # data encounter histories, eff counts
   # e vector of dates of first captures
   # garb vector of initial states 
   # km1 nb of recapture occasions (nb of capture occ - 1)
   # nh nb ind
   
   # OBSERVATIONS (+1)
   # 0 neither seen nor recovered
   # 1 seen alive
   # 2 recovered dead
   
   # STATES
   # 1 alive in study area
   # 2 alive outside study area
   # 3 recently dead and recovered
   # 4 recently dead, but not recovered, or dead 
   
   # PARAMETERS
   # ss: true survival probability
   # ff: fidelity probability
   # rr: recovery probability
   # pp: recapture/resighting probability
   
   # logit link for all parameters
   ss <- 1/(1+exp(-b[1]))
   ff <- 1/(1+exp(-b[2]))
   rr <- 1/(1+exp(-b[3]))
   pp <- 1/(1+exp(-b[4]))

   # prob of obs (rows) cond on states (col)
   B <- t(matrix(c(1 - pp, pp, 0,
                 1,  0, 0,
                 1 - rr,   0, rr, 
                  1, 0,  0),
                nrow = 4,
                ncol = 3,
                byrow = T))
   

   # first encounter
   BE <- t(matrix(c(0, 1, 0,
                  1, 0, 0,
                  1, 0, 0,
                  1, 0, 0),
               nrow = 4,
               ncol = 3,
               byrow = T))
   
   # prob of states at t+1 given states at t
   A <- matrix(c(ss * ff, ss * (1 - ff), 1 - ss, 0,
                 0      , ss           , 1 - ss, 0,
                 0      , 0            , 0     , 1,
                 0      , 0            , 0     , 1),
                 nrow = 4,
                 ncol = 4,
                 byrow = T)

   # init states
   PI <- c(1, 0, 0, 0)
   
   # likelihood
   l <- 0
   for (i in 1:nh) # loop on ind
   {
      ei <- e[i] # date of first det
      oe <- garb[i] + 1 # init obs
      evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
      ALPHA <- PI * BE[oe,]
      for (j in (ei+1):(km1+1)) # cond on first capture
      {
         if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!
         ALPHA <- (ALPHA %*% A) * B[evennt[j],]
      }
      l <- l + logprot(sum(ALPHA))*eff[i]
   }
   l <- -l
   l
}


# avoid explosion of log(v) for small values of v
logprot <- function(v){
   eps <- 2.2204e-016
   u <- log(eps) * (1+vector(length=length(v)))
   index <- (v>eps)
   u[index] <- log(v[index])
   u
}

# read in data
CH <- sim$CH
zsimul <- sim$CH.TRUE

# Compute date of first capture
get.first <- function(x) min(which(x!=0))
first <- apply(CH, 1, get.first)

# define various quantities
nh <- dim(CH)[1]
k <- dim(CH)[2]
km1 <- k - 1

# Compute state at first capture init.state
init.state <- NULL
for (i in 1:nh){
   init.state <- c(init.state, CH[i, first[i]])
}


# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
data <- CH  
data[data==3] <- 0

# counts
eff <- rep(1,nh)


# transpose data
data <- t(data)

# b = binit
# data = data
# eff = eff
# e = fc
# garb = init.state
# nh = nh
# km1 = km1
# 
# 
# # evaluate deviance
# devMULTIEVENT(b = binit,
#               data = data,
#               eff = eff,
#               e = fc,
#               garb = init.state,
#               nh = nh,
#               km1 = km1)

# init values

n <- 3
grid <- expand.grid(ss = seq(-5, 5, length = n),
            ff = seq(-5, 5, length = n),
            rr = seq(-5, 5, length = n),
            pp = seq(-5, 5, length = n))
nbMC <- nrow(grid)
res <- matrix(NA, nrow = nbMC, ncol = 5)

for (i in 1:nbMC){

binit <- grid[i,]

# fit model
deb <- Sys.time()
tmpmin <- optim(par = binit,
                fn = devMULTIEVENT,
                gr = NULL,
                hessian = FALSE,
                data,
                eff,
                first,
                init.state,
                nh,
                km1,
                method="BFGS",
                control=list(trace=1, REPORT=1))
fin <- Sys.time()
fin - deb 

# get estimates and back-transform
x <- tmpmin$par
ss <- 1/(1+exp(-x[1]))
ff <- 1/(1+exp(-x[2]))
rr <- 1/(1+exp(-x[3]))
pp <- 1/(1+exp(-x[4]))

res[i,1] <- ss 
res[i,2] <- ff 
res[i,3] <- rr 
res[i,4] <- pp 
res[i,5] <- tmpmin$value
}

res <- round(res,2)
colnames(res) <- c('ss_mle', 'ff_mle', 'rr_mle', 'pp_mle', 'dev')
init <- 1/(1+exp(-grid))
colnames(init) <- c('ss_init', 'ff_init', 'rr_init', 'pp_init')

library(tidyverse)
cbind(res, init) %>% 
   as_tibble() %>%
   arrange(desc(dev)) %>%
   round(2) %>%
   print(n = Inf)


# # A tibble: 81 x 9
# ss_mle ff_mle rr_mle pp_mle   dev ss_init ff_init rr_init pp_init
# <dbl>  <dbl>  <dbl>  <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#    1   1      1      0.95   0.08 2832.    0.01    0.01    0.01    0.99
# 2   1      1      0      0.08 2832.    0.01    0.01    0.5     0.99
# 3   1      1      0      0.08 2832.    0.01    0.01    0.99    0.99
# 4   1      1      0      0.08 2832.    0.5     0.01    0.99    0.99
# 5   1      1      0      0.08 2832.    0.01    0.5     0.99    0.99
# 6   1      0.47   0      0.48 2572.    0.01    0.01    0.99    0.01
# 7   1      0.47   0      0.48 2572.    0.01    0.5     0.99    0.01
# 8   1      0.47   0      0.48 2572.    0.01    0.99    0.99    0.01
# 9   1      0.47   0      0.48 2572.    0.5     0.5     0.99    0.5 
# 10   1      0.47   0      0.48 2572.    0.01    0.99    0.99    0.5 
# 11   1      0.47   0      0.48 2572.    0.01    0.99    0.99    0.99
# 12   0.54   1      0.06   0.36 1195.    0.99    0.01    0.01    0.01
# 13   0.54   1      0.06   0.36 1195.    0.01    0.01    0.99    0.5 
# 14   0.99   0.47   1      0.48 1163.    0.01    0.01    0.01    0.01
# 15   0.99   0.47   1      0.48 1163.    0.01    0.99    0.01    0.01
# 16   0.99   0.47   0.98   0.48 1163.    0.5     0.01    0.99    0.01
# 17   0.99   0.47   0.99   0.48 1163.    0.99    0.01    0.99    0.01
# 18   0.99   0.47   0.99   0.48 1163.    0.99    0.5     0.99    0.01
# 19   0.99   0.47   0.99   0.48 1163.    0.99    0.99    0.99    0.01
# 20   0.99   0.47   0.98   0.48 1163.    0.5     0.01    0.99    0.5 
# 21   0.99   0.47   0.99   0.48 1163.    0.99    0.01    0.99    0.5 
# 22   0.99   0.47   0.99   0.48 1163.    0.99    0.5     0.99    0.5 
# 23   0.99   0.47   1      0.48 1163.    0.99    0.99    0.99    0.5 
# 24   0.99   0.47   0.99   0.48 1163.    0.99    0.01    0.99    0.99
# 25   0.99   0.47   0.99   0.48 1163.    0.99    0.5     0.99    0.99
# 26   0.99   0.47   1      0.48 1163.    0.5     0.99    0.99    0.99
# 27   0.99   0.47   0.99   0.48 1163.    0.99    0.99    0.99    0.99
# 28   0.99   0.47   0.96   0.48 1163.    0.99    0.5     0.01    0.5 
# 29   0.99   0.47   0.94   0.48 1163.    0.5     0.99    0.99    0.01
# 30   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.01    0.01
# 31   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.01    0.01
# 32   0.83   0.55   0.09   0.48 1159.    0.5     0.5     0.01    0.01
# 33   0.83   0.55   0.1    0.48 1159.    0.99    0.5     0.01    0.01
# 34   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.01    0.01
# 35   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.01    0.01
# 36   0.83   0.55   0.1    0.48 1159.    0.01    0.01    0.5     0.01
# 37   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.5     0.01
# 38   0.83   0.55   0.1    0.48 1159.    0.99    0.01    0.5     0.01
# 39   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.5     0.01
# 40   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.5     0.01
# 41   0.83   0.55   0.1    0.48 1159.    0.99    0.5     0.5     0.01
# 42   0.83   0.55   0.09   0.48 1159.    0.01    0.99    0.5     0.01
# 43   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.5     0.01
# 44   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.5     0.01
# 45   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.99    0.01
# 46   0.83   0.55   0.1    0.48 1159.    0.01    0.01    0.01    0.5 
# 47   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.01    0.5 
# 48   0.83   0.55   0.1    0.48 1159.    0.99    0.01    0.01    0.5 
# 49   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.01    0.5 
# 50   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.01    0.5 
# 51   0.83   0.55   0.1    0.48 1159.    0.01    0.99    0.01    0.5 
# 52   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.01    0.5 
# 53   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.01    0.5 
# 54   0.83   0.55   0.1    0.48 1159.    0.01    0.01    0.5     0.5 
# 55   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.5     0.5 
# 56   0.83   0.55   0.1    0.48 1159.    0.99    0.01    0.5     0.5 
# 57   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.5     0.5 
# 58   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.5     0.5 
# 59   0.83   0.55   0.1    0.48 1159.    0.99    0.5     0.5     0.5 
# 60   0.83   0.55   0.1    0.48 1159.    0.01    0.99    0.5     0.5 
# 61   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.5     0.5 
# 62   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.5     0.5 
# 63   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.99    0.5 
# 64   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.99    0.5 
# 65   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.01    0.99
# 66   0.83   0.55   0.1    0.48 1159.    0.99    0.01    0.01    0.99
# 67   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.01    0.99
# 68   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.01    0.99
# 69   0.83   0.55   0.1    0.48 1159.    0.99    0.5     0.01    0.99
# 70   0.83   0.55   0.1    0.48 1159.    0.01    0.99    0.01    0.99
# 71   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.01    0.99
# 72   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.01    0.99
# 73   0.83   0.55   0.1    0.48 1159.    0.5     0.01    0.5     0.99
# 74   0.83   0.55   0.1    0.48 1159.    0.99    0.01    0.5     0.99
# 75   0.83   0.55   0.1    0.48 1159.    0.01    0.5     0.5     0.99
# 76   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.5     0.99
# 77   0.83   0.55   0.1    0.48 1159.    0.99    0.5     0.5     0.99
# 78   0.83   0.55   0.1    0.48 1159.    0.01    0.99    0.5     0.99
# 79   0.83   0.55   0.1    0.48 1159.    0.5     0.99    0.5     0.99
# 80   0.83   0.55   0.1    0.48 1159.    0.99    0.99    0.5     0.99
# 81   0.83   0.55   0.1    0.48 1159.    0.5     0.5     0.99    0.99



# the estimates we get from set of inits 14-29
# look similar to the posterior means we obtain 
# when using the inits from function ld.init()

# the lowest deviance is obtained for MLEs that look
# much alike the posterior means we get from a model
# fitted with the true latent states used as inits


