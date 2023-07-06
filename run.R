#Functions to generate data, calculate truncation threshold and calculate posterior for
#effective order
my_path = "." # path to directory containing R code and Stan program,
# starting in current working directory
source(paste(my_path, "/functions.R", sep=""))

#load cmdstanr
library(cmdstanr)
library(stringr)


## Option 1, generate new data set

#example data set
m = 3
p = 2
Sigma = diag(1, m)
mu = rep(0, m)

dataset = generate_data_set(N = 1000, m = m, p = p, Sigma = Sigma, mu = mu)

y = dataset$y
#mean centre data
y = t(apply(y, 1, function(x) x - colMeans(y)))


## Option 2, load one of existing data sets

#example with m = 3, p = 2, for the first experiment
m = 3
p = 2
experiment = 1

dataset = readRDS(paste0(my_path, "/data/Simulation_study_data_m", m, "_p", p, ".Rdata"))

y = dataset$y[[experiment]]
#mean centre data
y = t(apply(y, 1, function(x) x - colMeans(y)))

#correct dimensions when m = 1
if(m == 1) y = t(y)

#compile stan program
file=paste0(my_path, "/multiplicative_gamma.stan")
mod=cmdstan_model(file)

#create named list of everything from data block in stan program
data=list(m=m, p=8, N=1000, y=y, n_miss = 0, ind_miss = c(), df=m+4, S=diag(m), a1=2.5, a2=3, a = 6)

#run HMC with 4 chains on 4 parallel cores, with 1000 iterations of warmup and 4000 sampling iterations
output = mod$sample(data=data, chains=4, parallel_chains=4,
                    iter_warmup=1000,iter_sampling=4000)
#extract MCMC output as 3-dimensional array of form: [iterations, chains, parameters]
draws_arr = output$draws() 
draws_arr = unclass(draws_arr)

#extract samples of P matrices
nchains = 4
P = NULL
for(i in 1:nchains) {
  P = rbind(P, draws_arr[,i,which(stringr::str_detect(dimnames(draws_arr)[[3]], "P\\["))])
}

#calculate truncation threshold
threshold = calculate_threshold(m, N = 1000, beta = 0.99)

#calculate posterior mass for efffective order p*
posterior_mass = eff_order(P = P, pmax = 8, threshold = threshold, singval = FALSE)
posterior_mass$pmf
