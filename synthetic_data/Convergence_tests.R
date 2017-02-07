## Using cmdStan - my computer cluster (metworx) seems to have
## a bug with rStan 2.14.

# modelName <- "century_hier_me"
# modelName <- "century_hier_nonstiff"
modelName <- "century_hier"

## Relative paths (adjust for your configuration)
.libPaths("/data/svn-StanPmetrics/script/lib")
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
deliv <- file.path(scriptDir, paste0(modelName, "/deliv"))
toolsDir <- file.path(scriptDir, "tools")
cmdStanDir <- file.path(projectDir, "cmdstan-2.14.0")
modelDir <- file.path(scriptDir, modelName)

source(file.path(toolsDir, "cmdStanTools.R"))
source(file.path(toolsDir, "stanTools.R"))
source("simulate_data_century.R")

## Simulate Data
num_rep <- 5;
t_meas <- c(seq(from=1/360, to=7/360, by=1/360), seq(from=14/360, to=28/360, by=7/360),
            seq(from=60/360, to=360/360, by=30/360));
t_cap <- c(seq(from=0.5/360, to=6.5/360, by=1/360), seq(from=10.5/360, to=24.5/360, by=7/360),
           seq(from=45/360, to=345/360, by=30/360));
init_C <- 1E3*c(.1, .1, .8);
data <- simulate_data_century(t_meas, t_cap, init_C, num_rep);

## Create initial estimate
init <- function(){
    list(turnover = runif(3, 0, 2),
         gamma = rdirichlet(1, c(1, 1, 1))[1, ],
         sigma = runif(3, 0, 2),
         sigma_obs = runif(1, 0, 2),
         A1 = rdirichlet(num_rep, c(1, 1, 1)),
         A2 = rdirichlet(num_rep, c(1, 1, 1)),
         A3 = rdirichlet(num_rep, c(1, 1, 1)),
         A1_g = rdirichlet(1, c(1, 1, 1))[1, ],
         A2_g = rdirichlet(1, c(1, 1, 1))[1, ],
         A3_g = rdirichlet(1, c(1, 1, 1))[1, ],
         kappa = runif(1, 1, 2))
}

## run cmdStan
nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- nPost * nThin
nBurnin <- nBurn * nThin
RNGkind("L'Ecuyer-CMRG")
library(parallel)
mc.reset.stream()

compileModel(model = file.path(scriptDir, modelName), stanDir = cmdStanDir)

library(rstan)
with(data, stan_rdump(ls(data), file = file.path(modelDir, paste0(modelName, ".data.R"))))
inits <- init()
with(inits, stan_rdump(ls(inits), file = file.path(modelDir, paste0(modelName, ".init.R"))))

chains <- 1:nChains
mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init)
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = init, seed = sample(1:999999, 1),
                    chain = chain),
         ##                      adapt_delta = 0.95, stepsize = 0.01),
         model = file.path(modelDir),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = file.path(modelDir, paste0(modelName, ".init.R")), 
         iter = nIter, warmup = nBurn, thin = nThin,
         mc.cores = min(nChains, detectCores()))

glue <- function(...,sep='',collapse=NULL)paste(...,sep=sep,collapse=NULL)
fit <- read_stan_csv(file.path(modelDir, glue(modelName, chains, ".csv")))

######################################################################################
## Diagnostics + Plots
dir.create(deliv)
parametersToPlot <- c("lp__", "gamma", "turnover", "a21", "a31", "a12", "a32", "a13")

graphics.off()
pdf(file = file.path(deliv, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)

pairs(fit, pars = parametersToPlot)
ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(deliv, paste(modelName, "ParameterTable.csv", sep = "")))

dev.off()

## Predictive checks would be nice
