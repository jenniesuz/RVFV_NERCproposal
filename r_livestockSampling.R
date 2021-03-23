# script to estimate power to compare RVFV seroprevalence between predicted high and low risk areas
start.time <- Sys.time()

# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)

# no of data sets to simulate per species (1000 takes ~ 4 min)
nsim <- 1000

# sampling / design choices
n.village <- 50
n.hh <- 10         # per village
n <- 8             # no of animals per household

# effect size: OR representing difference in seroprevalence 
# between high and low risk villages
OR <- 3

# load RVFV serology data to get parameter (intercept and variance) estimates
# only need to do this once - then just load the estimates from file - so that
# I don't need to store confidential data locally.

if(!file.exists("parameter.estimates.csv")) {
  # load RVFV serology 
  dat <- read.csv("rvf_livestock_data_2719.csv")[, -1]
  # remove unknown species
  dat <- dat[dat$species != "dk_spe", ]
  # make unique village ID
  dat$village <- paste(dat$ward, dat$village, sep = ".")
  dat$species <- factor(dat$species)
  # loop over species, collecting model estimates
  par.tab <-
    sapply(levels(dat$species), function(sp) {
      datsp <- dat[dat$species == sp, ]
      fit <- glmer(result ~ (1 | barcode_hh) +(1 | village), family = binomial, data = datsp)
      fixef(fit)
      round(c(mean.hh.n = mean(table(datsp$barcode_hh)), 
              fixef(fit), unlist(VarCorr(fit))), 2)
    })
  # write results to file
  write.csv(par.tab, file = "parameter.estimates.csv")
}

# load parameter estimates
par.tab <- read.csv("parameter.estimates.csv", row.names = 1)
species <- colnames(par.tab)

# simulate RVFV serology data in each species

# function to simulate data and estimate p-value for null hypothesis 
# that high and low risk areas have the same seroprevalence
res.tab.fn <- function(...) {
  sapply(species, function(sp) {
    
    # create template data set
    dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
    print(sum(dat$n))
    # allocate villages to high and low prevalence in 1:1 ratio 
    dat$risk.level <- dat$village %% 2 - 0.5
    # simulate seropositives
    simdat <-
      sim.glmm(
        design.data = dat, 
        fixed.eff = 
          list(
            intercept = par.tab["(Intercept)", sp],
            risk.level = log(OR)),
        distribution = "binomial",
        rand.V = c(hh = par.tab["barcode_hh", sp], 
                   village = par.tab["village", sp]))
    
    fit <- glmer(cbind(response, n - response) ~ risk.level + (1 | hh) +(1 | village), family = binomial, data = simdat)
    fit0 <- update(fit, ~ . - risk.level)
    #coef(summary(fit))["risk.level", "Pr(>|z|)"] # Wald P not reliable - gives inflated type 1 error
    anova(fit, fit0)[2, "Pr(>Chisq)"]
  })
}


# repeat simulations many times and calculate p-value
sim.res <- mclapply(1:nsim, res.tab.fn, mc.cores = detectCores())
print(Sys.time() - start.time)

# estimate power
apply(do.call("rbind", sim.res) < 0.05, 2, mean)


