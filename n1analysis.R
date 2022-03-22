## Bayesian maladaptation analysis for N1
rm(list=ls())

library(R2jags)

# load data
dd <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data_16Sept14.csv")
names(dd)
# explore
dd$biomass
dd$a_chew_herb
hist(dd$a_chew_herb)
hist(dd$no_art_5)
hist(dd$biomass)
hist(dd$mal)
hist(dd$cn_ratio)
hist(dd$ln.vol_cub)
# housekeeping -- get mal and art only
bd <- data.frame(mal=dd$mal, art=dd$no_art_5, bio=dd$biomass, size=dd$ln.vol_cub, 
                 cn=dd$cn_ratio, host=ifelse(dd$host=="A", 0, 1), con=dd$con_art5)
bd <- bd[!is.na(bd$mal),]

# standardize covariates
bd$size.std <- as.numeric(scale(bd$size))
bd$mal.std <- as.numeric(scale(bd$mal))
bd$con.std <- as.numeric(scale(bd$con))

bd$con.host <- bd$con.std * bd$host

# frequentist approach

an1 <- lm(cn ~ size.std, data=bd)
summary(an1)

an2 <- glm(art ~ mal.std + size.std + host*con.std, family=quasipoisson, data=bd)
summary(an2)
# Bayesian

JAGS.data  <- list(Y      = bd$art,  #Response  
                   X      = X,               #Intercept + covariates
                   K      = K,               #Number of betas       
                   N      = nrow(bd))  # sample size


sink("n1_mal_art_pois.txt")
cat("
  model{
# priors

b0 ~ dnorm(0, 0.0001)
b.mal ~ dnorm(0, 0.0001)
b.host ~ dnorm(0, 0.0001)

# likelihood

for (i in 1:N) {
art[i] ~ dpois(mu[i])
log(mu[i]) <- b0 + b.mal * mal[i] + b.host * host[i]

}
# residuals

for (i in 1:N) {
E[i] <- art[i] - mu[i]
}

  }
    ", fill=TRUE)
sink()

# initial values

init <- function() {
  list(b0=rnorm(1, 0, .5), b.mal=rnorm(1, 0, .5), b.host=rnorm(1, 0, 0.5))
}

params <- c("b0", "b.mal", "b.host")

j1 <- jags(data=jags_data,
           inits=init,
           parameters=params,
           model="n1_mal_art_pois.txt", 
           n.thin=10, 
           n.chains=3, 
           n.burnin=4000,
           n.iter=5000)

j1 <- update(j1, n.iter=50000, n.thin=10)
j1


### negative binomial bayesian model

X <- model.matrix(~ 1 + mal.std + host + con.std +
                    con.std : host, 
                  data = bd)
K <- ncol(X) #Number of regression parameters

JAGS.data  <- list(Y      = bd$art,    #Response  
                   X      = X,         #Intercept + covariates
                   K      = K,         #Number of betas       
                   N      = nrow(bd))

sink("n1_mal_art_nb.txt")
cat("
    model{
    
    ## priors
    
    for (i in 1:K) {beta[i] ~ dnorm(0, 0.001)} # betas
    
    # half-Cauchy priors for size parameter

    num.s   ~ dnorm(0, 0.0016)            
    denom.s ~ dnorm(0, 1)                 
    size <- abs(num.s / denom.s) 
    
    ## likelihood
    
    for (i in 1:N) {
    art[i] ~ dnegbin(p[i], size)
    p[i] <- size / (size + mu[i]) 
    log(mu[i]) <- inprod(beta[], X[i,])

    ## residuals and simulated data
    
    Exp[i] <- mu[i] 
    Var[i] <- mu[i] + mu[i]^2 / size
    
    # Pearson residuals
    E[i]   <- (art[i]  - Exp[i]) / sqrt(Var[i])    
    
    # Simulated NB data 
    YNew[i] ~  dnegbin(p[i], size) 
    
    # Pearson residuals for simulated data
    ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
    D[i]    <- pow(E[i], 2)                      
    DNew[i] <- pow(ENew[i], 2)   
    } 

    Fit         <- sum(D[1:N])      #Sum of squared residuals  
    FitNew      <- sum(DNew[1:N])  #Sum of squared predicted residuals
    
}
    ", fill=TRUE)
sink()

cat("
    model{
    #1A. Priors beta and sigma
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    
    #1C. Prior for size
    num.s   ~ dnorm(0, 0.0016)            #<----from squirrel exercise
    denom.s ~ dnorm(0, 1)                 #<----from squirrel exercise
    size <- abs(num.s / denom.s)          #<----from squirrel exercise
    
    
    #2. Likelihood
    for (i in 1:N) {
    Y[i] ~  dnegbin(p[i], size)
    p[i] <- size / (size + mu[i])  
    
    log(mu[i]) <- eta[i] 
    eta[i]     <- inprod(beta[], X[i,]) 
    
    #3. New stuff: Discrepancy measures 
    Exp[i] <- mu[i] 
    Var[i] <- mu[i] + mu[i] * mu[i] / size
    
    #Pearson residuals
    E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
    
    #Simulated NB data 
    YNew[i] ~  dnegbin(p[i], size) 
    
    #Pearson residuals for simulated data
    ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
    D[i]    <- pow(E[i], 2)                      
    DNew[i] <- pow(ENew[i], 2)   
    }          
    Fit         <- sum(D[1:N])      #Sum of squared residuals  
    FitNew      <- sum(DNew[1:N])  #Sum of squared predicted residuals
    }
    ",fill = TRUE)
sink()

# initial values
inits  <- function () {
  list(
    beta  = rnorm(K, coef(an2), 0.1),
    num   = rnorm(1, 0, 25),    #<-- half-Cauchy(25)
    denom = rnorm(1, 0, 1),     #<-- half-Cauchy(25)
    num.s   = rnorm(1, 0, 25),  #for size
    denom.s = rnorm(1, 0, 1))   #for size
}

params <- c("beta", "E", 
            "mu", 
            "Fit", "FitNew",
            "size")

params <- c("beta")

j2 <- jags(data=JAGS.data,
           inits=inits,
           parameters=params,
           model="n1_mal_art_nb.txt", 
           n.thin=10, 
           n.chains=3, 
           n.burnin=4000,
           n.iter=5000)

j2 <- update(j2, n.iter=50000, n.thin=10)
j2

#5. Assess mixing
MyNames <- c(colnames(X), "sigma_Nest", "size")

MyBUGSChains(out, 
             c(uNames("beta", K), "sigma_Nest", "size"),
             PanelNames = MyNames)
MyBUGSACF(out, 
          c(uNames("beta", K), "sigma_Nest", "size"),
          PanelNames = MyNames)

#6. Present output
OUT1 <- MyBUGSOutput(out, 
                     c(uNames("beta", K), "sigma_Nest", "size"),
                     VarNames = MyNames)
print(OUT1, digits =5)

MyBUGSHist(out, 
           c(uNames("beta", K), "sigma_Nest", "size"),
           PanelNames = MyNames)

           