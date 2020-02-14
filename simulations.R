#Set-up

require(MASS)
require(ggplot2)

set.seed(147)

###########
# MODEL 1 #
###########

sim1 <- function(assort)
{
  
  output <- NULL ;
  tmp <- list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, assort), c(assort, 1)) #COVARIANCE MATRIX
    mu <- c(0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male) <- c("EnvM", "ChanceM") 
    male$XM <- 0.3 * male$EnvM + 0.7 * rnorm(1000, 0, 1)
    male$YM <- 0.3 * male$XM + 0.3 * male$EnvM + 0.4 * rnorm(1000, 0, 1)
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female) <- c("EnvF", "ChanceF")
    female$XF <- 0.3 * female$EnvF + 0.7 * rnorm(1000, 0, 1)
    female$YF <- 0.3 * female$XF + 0.3 * female$EnvF + 0.4 * rnorm(1000, 0, 1)
    
    om  <- male [order(male$ChanceM), ]
    of  <- female [order(female$ChanceF), ]
    
    pairings <- cbind(om, of)
    
    pairings$Xdiff <- pairings$XM - pairings$XF
    pairings$Ydiff <- pairings$YM - pairings$YF
    pairings$Envdiff <- pairings$EnvM - pairings$EnvF
    
    model1 < -lm(Ydiff ~ Xdiff, data = pairings)
    coef1 <- summary(model1)$coeff[2, 1]
    sd1 <- summary(model1)$coeff[2, 2]
    
    model2 <- lm(Ydiff ~ Xdiff + Envdiff, data = pairings)
    coef2 <- summary(model2)$coeff[2, 1]
    sd2 <- summary(model2)$coeff[2, 2]
    
    tmp[[i]] <- cbind(coef1, sd1, coef2, sd2)
    output <- data.frame(rbind(output, tmp[[i]]))
  }
  names(output) <- c("Beta1", "SE1", "Beta2", "SE2")
  output2 <- NULL;
  output2 <- c(mean(output$Beta1), mean(output$SE1), mean(output$Beta2), mean(output$SE2))
  
  return(output2)
}


###########
# MODEL 2 #
###########

sim2 <- function(assort1, assort2)
{
  
  output <- NULL;
  tmp <- list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, 0, assort1), c(0, 1, assort2), c(assort1, assort2, 1)) #COVARIANCE MATRIX
    mu <- c(0, 0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male) <- c("Trait1M", "Trait2M", "ChanceM") 
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female) <- c("Trait1F", "Trait2F", "ChanceF") 
    
    male$phenM <- male$Trait1M + male$Trait2M + rnorm(1000, 0, 1)
    female$phenF <- female$Trait1F + female$Trait2F + rnorm(1000, 0, 1)
    
    om  <- male [order(male$ChanceM), ]
    of  <- female [order(female$ChanceF), ]
    
    pairings <- cbind(om, of)
    pairings$Trait1Dif <- pairings$Trait1M - pairings$Trait1F
    pairings$Trait2Dif <- pairings$Trait2M - pairings$Trait2F
    pairings$phenDif <- pairings$phenM - pairings$phenF
    
    model1 <- lm(phenDif ~ Trait1Dif, data=pairings)
    coef <- summary(model1)$coeff[2, 1]
    sd <- summary(model1)$coeff[2, 2]
    
    tmp[[i]] <- cbind(coef, sd)
    output <- data.frame(rbind(output, tmp[[i]]))
  }
  names(output) <- c("Beta", "SE")
  output2 <- NULL;
  output2 <- c(mean(output$Beta), mean(output$SE))
  
  return(output2)
}


