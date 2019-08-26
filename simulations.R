require(MASS)
require(ggplot2)
set.seed(147)

###########
# MODEL 1 #
###########

sim1<-function(assort)
{
  
  output<-NULL;
  tmp<-list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, assort), c(assort, 1)) #COVARIANCE MATRIX
    mu <- c(0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male)<-c("EnvM", "ChanceM") 
    male$XM<-0.3*male$EnvM+0.7*rnorm(1000,0,1)
    male$YM<-0.3*male$XM+0.3*male$EnvM+0.4*rnorm(1000,0,1)
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female)<-c("EnvF", "ChanceF")
    female$XF<-0.3*female$EnvF+0.7*rnorm(1000,0,1)
    female$YF<-0.3*female$XF+0.3*female$EnvF+0.4*rnorm(1000,0,1)
    
    om  <- male [order(male $ChanceM),]
    of  <- female [order(female $ChanceF),]
    
    pairings <- cbind(om, of)
    
    pairings$Xdiff<-pairings$XM-pairings$XF
    pairings$Ydiff<-pairings$YM-pairings$YF
    pairings$Envdiff<-pairings$EnvM-pairings$EnvF
    
    model1<-lm(Ydiff ~ Xdiff,data=pairings)
    coef1<-summary(model1)$coeff[2,1]
    sd1<-summary(model1)$coeff[2,2]
    
    model2<-lm(Ydiff ~ Xdiff+Envdiff,data=pairings)
    coef2<-summary(model2)$coeff[2,1]
    sd2<-summary(model2)$coeff[2,2]
    
    tmp[[i]]<-cbind(coef1, sd1, coef2, sd2)
    output<-data.frame(rbind(output, tmp[[i]]))
  }
  names(output)<-c("Beta1", "SE1", "Beta2", "SE2")
  output2<-NULL;
  output2<-c(mean(output$Beta1), mean(output$SE1), mean(output$Beta2), mean(output$SE2))
  
  return(output2)
}

t1<-sim1(0)[1]
t2<-sim1(0.1)[1]
t3<-sim1(0.2)[1]
t4<-sim1(0.3)[1]
t5<-sim1(0.4)[1]
t6<-sim1(0.5)[1]
t7<-sim1(0.6)[1]
t8<-sim1(0.7)[1]
t9<-sim1(0.8)[1]
t10<-sim1(0.9)[1]



#Model 1 Plot


x<-c(0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0.6, 0.7, 0.8, 0.9, 1)

t11=0.3
bias<-c(t1, t2, t3, t4, t5, t6,
        t7, t8, t9, t10, t11)


df<-data.frame(x,bias)

png(filename="O:/Postdoc/Between spouse/Figures/Figure3A.png", width=20, height=10, units="cm", res=600)
c1<-ggplot(data=df,
           aes(x, bias)) +
  geom_line(color="darkblue") +
  xlim(0,1) +
  ylim(0.3, 0.5)+
  ggtitle('Figure 3A') +
  labs(subtitle="Simulations for Model 1: Spousal correlations controlling for confounding",
       y="Effect estimate of X on Y",
       x="Spousal correlation for E") +
  expand_limits(y=0.29)+ 
  geom_text(aes(0,0.3, label = "Unbiased estimate", vjust = -1, hjust=0.2))

c1 + geom_hline(yintercept=0.3, color="darkred")+geom_text(aes(0,0.3, label = "Between-spouse estimate", vjust = -20, hjust=0.13))
dev.off()


###########
# MODEL 2 #
###########

sim2<-function(assort1, assort2)
{
  
  output<-NULL;
  tmp<-list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, 0, assort1), c(0, 1, assort2), c(assort1, assort2, 1)) #COVARIANCE MATRIX
    mu <- c(0, 0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male)<-c("Trait1M", "Trait2M", "ChanceM") 
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female)<-c("Trait1F", "Trait2F", "ChanceF") 
    
    male$phenM<-male$Trait1M+male$Trait2M
    female$phenF<-female$Trait1F+female$Trait2F
    
    om  <- male [order(male $ChanceM),]
    of  <- female [order(female $ChanceF),]
    
    pairings <- cbind(om, of)
    pairings$Trait1Dif<-pairings$Trait1M-pairings$Trait1F
    pairings$Trait2Dif<-pairings$Trait2M-pairings$Trait2F
    pairings$phenDif<-pairings$phenM-pairings$phenF
    
    model1<-lm(phenDif ~ Trait1Dif,data=pairings)
    coef<-summary(model1)$coeff[2,1]
    sd<-summary(model1)$coeff[2,2]
    
    tmp[[i]]<-cbind(coef, sd)
    output<-data.frame(rbind(output, tmp[[i]]))
  }
  names(output)<-c("Beta", "SE")
  output2<-NULL;
  output2<-c(mean(output$Beta), mean(output$SE))
  
  return(output2)
}

a1<-(sim2(0.1,0)/1)[1]
a2<-(sim2(0.1,0.1)/1)[1]
a3<-(sim2(0.1,0.2)/1)[1]
a4<-(sim2(0.1,0.3)/1)[1]
a5<-(sim2(0.1,0.4)/1)[1]
a6<-(sim2(0.1,0.5)/1)[1]

b1<-(sim2(0.2,0)/1)[1]
b2<-(sim2(0.2,0.1)/1)[1]
b3<-(sim2(0.2,0.2)/1)[1]
b4<-(sim2(0.2,0.3)/1)[1]
b5<-(sim2(0.2,0.4)/1)[1]
b6<-(sim2(0.2,0.5)/1)[1]

c1<-(sim2(0.3,0)/1)[1]
c2<-(sim2(0.3,0.1)/1)[1]
c3<-(sim2(0.3,0.2)/1)[1]
c4<-(sim2(0.3,0.3)/1)[1]
c5<-(sim2(0.3,0.4)/1)[1]
c6<-(sim2(0.3,0.5)/1)[1]

d1<-(sim2(0.4,0)/1)[1]
d2<-(sim2(0.4,0.1)/1)[1]
d3<-(sim2(0.4,0.2)/1)[1]
d4<-(sim2(0.4,0.3)/1)[1]
d5<-(sim2(0.4,0.4)/1)[1]
d6<-(sim2(0.4,0.5)/1)[1]

e1<-(sim2(0.5,0)/1)[1]
e2<-(sim2(0.5,0.1)/1)[1]
e3<-(sim2(0.5,0.2)/1)[1]
e4<-(sim2(0.5,0.3)/1)[1]
e5<-(sim2(0.5,0.4)/1)[1]
e6<-(sim2(0.5,0.5)/1)[1]

#Model 2 Plot
library(ggplot2)
x<-as.factor(c("0.1", "0.1", "0.1", "0.1", "0.1", "0.1", 
               "0.2", "0.2", "0.2", "0.2", "0.2", "0.2", 
               "0.3", "0.3", "0.3", "0.3", "0.3", "0.3", 
               "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", 
               "0.5", "0.5", "0.5", "0.5", "0.5", "0.5" ))

y<-c(0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5)

bias<-c(a1, a2, a3, a4, a5, a6,
        b1, b2, b3, b4, b5, b6,
        c1, c2, c3, c4, c5, c6,
        d1, d2, d3, d4, d5, d6,
        e1, e2, e3, e4, e5, e6)

df<-data.frame(x,y,bias)

png(filename="O:/Postdoc/Between spouse/Figures/Figure3B.png", width=20, height=10, units="cm", res=600)
c2<-ggplot(data=df,
           aes(y, bias, colour=x)) +
  geom_line() +
  ggtitle('Figure 3B') +
  labs(subtitle="Simulations for Model 2: Between-spouse: assortment and collider bias",
       y="Between-spouse estimate of X1 on Y",
       x="Degree of assortment on X2",
       colour="Degree of assortment on X1") 
c2
dev.off()

###########
# MODEL 3 #
###########

sim3<-function(assort, h2)
{
  
  output<-NULL;
  tmp<-list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, sqrt(h2)), c(sqrt(h2), 1)) #COVARIANCE MATRIX
    mu <- c(0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male)<-c("genM", "phenM") 
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female)<-c("genF", "phenF") 
    
    #Fix assortment
    
    sqrtassort  <- sqrt(assort)                   
    theta <- acos(sqrtassort)             
    temp1     <- rnorm(1000, 0, 1)      
    temp2     <- rnorm(1000, 0, 1)
    
    X1     <- cbind(male$phenM  , temp1 )  
    X2    <- cbind(female$phenF  , temp2 )       
    Xctr1   <- scale(X1 , center=TRUE, scale=FALSE)   
    Xctr2   <- scale(X2 , center=TRUE, scale=FALSE)
    
    Id   <- diag(1000)                               
    Q1     <- qr.Q(qr(Xctr1 [ , 1, drop=FALSE]))
    Q2     <- qr.Q(qr(Xctr2 [ , 1, drop=FALSE]))        
    P1     <- tcrossprod(Q1 )          
    P2     <- tcrossprod(Q2 )       
    x2o1   <- (Id-P1 ) %*% Xctr1 [ , 2]   
    x2o2   <- (Id-P2 ) %*% Xctr2 [ , 2]              
    Xc21   <- cbind(Xctr1 [ , 1], x2o1 )  
    Xc22   <- cbind(Xctr2 [ , 1], x2o2 )               
    Y1     <- Xc21  %*% diag(1/sqrt(colSums(Xc21 ^2)))  
    Y2     <- Xc22  %*% diag(1/sqrt(colSums(Xc22 ^2))) 
    
    male $ChanceM <- Y1 [ , 2] + (1 / tan(theta)) * Y1 [ , 1]     
    female $ChanceF <- Y2 [ , 2] + (1 / tan(theta)) * Y2 [ , 1] 
    om  <- male [order(male $ChanceM),]
    of  <- female [order(female $ChanceF),]
    
    pairings <- cbind(om, of)
    pairings$GenDif<-pairings$genM-pairings$genF
    pairings$PhenDif<-pairings$phenM-pairings$phenF
    
    model1<-lm(PhenDif ~ GenDif,data=pairings)
    coef<-summary(model1)$coeff[2,1]
    sd<-summary(model1)$coeff[2,2]
    
    tmp[[i]]<-cbind(coef, sd)
    output<-data.frame(rbind(output, tmp[[i]]))
  }
  names(output)<-c("Beta", "SE")
  output2<-NULL;
  output2<-c(mean(output$Beta), mean(output$SE))
  
  return(output2)
}

l1<-sim3(0,0.2)
a1<-sim3(0.1, 0.2)[1]/l1
a2<-sim3(0.2, 0.2)[1]/l1
a3<-sim3(0.3, 0.2)[1]/l1
a4<-sim3(0.4, 0.2)[1]/l1
a5<-sim3(0.5, 0.2)[1]/l1
a6<-sim3(0.6, 0.2)[1]/l1
a7<-sim3(0.7, 0.2)[1]/l1
a8<-sim3(0.8, 0.2)[1]/l1

l2<-sim3(0, 0.4)
b1<-sim3(0.1, 0.4)[1]/l2
b2<-sim3(0.2, 0.4)[1]/l2
b3<-sim3(0.3, 0.4)[1]/l2
b4<-sim3(0.4, 0.4)[1]/l2
b5<-sim3(0.5, 0.4)[1]/l2
b6<-sim3(0.6, 0.4)[1]/l2
b7<-sim3(0.7, 0.4)[1]/l2
b8<-sim3(0.8, 0.4)[1]/l2

l3<-sim3(0, 0.6)
c1<-sim3(0.1, 0.6)[1]/l3
c2<-sim3(0.2, 0.6)[1]/l3
c3<-sim3(0.3, 0.6)[1]/l3
c4<-sim3(0.4, 0.6)[1]/l3
c5<-sim3(0.5, 0.6)[1]/l3
c6<-sim3(0.6, 0.6)[1]/l3
c7<-sim3(0.7, 0.6)[1]/l3
c8<-sim3(0.8, 0.6)[1]/l3

l4<-sim3(0, 0.8)
d1<-sim3(0.1, 0.8)[1]/l4
d2<-sim3(0.2, 0.8)[1]/l4
d3<-sim3(0.3, 0.8)[1]/l4
d4<-sim3(0.4, 0.8)[1]/l4
d5<-sim3(0.5, 0.8)[1]/l4
d6<-sim3(0.6, 0.8)[1]/l4
d7<-sim3(0.7, 0.8)[1]/l4
d8<-sim3(0.8, 0.8)[1]/l4


#Model 3 Plot

h2<-as.factor(c("0.2", "0.2", "0.2", "0.2", "0.2", "0.2", "0.2", "0.2", 
                "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", 
                "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", 
                "0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8" ))

p<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

bias<-c(a1[1],a2[1],a3[1],a4[1],a5[1],a6[1],a7[1],a8[1],
        b1[1],b2[1],b3[1],b4[1],b5[1],b6[1],b7[1],b8[1],
        c1[1],c2[1],c3[1],c4[1],c5[1],c6[1],c7[1],c8[1],
        d1[1],d2[1],d3[1],d4[1],d5[1],d6[1],d7[1],d8[1])


df2<-data.frame(h2,p,bias)

png(filename="O:/Postdoc/Between spouse/Figures/Figure3C.png", width=20, height=10, units="cm", res=600)
c3<-ggplot(data=df2,
           aes(p, bias, colour=h2)) +
  geom_line() +
  ggtitle('Figure 3C') +
  labs(subtitle="Simulations for Model 3: Between-spouse: genetic associations and collider bias",
       y="Between-spouse estimate of G on X",
       x="Degree of assortment on X",
       colour="Heritability of X") 

c3
dev.off()

###########
# MODEL 4 #
###########



sim4<-function(assort1, assort2)
{
  
  output<-NULL;
  tmp<-list()
  
  for (i in c(1:1000)) {
    
    rho <- cbind(c(1, 0, assort1), c(0, 1, assort2), c(assort1, assort2, 1)) #COVARIANCE MATRIX
    mu <- c(0, 0, 0) #STANDARDISED MEANS
    
    male <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(male)<-c("Trait1M", "Trait2M", "ChanceM") 
    
    female <- data.frame(mvrnorm(1000, mu=mu, Sigma=rho))
    names(female)<-c("Trait1F", "Trait2F", "ChanceF") 
    
    male$phenM<-male$Trait1M+male$Trait2M
    female$phenF<-female$Trait1F+female$Trait2F
    
    om  <- male [order(male $ChanceM),]
    of  <- female [order(female $ChanceF),]
    
    pairings <- cbind(om, of)
    pairings$Trait1Dif<-pairings$Trait1M-pairings$Trait1F
    pairings$Trait2Dif<-pairings$Trait2M-pairings$Trait2F
    pairings$phenDif<-pairings$phenM-pairings$phenF
    
    pairings$Child<-0.2*pairings$Trait1F+0.2*pairings$Trait2F+0.6*rnorm(1000,0, 1)
    model1<-lm(Child ~ Trait1M,data=pairings)
    coef<-summary(model1)$coeff[2,1]
    sd<-summary(model1)$coeff[2,2]
    
    tmp[[i]]<-cbind(coef, sd)
    output<-data.frame(rbind(output, tmp[[i]]))
  }
  names(output)<-c("Beta", "SE")
  output2<-NULL;
  output2<-c(mean(output$Beta), mean(output$SE))
  
  return(output2)
}

##Save coefficients of Model 4


a1<-(sim4(0.1,0)/0.2)[1]
a2<-(sim4(0.1,0.1)/0.2)[1]
a3<-(sim4(0.1,0.2)/0.2)[1]
a4<-(sim4(0.1,0.3)/0.2)[1]
a5<-(sim4(0.1,0.4)/0.2)[1]
a6<-(sim4(0.1,0.5)/0.2)[1]

b1<-(sim4(0.2,0)/0.2)[1]
b2<-(sim4(0.2,0.1)/0.2)[1]
b3<-(sim4(0.2,0.2)/0.2)[1]
b4<-(sim4(0.2,0.3)/0.2)[1]
b5<-(sim4(0.2,0.4)/0.2)[1]
b6<-(sim4(0.2,0.5)/0.2)[1]

c1<-(sim4(0.3,0)/0.2)[1]
c2<-(sim4(0.3,0.1)/0.2)[1]
c3<-(sim4(0.3,0.2)/0.2)[1]
c4<-(sim4(0.3,0.3)/0.2)[1]
c5<-(sim4(0.3,0.4)/0.2)[1]
c6<-(sim4(0.3,0.5)/0.2)[1]

d1<-(sim4(0.4,0)/0.2)[1]
d2<-(sim4(0.4,0.1)/0.2)[1]
d3<-(sim4(0.4,0.2)/0.2)[1]
d4<-(sim4(0.4,0.3)/0.2)[1]
d5<-(sim4(0.4,0.4)/0.2)[1]
d6<-(sim4(0.4,0.5)/0.2)[1]

e1<-(sim4(0.5,0)/0.2)[1]
e2<-(sim4(0.5,0.1)/0.2)[1]
e3<-(sim4(0.5,0.2)/0.2)[1]
e4<-(sim4(0.5,0.3)/0.2)[1]
e5<-(sim4(0.5,0.4)/0.2)[1]
e6<-(sim4(0.5,0.5)/0.2)[1]

#Model 4 Plot

x<-as.factor(c("0.1", "0.1", "0.1", "0.1", "0.1", "0.1", 
               "0.2", "0.2", "0.2", "0.2", "0.2", "0.2", 
               "0.3", "0.3", "0.3", "0.3", "0.3", "0.3", 
               "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", 
               "0.5", "0.5", "0.5", "0.5", "0.5", "0.5" ))

y<-c(0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5,
     0, 0.1, 0.2, 0.3, 0.4, 0.5)

bias<-c(a1, a2, a3, a4, a5, a6,
        b1, b2, b3, b4, b5, b6,
        c1, c2, c3, c4, c5, c6,
        d1, d2, d3, d4, d5, d6,
        e1, e2, e3, e4, e5, e6)

df<-data.frame(x,y,bias)

png(filename="O:/Postdoc/Between spouse/Figures/Figure3D.png", width=20, height=10, units="cm", res=600)
c4<-ggplot(data=df,
           aes(y, bias, colour=x)) +
  geom_line() +
  ggtitle('Figure 3D') +
  labs(subtitle="Simulations for Model 4: Assortment and parental comparisons",
       y="Bias in estimate of Paternal X1 on Y",
       x="Degree of assortment on X2",
       colour="Degree of assortment on X1") 

c4
dev.off()
