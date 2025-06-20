# ---------------------------------------------------------------
# Author: Rafael Prieto-Curiel
# Date: May 2025
# Title: Cartel Participation Model – Individual and Collective Simulation
#
# Description:
# This R script simulates a cartel participation model at both the 
# individual and collective levels. It includes routines to generate 
# figures that illustrate the dynamics and outcomes of the simulation.
# It also captures the optimal size of cartel depending on the personal costs
# ---------------------------------------------------------------

#### logic
# M people
# binary variable for A being part of a cartel

#### parameters
{
require(scales)
N <- 25000 #### number of initial people
sat <- 50/(N^2) #### saturation parameter 
rho = 150/25000 #### rate of recruitment each week with base = 1/2
pi_v = 0.8 ### probability of voluntary recruitment
pi_f = 0.1 ### probability of forced recruitment
pi_c = 1 - pi_v - pi_f ### probability of being ignored during recruitment
B = 1/1000 #### sensitivity to money, to get rid of this parameter
delta = 5000 #### sensitivity to incapacitation
eta = 10000 #### sensitivity to lethality
#### recruitment rate 0.0026
#### max cartel size 50000
### yearly data in US$
w = 14.7*300 ### minimum yearly salary working 300 days in a year
v = 832 ### social program
r = 250*52 #### U$13000 yearly cartel salary
omega <- 0.2*r
kappaP = 800/(52*25000) #### probability of arrest
kappaK = 700/(52*25000) #### probability of killed 
alpha = -log(pi_v/(0.8-pi_f) - 1) - B* (r - w - v) +delta*kappaP + eta*kappaK#### baseline for the function of recruitment
}

#### functions
{
IsRecruited <- function(N, alpha, B, delta, eta, w, v,r, kappaP, kappaK, pi_v, pi_f){
  #### a function that takes all values and returns a 1 if the person is recruited
  s <- alpha  + B*(r - w - v) - delta * kappaP - eta * kappaK
  return(sum(runif(N) < (exp(s)/(1+exp(s))*pi_v+pi_f)))
}
  
PlotIsRecruited <- function(alpha, B, delta, eta, w, v,r, kappaP, kappaK, pi_v, pi_f){
      #### a function that takes all values and returns a 1 if the person is recruited
      s <- alpha  + B*(r - w - v) - delta * kappaP - eta * kappaK
      return((exp(s)/(1+exp(s))*pi_v+pi_f))
} 

#### run simulations and return the whole time series
EffectTimeSeries <- function(N, alpha, B, delta, eta, w, 
                     v, #### social program
                     r, #### cartel salary
                     kappaP, kappaK, pi_v, pi_f, rho){
  steps <- 1500 ### so 30 years
  
  #### Initial members of a cartel  
  Active <- N
  
  #### not recruited
  NotRecruited <- 0*Active
  Killed <- 0*Active #### killed
  Incapacitated <- 0*Active
  Recovered <- Killed + Incapacitated
  Exposed <- 0
  Saturated <- 0
  
  #### time vars
  At <- Active
  Kt <- Killed
  It <- Incapacitated
  Rt <- Recovered
  Nt <- NotRecruited
  Et <- Exposed
  St <- Saturated
  
  #### is there people to recruit?
  for(k in 1:steps){
    Exposed <- round(Active* rho)
    I <- IsRecruited(N = Exposed,
                     alpha = alpha, 
                     B = B, 
                     delta = delta, 
                     eta = eta, 
                     w = w, 
                     v = v,
                     r = r, 
                     kappaP = kappaP, 
                     kappaK = kappaK, 
                     pi_v = pi_v, 
                     pi_f = pi_f)
    NotRecruited <- Exposed - I ##### they were exposed and were not recruited
    Active <- Active + I ##### new members of a cartel
    
    #### suimualte arrest
    Incapacitated <- sum(runif(Active) < kappaP)
    Active <- Active - Incapacitated
    
    #### suimualte kills
    Killed <- sum(runif(Active) < kappaK)
    Active <- Active - Killed
    Recovered <- Killed + Incapacitated
    
    #### simulate saturation
    Saturated <- round(sat*Active^2)
    Active <- Active - Saturated
    
    At <- c(At, Active)
    Kt <- c(Kt, Killed)
    It <- c(It, Incapacitated)
    Rt <- c(Rt, Recovered)
    Nt <- c(Nt, NotRecruited)
    Et <- c(Et, Exposed)
    St <- c(St, Saturated)
  }
  return(data.frame(At = At, 
                    Kt = Kt, 
                    It = It, 
                    Rt = Rt, 
                    Nt = Nt, 
                    Et = Et, 
                    St = St))
}

#### run simulations and return the whole time series
ImpactSP <- function(N, alpha, B, delta, eta, w, 
                             v, #### social program
                             r, #### cartel salary
                             kappaP, kappaK, pi_v, pi_f, rho){
  steps <- 500 ### so 10 years
  
  #### Initial members of a cartel  
  Active <- N
  
  #### not recruited
  NotRecruited <- 0*Active
  Killed <- 0*Active #### killed
  Incapacitated <- 0*Active
  Recovered <- Killed + Incapacitated
  Exposed <- 0
  Saturated <- 0
  
  #### time vars
  At <- Active
  Kt <- Killed
  It <- Incapacitated
  Rt <- Recovered
  Nt <- NotRecruited
  Et <- Exposed
  St <- Saturated
  
  #### is there people to recruit?
  for(k in 1:steps){
    Exposed <- round(Active* rho)
    I <- IsRecruited(N = Exposed,
                     alpha = alpha, 
                     B = B, 
                     delta = delta, 
                     eta = eta, 
                     w = w, 
                     v = v,
                     r = r, 
                     kappaP = kappaP, 
                     kappaK = kappaK, 
                     pi_v = pi_v, 
                     pi_f = pi_f)
    NotRecruited <- Exposed - I ##### they were exposed and were not recruited
    Active <- Active + I ##### new members of a cartel
    
    #### suimualte arrest
    Incapacitated <- sum(runif(Active) < kappaP)
    Active <- Active - Incapacitated
    
    #### suimualte kills
    Killed <- sum(runif(Active) < kappaK)
    Active <- Active - Killed
    Recovered <- Killed + Incapacitated
    
    #### simulate saturation
    Saturated <- round(sat*Active^2)
    Active <- Active - Saturated
    
    At <- c(At, Active)
    Kt <- c(Kt, Killed)
    It <- c(It, Incapacitated)
    Rt <- c(Rt, Recovered)
    Nt <- c(Nt, NotRecruited)
    Et <- c(Et, Exposed)
    St <- c(St, Saturated)
  }
  
  return(data.frame(Active = Active, 
                    Killed = sum(Kt), 
                    Incapacitated = sum(It), 
                    Recovered = sum(Rt), 
                    NotRecruited = sum(Nt), 
                    Exposed = sum(Et), 
                    Saturated = sum(St)))
}

#### smoothing curve with intervals
IV <- function(x, y, span = 0.3){
  M <- loess.smooth(x, y,
                    evaluation = 1000,
                    span = span, family = "gaussian")
  Res <- M$y[floor(1+999*x)]
#  MRes <- loess.smooth(x,abs(Res- y),
#                       evaluation = 1000,
#                       span = .3, family = "gaussian")
  return(data.frame(x = M$x,
                    mean = M$y))
}


}

#### check parameters and how to call functions
{
IsRecruited(N = 1000000,
            alpha = alpha, 
            B = B, 
            delta = delta, 
            eta = eta, 
            w = w, 
            v = v,
            r = r, 
            kappaP = kappaP, 
            kappaK = kappaK, 
            pi_v = pi_v, 
            pi_f = pi_f)

Efs <- EffectTimeSeries(N = 25000,
                        alpha = alpha, 
                        B = B, 
                        delta = delta, 
                        eta = eta, 
                        w = w, 
                        v = v, #### social program with 832 as baseline
                        r = r*1, #### cartel salary 
                        kappaP = kappaP, 
                        kappaK = kappaK, 
                        pi_v = pi_v, 
                        pi_f = pi_f,
                        rho = rho)

Res <- ImpactSP(N = 25000,
         alpha = alpha, 
         B = B, 
         delta = delta, 
         eta = eta, 
         w = w, 
         v = v, #### social program with 832 as baseline
         r = r*1, #### cartel salary 
         kappaP = kappaP, 
         kappaK = kappaK, 
         pi_v = pi_v, 
         pi_f = pi_f,
         rho = rho)

#### run for a huge social program
IRes <- ImpactSP(N = 25000,
                 alpha = alpha, 
                 B = B, 
                 delta = delta, 
                 eta = eta, 
                 w = w, 
                 v = 10000, #### social program with 832 as baseline
                 r = r, #### cartel salary 
                 kappaP = kappaP, 
                 kappaK = kappaK, 
                 pi_v = pi_v, 
                 pi_f = pi_f,
                 rho = rho)
}

#### impact of social programs. 250 simulations with different value of social program
{
Res <- data.frame(socialProgram = c(), 
                  Active = c(),
                  Killed = c(),
                  Incapacitated = c(),
                  NotRecruited = c(),
                  Exposed = c(),
                  Saturated = c())
for(k in 1:250){
  v_sim <- 20000*runif(1)
  IRes <- ImpactSP(N = 25000,
                  alpha = alpha, 
                  B = B, 
                  delta = delta, 
                  eta = eta, 
                  w = w, 
                  v = v_sim, #### social program with 832 as baseline
                  r = r, #### cartel salary 
                  kappaP = kappaP, 
                  kappaK = kappaK, 
                  pi_v = pi_v, 
                  pi_f = pi_f,
                  rho = rho)
  Res <- rbind(Res, data.frame(socialProgram = v_sim, 
                    Active = IRes$Active,
                    Killed = IRes$Killed,
                    Incapacitated = IRes$Incapacitated,
                    NotRecruited = IRes$NotRecruited,
                    Exposed = IRes$Exposed,
                    Saturated = IRes$Saturated))
}
save(Res, file = "output/ImpactSocialPrograms20250430.RData")
}

#### the sublinear income of a cartel, calculations
{
betaIncome <- (2/3 + 3/4)/2
AlphaIncome <- (r +omega)/ (betaIncome * N^(betaIncome-1))
CMax <- (AlphaIncome/(r+omega))^(1/(1-betaIncome))
CartelSize <- seq(from = 0, to = CMax, by = 100)
CIncome <- AlphaIncome* CartelSize^(betaIncome)
CCost <- CartelSize * (r+omega)
SocialProg <- seq(from = 0, to = 8000, by = 1)
Cstar <- ((r+omega + SocialProg-v)/(AlphaIncome*betaIncome))^(1/(betaIncome-1))
Revenue <-  (AlphaIncome* Cstar^betaIncome )- (r+omega + SocialProg-v)*Cstar

betaIncome1 <- 2/3
AlphaIncome1 <- (r+omega) / (betaIncome1 * N^(betaIncome1-1))
betaIncome2 <- 3/4
AlphaIncome2 <- (r+omega) / (betaIncome2 * N^(betaIncome2-1))
betaIncome3 <- (betaIncome1 + betaIncome2)/2
AlphaIncome3 <- (r+omega) / (betaIncome3 * N^(betaIncome3-1))

Cstar1 <- ((r+omega + SocialProg-v)/(AlphaIncome1*betaIncome1))^(1/(betaIncome1-1))
Cstar2 <- ((r+omega + SocialProg-v)/(AlphaIncome2*betaIncome2))^(1/(betaIncome2-1))
Cstar3 <- ((r+omega + SocialProg-v)/(AlphaIncome3*betaIncome3))^(1/(betaIncome3-1))

CartEff <- function(betaIncome, 
                    N,
                    r,
                    omega,
                    SocialProg,
                    v){
  AlphaIncome <- (r +omega)/ (betaIncome * N^(betaIncome-1))
  CIncome <- AlphaIncome* N^(betaIncome)
  CCost <- N * (r+omega)
  Revenue <-  (AlphaIncome* N^betaIncome )- (r+omega + SocialProg-v)*N
  return(list(c(CIncome, CCost, Revenue)/1000000))
}

CartEff(betaIncome = (2/3+3/4)/2,
        N = 25000,
        r,
        omega,
        SocialProg = 0*v,
        v)
}

#### FIGURES
#### colors
{
cols <- c("#00CC99", "#666699", "#E8686D", "#FFC007", "#FF6799", "#7F7F7F", "#6EC2E8")
}

#### 2 by 2 with individual model and collective model
{
  pdf("IndividualCollectiveImpacts.pdf", width = 8, height = 6)
  load(file = "output/ImpactSocialPrograms20250430.RData")
  par(mfrow = c(2, 2),        # 2 rows, 2 columns
      mar = c(4, 4, 1.5, 1),  # Increase bottom/top margins for more vertical space
      oma = c(0, 0, 0.5, 0),  # Reduce top outer margin to cut whitespace
      mgp = c(2, 0.5, 0))     # Label positioning
  #### piled figure for the effect of a program
  {
    RangeSocial <- seq(from = 0, to = 8000, by = 1)
    EffectSocial <- PlotIsRecruited(alpha = alpha, 
                                    B = B, 
                                    delta = delta, 
                                    eta = eta, 
                                    w = w, 
                                    v = RangeSocial,
                                    r = r, 
                                    kappaP = kappaP, 
                                    kappaK = kappaK, 
                                    pi_v = pi_v, 
                                    pi_f = pi_f)
    plot(RangeSocial, EffectSocial, type = "l",
         xlab = "Social program (U$)",
         ylab = "Probability of recruitment",
         col = NA, ylim = c(0, 1))
    polygon(c(RangeSocial, rev(RangeSocial)), 
            c(EffectSocial, 0*EffectSocial), col = cols[3])
    polygon(c(RangeSocial, rev(RangeSocial)), 
            c(EffectSocial, 0*EffectSocial+1), col = cols[7])
    points(RangeSocial, 0*RangeSocial + pi_f, type = "l", lty = 2)
    points(RangeSocial, 0*RangeSocial + pi_f+pi_v, type = "l", lty = 2)
    text(v+100, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 0)
    text(v+100, 0.205, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 0)
    text(v+100, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 0)
    points(c(v,v), c(0, 10000000), type = "l", lty = 2)
    
  }
  
  #### piled figure for the effect of the wage of cartels
  {
    RangeCartelWage <- seq(from = 0, to = 17000, by = 1)
    EffectCartelWage <- PlotIsRecruited(alpha = alpha, 
                                        B = B, 
                                        delta = delta, 
                                        eta = eta, 
                                        w = w, 
                                        v = v,
                                        r = RangeCartelWage, 
                                        kappaP = kappaP, 
                                        kappaK = kappaK, 
                                        pi_v = pi_v, 
                                        pi_f = pi_f)
    plot(RangeCartelWage, EffectCartelWage, type = "l",
         xlab = "Cartel wage (U$)",
         ylab =  "Probability of recruitment", 
         col = NA, ylim = c(0, 1))
    polygon(c(RangeCartelWage, rev(RangeCartelWage)), 
            c(EffectCartelWage, 0*EffectCartelWage), col = cols[3])
    polygon(c(RangeCartelWage, rev(RangeCartelWage)), 
            c(EffectCartelWage, 0*EffectCartelWage+1), col = cols[7])
    points(RangeCartelWage, 0*RangeCartelWage + pi_f, type = "l", lty = 2)
    points(RangeCartelWage, 0*RangeCartelWage + pi_f+pi_v, type = "l", lty = 2)
    text(15990, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 1)
    text(15990, 0.205, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 1)
    text(15990, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 1)
  }
  
  
  #### smooth curve - exposed 
  {
    SmExposed <- IV(Res$socialProgram, Res$Exposed, span = .1)
    f <- SmExposed$x < 8000
    plot(SmExposed$x, SmExposed$mean, col = NA,
         ylim = c(0, max(SmExposed$mean)),
         xlim = c(0, 8000),
         xlab = "Social program  (U$)",
         yaxt = "n",
         ylab = "Exposed individuals (thousands)")
    axis(2, at = (0:5)*25000, 
         labels = (0:5)*25)
    polygon(c(SmExposed$x[f], rev(SmExposed$x[f])),
            c(SmExposed$mean[f], 0* SmExposed$mean[f]),     
            col = cols[2])
    for(k in 1:10){points(SmExposed$x[f], 0*SmExposed$x[f] + k*25000, type = "l", col = "gray30", lty = 2)}
    points(c(v,v), c(0, 10000000), type = "l", lty = 2)
    text(900, 12500, "Exposed to C1 recruitment", col = "white", adj = 0,cex = 1)
  }
  
  #### smooth curve - exposed and recruited
  {
    SmExposed <- IV(Res$socialProgram, Res$Exposed - Res$NotRecruited, span = .1)
    plot(SmExposed$x[f], SmExposed$mean[f], col = NA,
         ylim = c(0, max(SmExposed$mean+SmNotR$mean)),
         xlim = c(0, 8000),
         yaxt = "n",
         xlab = "Social program  (U$)",
         ylab = "Exposed and recruited people")
    axis(2, at = (0:5)*25000, 
         labels = (0:5)*25)
    polygon(c(SmExposed$x[f], rev(SmExposed$x[f])),
            c(SmExposed$mean[f], 0* SmExposed$mean[f]),     
            col = cols[3])
    polygon(c(SmNotR$x[f], rev(SmNotR$x[f])),
            c(SmExposed$mean[f] + SmNotR$mean[f], 
              rev(SmExposed$mean[f])),     
            col = cols[1])
    for(k in 1:10){points(SmExposed$x[f], 0*SmExposed$x[f] + k*25000, type = "l", col = "gray30", lty = 2)}
    points(c(v,v), c(0, 10000000), type = "l", lty = 2)
    
    text(7900, 35000, "Not recruited", col = "white", adj = 1,cex = 1)
    text(900, 15000, "Recruited", col = "white", adj = 0,cex = 1)
    
  }
  
  dev.off()
}

#### 2 by 2 with individual model
{
  pdf("IndividualImpacts.pdf", width = 8, height = 6)
  par(mfrow = c(2, 2),        # 2 rows, 2 columns
      mar = c(4, 4, 1.5, 1),  # Increase bottom/top margins for more vertical space
      oma = c(0, 0, 0.5, 0),  # Reduce top outer margin to cut whitespace
      mgp = c(2, 0.5, 0))     # Label positioning
  #### piled figure for the effect of a program
  {
    RangeSocial <- seq(from = 0, to = 200*52, by = 1)
    EffectSocial <- PlotIsRecruited(alpha = alpha, 
                                    B = B, 
                                    delta = delta, 
                                    eta = eta, 
                                    w = w, 
                                    v = RangeSocial,
                                    r = r, 
                                    kappaP = kappaP, 
                                    kappaK = kappaK, 
                                    pi_v = pi_v, 
                                    pi_f = pi_f)
    plot(RangeSocial, EffectSocial, type = "l",
         xlab = "value of a social program",
         ylab = "Prob of recruitment",
         col = NA, ylim = c(0, 1))
    polygon(c(RangeSocial, rev(RangeSocial)), 
            c(EffectSocial, 0*EffectSocial), col = cols[3])
    polygon(c(RangeSocial, rev(RangeSocial)), 
            c(EffectSocial, 0*EffectSocial+1), col = cols[7])
    points(RangeSocial, 0*RangeSocial + pi_f, type = "l", lty = 2)
    points(RangeSocial, 0*RangeSocial + pi_f+pi_v, type = "l", lty = 2)
    text(10, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 0)
    text(10, 0.45, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 0)
    text(10, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 0)
  }
  
  #### piled figure for the effect of the wage of cartels
  {
    RangeCartelWage <- seq(from = 0, to = 20000, by = 1)
    EffectCartelWage <- PlotIsRecruited(alpha = alpha, 
                                        B = B, 
                                        delta = delta, 
                                        eta = eta, 
                                        w = w, 
                                        v = v,
                                        r = RangeCartelWage, 
                                        kappaP = kappaP, 
                                        kappaK = kappaK, 
                                        pi_v = pi_v, 
                                        pi_f = pi_f)
    plot(RangeCartelWage, EffectCartelWage, type = "l",
         xlab = "cartel wage",
         ylab = "Prob of recruitment",
         col = NA, ylim = c(0, 1))
    polygon(c(RangeCartelWage, rev(RangeCartelWage)), 
            c(EffectCartelWage, 0*EffectCartelWage), col = cols[3])
    polygon(c(RangeCartelWage, rev(RangeCartelWage)), 
            c(EffectCartelWage, 0*EffectCartelWage+1), col = cols[7])
    points(RangeCartelWage, 0*RangeCartelWage + pi_f, type = "l", lty = 2)
    points(RangeCartelWage, 0*RangeCartelWage + pi_f+pi_v, type = "l", lty = 2)
    text(19990, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 1)
    text(19990, 0.45, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 1)
    text(19990, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 1)
  }
  
  #### piled figure for the perceived cost of incapacitation
  {
    PerceivedInc <- seq(from = 0, to = 20000, by = 10)
    EffectPerceivedInc <- PlotIsRecruited(alpha = alpha, 
                                          B = B, 
                                          delta = PerceivedInc, 
                                          eta = eta, 
                                          w = w, 
                                          v = v,
                                          r = r, 
                                          kappaP = kappaP, 
                                          kappaK = kappaK, 
                                          pi_v = pi_v, 
                                          pi_f = pi_f)
    plot(PerceivedInc, EffectPerceivedInc, type = "l",
         xlab = "Perceived cost of incapacitation",
         ylab = "Prob. of recruitment",
         col = NA, ylim = c(0, 1))
    polygon(c(PerceivedInc, rev(PerceivedInc)), 
            c(EffectPerceivedInc, 0*EffectPerceivedInc), col = cols[3])
    polygon(c(PerceivedInc, rev(PerceivedInc)), 
            c(EffectPerceivedInc, 0*EffectPerceivedInc+1), col = cols[7])
    points(PerceivedInc, 0*PerceivedInc + pi_f, type = "l", lty = 2)
    points(PerceivedInc, 0*PerceivedInc + pi_f+pi_v, type = "l", lty = 2)
    text(90, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 0)
    text(90, 0.45, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 0)
    text(90, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 0)
  }
  
  #### piled figure for the perceived risk of being killed
  {
    PerceivedKappaK <- seq(from = 0, to = 5*kappaK, length.out = 1000)
    EffectPerceivedKappaK <- PlotIsRecruited(alpha = alpha, 
                                          B = B, 
                                          delta = delta, 
                                          eta = eta, 
                                          w = w, 
                                          v = v,
                                          r = r, 
                                          kappaP = kappaP, 
                                          kappaK = PerceivedKappaK, 
                                          pi_v = pi_v, 
                                          pi_f = pi_f)
    plot(PerceivedKappaK, EffectPerceivedKappaK, type = "l",
         xlab = "Perceived risks of being killed",
         ylab = "Prob. of recruitment",
         col = NA, ylim = c(0, 1))
    polygon(c(PerceivedKappaK, rev(PerceivedKappaK)), 
            c(EffectPerceivedKappaK, 0*EffectPerceivedKappaK), col = cols[3])
    polygon(c(PerceivedKappaK, rev(PerceivedKappaK)), 
            c(EffectPerceivedKappaK, 0*EffectPerceivedKappaK+1), col = cols[7])
    points(PerceivedKappaK, 0*PerceivedKappaK + pi_f, type = "l", lty = 2)
    points(PerceivedKappaK, 0*PerceivedKappaK + pi_f+pi_v, type = "l", lty = 2)
    text(.00010, 0.95, expression(paste("Not recruited - ", pi[i])), col = "white", adj = 0)
    text(.00010, 0.45, expression(paste("Voluntary - ", pi[v])), col = "white", adj = 0)
    text(.00010, 0.05, expression(paste("Forced recruitment - ", pi[f])), col = "white", adj = 0)
  }
  
  dev.off()
}
  
#### 2 by 2 for how cartel reacts PDF
{
pdf("CartelStrategy.pdf", width = 8, height = 6)
  par(mfrow = c(2, 2),        # 2 rows, 2 columns
      mar = c(4, 4, 1.5, 1),  # Increase bottom/top margins for more vertical space
      oma = c(0, 0, 0.5, 0),  # Reduce top outer margin to cut whitespace
      mgp = c(2, 0.5, 0))     # Label positioning

#### costs and income for cartel
{
  plot(CartelSize, CCost, type = "l", 
       xlab = "C1 size (thousands)",
       ylab = "Yearly U$ in millions",
       yaxt = "n", xaxt = "n",
       col = NA, lwd = 3)
  axis(2, at = 400*(0:15)*1e6, 
       labels = 400*(0:15)*1)
  axis(1, at = 10000*(0:40), 
       labels = 10*(0:40)*1)
  
  points(c(N,N), c(0, 8e9), type = "l", lty = 2)
  I <- CartelSize <= N
  polygon(c(CartelSize[I], rev(CartelSize[I])),
          c(CIncome[I], rev(CCost[I])), 
          col = cols[1])
  I <- CartelSize >= N
  polygon(c(CartelSize[I], rev(CartelSize[I])),
          c(CIncome[I], rev(CCost[I])), 
          col = cols[3])
  points(CartelSize, CCost, 
         col = cols[7],
         type = "l", lwd = 3)
  points(CartelSize, CIncome, type = "l", col = 1, lwd = 3)
  points(CartelSize, CCost+(AlphaIncome*N^betaIncome - (r+omega)*N), type = "l", lty = 2)
  arrows(x0 = 26000, y0 = 1280e6, 
         x1 = 38000, y1 = 1280e6,  # Change end point as needed
         col = cols[3], 
         length = 0.1) 
  
  arrows(x0 = 24000, y0 = 1280e6, 
         x1 = 14000, y1 = 1280e6,  # Change end point as needed
         col = cols[1], 
         length = 0.1) 
  text(26000, 1200e6, "less social programs", col =  cols[3], adj = 0)
  text(24000, 1200e6, "more", col =  cols[1], adj = 1)
  text(43000, 500e6, 
       col = cols[7],
       expression("Costs:   " ~ X(N) == (r + omega)*N), adj = 0)
  text(40550, 850e6, expression("Income:  " ~ Y(N) == alpha * N^beta), adj = 1)
  

  }

#### opt size for a cartel depending on social program 
{
  plot(SocialProg, Cstar, type = "l", 
       xlab = "Social program (U$)",
       ylab = "Optimal size (thousands)",
       yaxt = "n",
       ylim = c(0, max(Cstar)),
       col = NA, lwd = 3)
  axis(2, at = (0:5)*10000, 
       labels = (0:5)*10)
  I <- SocialProg <= v
  polygon(c(SocialProg[I], rev(SocialProg[I])),
          c(Cstar[I], 0*rev(Cstar[I])), 
          col = cols[3])
  I <- SocialProg >= v
  polygon(c(SocialProg[I], rev(SocialProg[I])),
          c(Cstar[I], 0*rev(Cstar[I])), 
          col = cols[1])
  points(c(v,v), c(0, max(Cstar)), type = "l", lty = 2)
  arrows(x0 = v*1.1, y0 = 30000, 
         x1 = v*3, y1 = 30000,  # Change end point as needed
         col = cols[1], 
         length = 0.1) 
  
  text(v*1.1, 28000, "more social programs", col =  cols[1], adj = 0)

}

#### cartel revenue
{
  plot(SocialProg, Revenue, type = "l", 
       xlab = "Social program U$ by year",
       ylab = "C1 profit (US$ million)",
       yaxt = "n",
       ylim = c(0, max(Revenue)),
       col = NA, lwd = 3)
  axis(2, at = 50*(0:5)*1e6, 
       labels = 50*(0:5)*1)
  I <- SocialProg <= v
  polygon(c(SocialProg[I], rev(SocialProg[I])),
          c(Revenue[I], 0*rev(Revenue[I])), 
          col = cols[3])
  I <- SocialProg >= v
  polygon(c(SocialProg[I], rev(SocialProg[I])),
          c(Revenue[I], 0*rev(Revenue[I])), 
          col = cols[1])
  points(c(v,v), c(0, max(Revenue)), type = "l", lty = 2)
  arrows(x0 = v*1.1, y0 = 1.8e8, 
         x1 = v*3, y1 = 1.8e8,  # Change end point as needed
         col = cols[1], 
         length = 0.1) 
  
  text(v*1.1, 1.7e8, "more social programs", col =  cols[1], adj = 0)
  
}

#### opt size for a cartel depending on social program for different beta
{
  plot(SocialProg, Cstar, type = "l", 
       xlab = "Social program U$ by year",
       ylab = "Optimal size (thousands)",
       yaxt = "n",
       ylim = c(0, 31000),
       col = NA, lwd = 3)
  axis(2, at = (0:5)*10000, 
       labels = (0:5)*10)
  points(SocialProg, Cstar1, type = "l", lwd = 2, 
         col = cols[4])
  points(SocialProg, Cstar2, type = "l", lwd = 2,
         col = cols[7])
  points(SocialProg, Cstar3, type = "l", lwd = 2,
         col = "black")
  points(c(v,v), c(0, 1.2*max(Cstar)), type = "l", lty = 2)
  arrows(x0 = v*1.1, y0 = 30000, 
         x1 = v*3, y1 = 30000,  # Change end point as needed
         col = cols[1], 
         length = 0.1) 
  
  text(v*1.1, 28000, "more social programs", col =  cols[1], adj = 0)
  
  text(6000, 14000, expression(beta == 2/3), col = cols[4])
  #text(6000, 25000, expression(beta == 0.83), col = cols[1])
  text(6000, 4500, expression(beta == 3/4), col = cols[7])
}
dev.off()
}

#### costs and income for cartel per person
{
  plot(CartelSize, CIncome/CartelSize, type = "l", 
       xlab = "C1 size",
       ylab = "Yearly U$ in millions",
       yaxt = "n",
       ylim = c(0, 20000),
       col = NA, lwd = 3)
  axis(2, at = 100*(0:15)*1e6, 
       labels = 100*(0:15)*1)
  
  points(c(N,N), c(0, 8e9), type = "l", lty = 2)
  I <- CartelSize <= N
  y <- AlphaIncome*betaIncome*(CartelSize^(betaIncome-1))
  polygon(c(CartelSize[I], rev(CartelSize[I])),
          c(y[I], 0*rev(CCost[I])), 
          col = cols[3])
  I <- CartelSize >= N
  polygon(c(CartelSize[I], rev(CartelSize[I])),
          c(y[I], 0*rev(CCost[I])), 
          col = cols[1])
  points(CartelSize, CCost/CartelSize, type = "l", lty = 2)
  points(CartelSize, CIncome/CartelSize, type = "l", col = 1, lwd = 3)
  points(CartelSize, CCost+(AlphaIncome*N^betaIncome - r*N), type = "l", lty = 2)
}

#### impact of Social programs on exposed, recruited and etc
{
load(file = "output/ImpactSocialPrograms20250430.RData")
plot(Res$socialProgram, Res$Exposed, type = "h", lwd = 3, col = "tomato", ylim = c(0, max(Res$Exposed)))
points(Res$socialProgram, Res$NotRecruited, type = "h", lwd = 3, 
       col = "blue", ylim = c(0, max(Res$Exposed)))

plot(Res$socialProgram, Res$NotRecruited, type = "h", lwd = 3, col = "tomato", ylim = c(0, max(Res$Exposed)))
plot(Res$socialProgram, Res$NotRecruited/Res$Exposed)
plot(Res$socialProgram, Res$Exposed - Res$NotRecruited, 
     ylab = "Recruited",
     type = "h", lwd = 3, col = "green3",
     ylim = c(0, max(Res$Exposed - Res$NotRecruited)))

plot(Res$socialProgram, Res$Killed, 
     type = "h", lwd = 3, col = "tomato")
}

#### 2 by 2 with smoothed results
{
  pdf("RecruitmentModel.pdf", width = 8, height = 6)
  load(file = "output/ImpactSocialPrograms20250430.RData")
  par(mfrow = c(2, 2),        # 2 rows, 2 columns
      mar = c(3, 3, 2, 1),    # Margins around individual plots: bottom, left, top, right
      oma = c(0, 0, 2, 0),    # Outer margins (can adjust if you use outer titles)
      mgp = c(2, 0.5, 0)) 

#### smooth curve - exposed 
{
  SmExposed <- IV(Res$socialProgram, Res$Exposed, span = .1)
  plot(SmExposed$x, SmExposed$mean, col = NA,
       ylim = c(0, max(SmExposed$mean)),
       xlim = c(0, 10000),
       xlab = "Social program",
       yaxt = "n",
       ylab = "Exposed (thousands)")
  axis(2, at = (0:5)*20000, 
       labels = (0:5)*20)
  polygon(c(SmExposed$x, rev(SmExposed$x)),
          c(SmExposed$mean, 0* SmExposed$mean),     
          col = cols[2])
  for(k in 1:10){points(SmExposed$x, 0*SmExposed$x + k*10000, type = "l", col = "gray30", lty = 2)}
  points(c(v,v), c(0, 10000000), type = "l", lwd = 2)
  text(900, 20000, "Exposed to C1 recruitment", col = "white", adj = 0,cex = 1.3)
  }

#### smooth curve - exposed and NotRecruited
{
  SmNotR <- IV(Res$socialProgram, Res$NotRecruited, span = .1)
  plot(SmNotR$x, SmNotR$mean, col = NA,
       ylim = c(0, max(SmNotR$mean)),
       xlim = c(0, 10000),
       yaxt = "n",
       xlab = "Social program",
       ylab = "Exposed and not recruited (thousands)")
  axis(2, at = (0:5)*20000, 
       labels = (0:5)*20)
  polygon(c(SmNotR$x, rev(SmNotR$x)),
          c(SmNotR$mean, 0* SmNotR$mean),     
          col = cols[3])
  for(k in 1:10){points(SmExposed$x, 0*SmExposed$x + k*10000, type = "l", col = "gray30", lty = 2)}
  points(c(v,v), c(0, 10000000), type = "l", lwd = 2)
  text(900, 10000, "Exposed and not recruited", col = "white", adj = 0,cex = 1.3)
  
}

#### smooth curve - exposed and recruited
{
  SmExposed <- IV(Res$socialProgram, Res$Exposed - Res$NotRecruited, span = .1)
  plot(SmExposed$x, SmExposed$mean, col = NA,
       ylim = c(0, max(SmExposed$mean+SmNotR$mean)),
       xlim = c(0, 10000),
       yaxt = "n",
       xlab = "Social program",
       ylab = "Exposed and recruited")
  axis(2, at = (0:5)*20000, 
       labels = (0:5)*20)
  polygon(c(SmExposed$x, rev(SmExposed$x)),
          c(SmExposed$mean, 0* SmExposed$mean),     
          col = cols[1])
  polygon(c(SmNotR$x, rev(SmNotR$x)),
          c(SmExposed$mean + SmNotR$mean, rev(SmExposed$mean)),     
          col = cols[3])
  for(k in 1:10){points(SmExposed$x, 0*SmExposed$x + k*10000, type = "l", col = "gray30", lty = 2)}
  points(c(v,v), c(0, 10000000), type = "l", lwd = 2)
  
  text(9900, 30000, "Not recruited", col = "white", adj = 1,cex = 1.3)
  text(900, 20000, "Recruited", col = "white", adj = 0,cex = 1.3)
  
}

#### smooth curve - killed
{
SmKilled <- IV(Res$socialProgram, Res$Killed, span = .1)
SmIncapacitated <- IV(Res$socialProgram, Res$Incapacitated, span = .1)
  
plot(SmKilled$x, SmKilled$mean, col = NA,
     ylim = c(0, max(SmKilled$mean + SmIncapacitated$mean)),
     xlim = c(0, 10000),
     yaxt = "n",
     xlab = "Social program",
     ylab = "C1 killed and incapacitated members")
axis(2, at = (0:5)*20000, 
     labels = (0:5)*20)
polygon(c(SmKilled$x, rev(SmKilled$x)),
        c(SmKilled$mean, 0* SmKilled$mean),     
        col = cols[4])
polygon(c(SmKilled$x, rev(SmKilled$x)),
        c(SmKilled$mean+SmIncapacitated$mean, rev(SmKilled$mean)),     
        col = cols[5])
for(k in 1:10){points(SmExposed$x, 0*SmExposed$x + k*10000, type = "l", col = "gray30", lty = 2)}
points(c(v,v), c(0, 10000000), type = "l", lwd = 2)
text(9900, 2500, "Killed", col = "white", adj = 1,cex = 1.3)
text(9900, 7000, "Incapacitated", col = "white", adj = 1,cex = 1.3)

}
dev.off()
}

#### piled figure for all status
{
  t <- 1:length(Efs$At)
  y1 <- Efs$At
  y2 <- cumsum(Efs$Kt)
  y3 <- cumsum(Efs$It)
  y4 <- cumsum(Efs$Nt)
  y5 <- 0*Efs$It
  plot(y1, type = "l", 
       xlab = "time (weeks)",
       ylab = "People",
       col = NA, ylim = c(0, 1.1*max(y1+y2+y3+y4+y5)))
  
  polygon(c(t, rev(t)), c(y1, 0*y1), col = cols[3])
  polygon(c(t, rev(t)), c(y1 + y2, rev(y1)+0*y2), col = cols[4])
  polygon(c(t, rev(t)), c(y3 + y2 + y1, rev(y1 + y2)+0*y3), col = cols[2])
  polygon(c(t, rev(t)), c(y4 + y3 + y2 + y1, rev(y1+y2+y3)+0*y4), col = cols[1])
  polygon(c(t, rev(t)), c(y5 + y4 + y3 + y2 + y1, rev(y1+y2+y3+y4)+0*y5), col = cols[5])
  text(length(y1)-10, y1[length(y1)]/2, "Cartel members", col = "white", adj = 1)
  text(length(y1)-10, y1[length(y1)] + y2[length(y1)]/2, "Killed", col = "white", adj = 1)
  text(length(y1)-10, y1[length(y1)] + y2[length(y1)] + y3[length(y1)]/2, "Incapacitated", col = "white", adj = 1)
  text(length(y1)-10, y1[length(y1)] + y2[length(y1)] + y3[length(y1)] + y4[length(y1)]/2, "Exposed and not recruited", col = "black", adj = 1)
}
