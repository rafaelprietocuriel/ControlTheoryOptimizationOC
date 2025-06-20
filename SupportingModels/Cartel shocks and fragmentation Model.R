# ---------------------------------------------------------------
# Author: Rafael Prieto-Curiel
# Date: May 2025
# Title: Conflict Shocks and Fragmentation – Time Series Simulation
#
# Description:
# This R script models shocks and fragmentation by varying the conflict 
# parameter over time as a time series. It enables the quantification 
# of shock magnitude and its impact on system fragmentation. 
# The script produces visual outputs to support the analysis.
# ---------------------------------------------------------------

#### parameters
{
require(scales)
N1 <- 25000 #### number of initial people
N2<- 25000 
sat <- 50/(N1^2)
rho = 150/25000 #### rate of recruitment each week with base = 1/2
kappaP = 800/(52*25000) #### probability of arrest
kappaK = 700/(52*25000) #### probability of killed 
}

#### functions
{
OneWeek <- function(C1, C2, 
                    rho = 150 / 25000, 
                    theta = 2*17 / (25000*25000),
                    eta = 80 / (2*25000), 
                    omega = 76 / (25000^2)) {
  Kills <- 2*theta*C1*C2
  Arrests <- eta*(C1+C2)
  C1N <- C1 + rho*C1 - theta*C1*C2 - eta*C1 - omega*C1*C1
  C2N <- C2 + rho*C2 - theta*C1*C2 - eta*C2 - omega*C2*C2
  C1 <- C1N
  C2 <- C2N
  return(data.frame(C1 = C1, 
                    C2 = C2, 
                    theta = theta,
                    Kills = Kills, 
                    Arrests = Arrests))
  
}

  OneWeekCollapse <- function(C1, C2, 
                      rho = 150 / 25000, 
                      theta = 2*17 / (25000*25000),
                      eta = 80 / (2*25000), 
                      omega = 76 / (25000^2)) {
    Kills <- 2*theta*C1*C2
    Arrests <- eta*(C1+C2)
    C1N <- C1 + rho*C1 - theta*C1*C2 - eta/2 - omega*C1*C1
    C2N <- C2 + rho*C2 - theta*C1*C2 - eta/2 - omega*C2*C2
    C1 <- max(C1N, 0)
    C2 <- max(C2N, 0)
    return(data.frame(C1 = C1, 
                      C2 = C2, 
                      theta = theta,
                      Kills = Kills, 
                      Arrests = Arrests))
    
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

#### 1500 simulations for steady state
{
ResW <- OneWeek(C1 = N1, C2 = N2)  
Res <- data.frame(C1 = ResW$C1, 
                  C2 = ResW$C2, 
                  theta = ResW$theta,
                  Kills = ResW$Kills, 
                  Arrests = ResW$Arrests)
for(k in 1:1500){ ### 30 years
  ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2)
  Res <- rbind(Res, data.frame(C1 = ResW$C1,
                               C2 = ResW$C2, 
                               theta = ResW$theta,
                               Kills = ResW$Kills, 
                               Arrests = ResW$Arrests))
}
ResBaseline <- Res
save(ResBaseline, file = "output/Retaliation20250517.RData")
}

#### 1500 simulations for steady state with 6 theta
{
  thetaM <- rep(2*17 / (25000*25000), 1500)
  thetaM[100:125] <- thetaM[100:125]*6 ### so 26 weeks or 6 months
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = thetaM[1],
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  for(k in 1:1500){ ### 30 years
    ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2,
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResSixTheta <- Res
  save(ResSixTheta, file = "output/RetaliationSixTheta20250517.RData")
}

#### 1500 simulations for steady state with 3 theta
{
  thetaM <- rep(2*17 / (25000*25000), 1500)
  thetaM[100:125] <- thetaM[100:125]*3 ### so 26 weeks or 6 months
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = thetaM[1],
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  for(k in 1:1500){ ### 30 years
    ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2,
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResThreeTheta <- Res
  save(ResThreeTheta, file = "output/RetaliationThreeTheta20250517.RData")
}

#### 1500 simulations for steady state with 2 theta
{
  thetaM <- rep(2*17 / (25000*25000), 1500)
  thetaM[100:125] <- thetaM[100:125]*2 ### so 26 weeks or 6 months
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = thetaM[1],
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  for(k in 1:1500){ ### 30 years
    ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2,
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResTwoTheta <- Res
  save(ResTwoTheta, file = "output/RetaliationTwoTheta20250517.RData")
}

#### 1500 simulations for steady state with 0 theta
{
  thetaM <- rep(2*17 / (25000*25000), 1500)
  thetaM[100:125] <- thetaM[100:125]*0 ### so 26 weeks or 6 months
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = thetaM[1],
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  for(k in 1:1500){ ### 30 years
    ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2,
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResZeroTheta <- Res
  save(ResZeroTheta, file = "output/RetaliationZeroTheta20250517.RData")
}

#### vary theta simulations
{
VaryT <- data.frame(ThetaEff = c(),
                    NewKills = c(),
                    oldKills = c(),
                    changeKills = c())

for(j in 1:250){
  ThetaEff <- runif(1)*11
  thetaM <- rep(2*17 / (25000*25000), 1500)
  thetaM[100:125] <- thetaM[100:125]*ThetaEff
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = thetaM[1],
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  for(k in 1:1500){ ### 30 years
    ResW <- OneWeek(C1 = ResW$C1, C2 = ResW$C2,
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResThreeTheta <- Res
  
  VaryT <- rbind(VaryT, data.frame(ThetaEff = ThetaEff,
                      NewKills = sum(ResThreeTheta$Kills),
                      oldKills = sum(ResBaseline$Kills),
                      changeKills = (sum(ResThreeTheta$Kills)-sum(ResBaseline$Kills))/sum(ResBaseline$Kills)))
  
}
save(VaryT, file = "output/VaryTheta20250517.RData")
}

#### 1500 simulations for collapse state
{
  ResW <- OneWeek(C1 = N1, C2 = N2)  
  Res <- data.frame(C1 = ResW$C1, 
                    C2 = ResW$C2, 
                    theta = ResW$theta,
                    Kills = ResW$Kills, 
                    Arrests = ResW$Arrests)
  thetaM <- rep(10 / (25000*25000), 1500)
  #thetaM[100:125] <- thetaM[100:125]*38
  thetaM[100:125] <- thetaM[100:125]*68
  for(k in 1:1500){ ### 30 years

    ResW <- OneWeekCollapse(C1 = ResW$C1, C2 = ResW$C2,
                    rho = 120/25000,
                    eta = 100 ,
                    omega = 60 / (25000^2),
                    theta = thetaM[k])
    Res <- rbind(Res, data.frame(C1 = ResW$C1,
                                 C2 = ResW$C2, 
                                 theta = ResW$theta,
                                 Kills = ResW$Kills, 
                                 Arrests = ResW$Arrests))
  }
  ResCollapse <- Res
  plot(ResCollapse$C1, type = "l", ylim = c(0, 25000))
  points(ResCollapse$C2, type = "l", ylim = c(0, 25000))
  save(ResCollapse, file = "output/RetaliationCollapse20250517.RData")
}

#### FIGURES
#### colors
{
cols <- c("#00CC99", "#666699", "#E8686D", "#FFC007", "#FF6799", "#7F7F7F", "#6EC2E8")
}

#### 2 by 2 with individual model and collective model
{
  pdf("RetaliationDynamics.pdf", width = 8, height = 6)
  load(file = "output/Retaliation20250517.RData")
  load(file = "output/VaryTheta20250517.RData")
  load(file = "output/RetaliationThreeTheta20250517.RData")
  load(file = "output/RetaliationSixTheta20250517.RData")
  load(file = "output/RetaliationTwoTheta20250517.RData")
  load(file = "output/RetaliationZeroTheta20250517.RData")
  par(mfrow = c(2, 2),        # 2 rows, 2 columns
      mar = c(4, 4, 1.5, 1),  # Increase bottom/top margins for more vertical space
      oma = c(0, 0, 0.5, 0),  # Reduce top outer margin to cut whitespace
      mgp = c(2, 0.5, 0))     # Label positioning
  #### piled figure kills
  {
    t <- 1:length(ResSixTheta$C1)
    f <- t<=520
    plot(t, ResThreeTheta$Kills, 
         xlab = "time (weeks)",
         ylab = "C1 - C2 homicides",
         col = NA, 
         type = "l",
         xlim = c(0, 520),
         ylim = c(0, 210))
    polygon(c(t[f], rev(t[f])), 
            c(ResThreeTheta$Kills[f], 0*ResSixTheta$Kills[f]), col = cols[3])
    points(t[f], ResBaseline$Kills[f], type = "l", lty = 2)
    text(90, 150, "shock", col = "black", adj = 1)
    text(90, 130,  expression(kappa == 2), col = "black", adj = 1)
    text(130, 150, "six months later", col = "black", adj = 0)
    text(230, ResBaseline$Kills[1]+20, "baseline", col = "black", adj = 0)
  }
  
  #### cumulative number of kills (%)
  {
    t <- 1:length(ResThreeTheta$C1)
    yBASE <- cumsum(ResBaseline$Kills)/sum(ResBaseline$Kills[1:500])
    yINCREASE3 <- cumsum(ResThreeTheta$Kills)/sum(ResBaseline$Kills[1:500])
    yINCREASE6 <- cumsum(ResSixTheta$Kills)/sum(ResBaseline$Kills[1:500])
    yINCREASE2 <- cumsum(ResTwoTheta$Kills)/sum(ResBaseline$Kills[1:500])
    yINCREASE0 <- cumsum(ResZeroTheta$Kills)/sum(ResBaseline$Kills[1:500])
    f <- t<520
    plot(t, yINCREASE, 
         xlab = "time (weeks)",
         ylab = "change in homicides (%)",
         col = NA, 
         type = "l",
         yaxt = "n",
         xlim = c(0, 520),
         ylim = c(0, 1.2))
axis(2, at = (0:10)*.2, 
     labels = (0:10)*20)
    
    points(t[f], yBASE[f], type = "l", lty = 2)
    points(t[f], yINCREASE0[f], type = "l", lwd = 3, col = cols[1])
    points(t[f], yINCREASE2[f], type = "l", lwd = 3, col = cols[2])
    points(t[f], yINCREASE3[f], type = "l", lwd = 3, col = cols[3])
    points(c(100,100),c(0,25000), type = "l", lty = 2)
    points(c(126,126),c(0,25000), type = "l", lty = 2)
    
    
    text(90, 0.8, "shock", col = "black", adj = 1)
    text(130, 0.8, "six months later", col = "black", adj = 0)
    text(400, 1,  expression(kappa == 2), col = cols[3], adj = 1)
    text(400, 0.91,  expression(kappa == 1), col = cols[2], adj = 1)
    text(400, 0.6,  expression(kappa == 0), col = cols[1], adj = 1)
 
  }
  
  
  #### C1 size
  {
    t <- 1:length(ResThreeTheta$C1)
    f <- t<= 520
    plot(t, ResThreeTheta$Kills, 
         xlab = "time (weeks)",
         ylab = "C1 size (thousands)",
         col = NA, 
         yaxt = "n",
         type = "l",
         xlim = c(0, 520),
         ylim = c(0, 28000))
    axis(2, at = (0:10)*5000, 
         labels = (0:10)*5)
    polygon(c(t[f], rev(t[f])), 
            c(ResThreeTheta$C1[f], 0*ResThreeTheta$Kills[f]), col = cols[1])
    points(t[f], ResBaseline$C1[f], type = "l", lty = 2)
    points(c(100,100),c(0,25000), type = "l", lty = 2)
    points(c(126,126),c(0,25000), type = "l", lty = 2)
    text(90, 8000, "shock", col = "black", adj = 1)
    text(90, 6000,  expression(kappa == 2), col = "black", adj = 1)
    text(130, 8000, "six months later", col = "black", adj = 0)
    text(230, ResBaseline$C1[1]+2000, "baseline", col = "black", adj = 0)
  }
  
  #### smooth kills change 
  {
    SmEffect <- IV(VaryT$ThetaEff, VaryT$NewKills/sum(ResBaseline$Kills), span = .1)
    plot(SmEffect$x-1, SmEffect$mean, col = NA,
         #ylim = c(0, max(SmEffect$mean)),
         xlim = c(-1, 10),
         xlab = expression("Shock size " * kappa),         
         yaxt = "n",
         ylim = c(0.85,1.05),
         ylab = "Change in C1 - C2 homicides (%)")
    axis(2, at = 0.8+(0:10)*.05, 
         labels = (0.8+(0:10)*0.05)*100)
    f <- SmEffect$x>1
    polygon(c(SmEffect$x[f], rev(SmEffect$x[f]))-1,
            c(SmEffect$mean[f], 0* SmEffect$mean[f]),     
            col = cols[6])
    f <- SmEffect$x<1
    polygon(c(SmEffect$x[f], rev(SmEffect$x[f]))-1,
            c(SmEffect$mean[f], 0* SmEffect$mean[f]),     
            col = cols[7])
    points(c(-1,10), c(1,1), type = "l", col = "gray30", lty = 2)
    text(9.3, 0.99, "baseline - 100%", col = "white", adj = 1,cex = 1)
  }
  dev.off()
}

