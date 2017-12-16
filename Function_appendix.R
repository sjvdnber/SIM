### Example data ###
setwd("~/Documents/PHD/Paper Mediation 1/FINAL")

library(foreign)

full <- read.spss("DataExampleTom.sav",to.data.frame = TRUE)
data <- as.data.frame(subset(full,filter_.=="Selected"))

data$educ <- ifelse(data$education=="secondary school (six years)" |
                    data$education=="secondary school (first four years)",0,1)
data$prof <- ifelse(data$profession=="houseman" | data$profession=="unemployed" |
                    data$profession=="retired",0,ifelse(data$profession=="student",1,2))
data$prof1 <- ifelse(data$prof==1,1,0)
data$prof2 <- ifelse(data$prof==2,1,0)
data$sex <- as.numeric(data$sex)-1
data$rel <- ifelse(data$MaritalSatus=="single",0,1)
data$cage <- data$age - mean(data$age)
data$cage2 <- data$cage*data$cage


### Specifying variables and models ###

Y <- data$PAIN # outcome
M <- data$EVALUATION # mediator
E <- as.numeric(data$condition)-1 # exposure

# Predictors outcome model in the treated (first model)
Xoutc1 <- cbind(data$cage,data$cage2)
# Xoutc1 <- cbind(data$sex,data$cage,data$cage2,data$sex*data$cage,data$sex*data$cage2)

# Predictors outcome model in the untreated (first model)
Xoutc0 <- cbind(data$prof1,data$prof2,data$rel)
# Xoutc0 <- cbind(data$educ,data$prof1,data$prof2,data$cage,data$rel,data$cage2,
#                 data$educ*data$cage,data$educ*data$cage2,data$prof1*data$cage,data$prof1*data$rel,
#                 data$prof1*data$cage2,data$prof2*data$cage,data$prof2*data$rel,data$prof2*data$cage2)

# Predictors third model in the untreated
Xthird0 <- cbind(data$prof1,data$prof2,data$cage)
# Xthird0 <- cbind(data$educ,data$sex,data$cage,data$cage2,data$educ*data$sex,data$sex*data$cage,
#            data$sex*data$cage2) 

# Predictors second model in the treated
Xsec1 <- cbind(data$educ,data$sex,data$rel)
# Xsec1 <- cbind(data$educ,data$sex,data$cage,data$cage2,data$educ*data$cage,data$educ*data$cage2,
#          data$sex*data$cage,data$sex*data$cage2,data$cage*data$cage2)

# Predictors second model in the untreated
Xsec0 <- cbind(data$prof1,data$prof2,data$rel)
# Xsec0 <- cbind(data$educ,data$prof1,data$prof2,data$sex,data$cage,data$rel,data$cage2,data$educ*data$sex,
#          data$educ*data$cage,data$educ*data$rel,data$educ*data$cage2,data$prof1*data$rel,
#          data$prof1*data$cage2,data$prof2*data$rel,data$prof2*data$cage2)

## mediation formula:
library('mediation')
m.model <- lm(EVALUATION~E+sex+educ+rel+cage,data=data)
y.model <- lm(PAIN~E+EVALUATION+prof1+prof2+rel,data=data)
summary(mediate(m.model,y.model,sims=1000,boot=T,treat="E",mediator="EVALUATION"))

### Function ###


eff.mediate <- function(Y,M,E,Xoutc1,Xoutc0,Xsec1,Xsec0,Xthird0,typ)
                    {
                      if (typ=="continuous") {
                      
                      d <- data.frame(Y=Y,M=M,E=E,Xoutc1=I(Xoutc1),Xoutc0=I(Xoutc0),Xsec1=I(Xsec1),
                           Xsec0=I(Xsec0),Xthird0=I(Xthird0))
                        
                      ## Restricted MLE:
                      ymfit1 <- lm(Y~Xoutc1+M,data=d[d$E==1,])
                      ymfit0 <- lm(Y~Xoutc0+M,data=d[d$E==0,])
                      d$pred1 <- predict(ymfit1,newdata=d,type="response")
                      d$pred0 <- predict(ymfit0,newdata=d,type="response")
                          
                      y1m0.rmle <- mean(predict(ymfit1,newdata=d[d$E==0,],type='response'))
                      y0m0.rmle <- mean(predict(ymfit0,newdata=d[d$E==0,],type='response'))
                      y1m1.rmle <- mean(predict(ymfit1,newdata=d[d$E==1,],type='response'))
                      direct.rmle <- mean(predict(ymfit1,newdata=d[d$E==0,],type='response')) - 
                                     mean(predict(ymfit0,newdata=d[d$E==0,],type='response'))
                      indirect.rmle <- mean(predict(ymfit1,newdata=d[d$E==1,],type='response')) - 
                                       mean(predict(ymfit1,newdata=d[d$E==0,],type='response'))
                      
                      ## Efficient estimators:
                      d$P <- mean(d$E)
                      
                      
                      ## Regression Model for Covariates on Outcome in the Control Arm:
                      
                      ## Approach 1:
                      yfit0 <- predict(lm(Y~Xsec0,data=d[d$E==0,],weights=(1/(1-d$P[d$E==0]))),
                               newdata=d,type="response")
                      
                      ## Approach 2:
                      d$outc_1 <- (Y*(1-d$E)/(1-d$P)) + ((d$E-d$P)/(1-d$P))*y0m0.rmle
                      d$x_vect <- t(t(Xsec0)-apply(Xsec0,2,mean))*((d$P-d$E)/(1-d$P))
                      lm0a <- as.vector(coef(lm(outc_1~x_vect,data=d)))
                      yfit0a <- apply((t(Xsec0)-apply(Xsec0,2,mean))*lm0a[-1],2,sum)
                      
                      
                      ## Regression Model for Covariates on Outcome in the Experimental Arm:
                      
                      ## Approach 1:
                      yfit1 <- predict(lm(Y~Xsec1,data=d[d$E==1,],weights=(1/(d$P[d$E==1]))),
                               newdata=d,type="response")
                      
                      ## Approach 2:
                      d$outc_2 <- Y*(d$E)/(d$P) - ((d$E-d$P)/(d$P))*y1m1.rmle
                      d$x_vect2 <- t(t(Xsec1)-apply(Xsec1,2,mean))*((d$E-d$P)/(d$P))
                      lm1a <- as.vector(coef(lm(outc_2~x_vect2,data=d)))
                      yfit1a <- apply((t(Xsec1)-apply(Xsec1,2,mean))*lm1a[-1],2,sum)
                
                      
                      ## Regress Fitted Values from Outcome Model on Covariates in Control Arm:
                      
                      ## Approach 1:
                      z1 <- predict(lm(pred1~Xthird0,data=d[d$E==0,],weights=(1/(1-d$P[d$E==0]))),
                            newdata=d,type="response")
                    
                      ## Approach 2:
                      d$o1 <- d$pred1*(1-d$E)/(1-d$P) + ((d$E-d$P)/(1-d$P))*y1m0.rmle
                      d$xv <- t(t(Xthird0)-apply(Xthird0,2,mean))*((d$P-d$E)/(1-d$P))
                      B <- as.vector(coef(lm(o1~xv,data=d)))
                      z2 <- apply((t(Xthird0)-apply(Xthird0,2,mean))*B[-1],2,sum)
                      
                      
                      ## Proposal I: locally efficient estimators
                      y0m0.le <- mean(yfit0)
                      y1m1.le <- mean(yfit1)
                      y1m0.le <- mean(z1)
                      direct.le <- mean(z1) - mean(yfit0)
                      indirect.le <- mean(yfit1) - mean(z1)
                      
                      
                      ## Proposal II: restricted efficient estimators
                      y0m0.re <- mean((1/(1-d$P))*(((1-d$E)*Y)-(yfit0a*(d$P-d$E))))
                      y1m1.re <- mean((1/(d$P))*(((d$E)*Y)-(yfit1a*(d$E-d$P))))
                      y1m0.re <- mean((1/(1-d$P))*(((1-d$E)*d$pred1)-(z2*(d$P-d$E))))
                      direct.re <- y1m0.re - y0m0.re
                      indirect.re <- y1m1.re - y1m0.re
                      
                      
                      ## Standard errors restricted MLE:
                      col_Xym1 <- cbind(1,d$M,Xoutc1)
                      col_Xym0 <- cbind(1,d$M,Xoutc0)
                      
                      row_Xym1 <- rbind(1,d$M,t(Xoutc1))
                      row_Xym0 <- rbind(1,d$M,t(Xoutc0))
                      
                      # Constructing U:
                      r_y1m0 <- d$Y-d$pred1
                      U_eta_y1m0 <- (d$E*col_Xym1)*r_y1m0
                      U_gamma_y1m0 <- (d$pred1 - y1m0.rmle)*(1-d$E)
                      U_y1m0 <- cbind(U_eta_y1m0,U_gamma_y1m0)
                      
                      r_y0m0 <- d$Y-d$pred0
                      U_eta_y0m0 <- ((1-d$E)*col_Xym0)*r_y0m0
                      U_gamma_y0m0 <- (d$pred0 - y0m0.rmle)*(1-d$E)
                      U_y0m0 <- cbind(U_eta_y0m0,U_gamma_y0m0)
                      
                      r_y1m1 <- d$Y-d$pred1
                      U_eta_y1m1 <- (d$E*col_Xym1)*r_y1m1
                      U_gamma_y1m1 <- (d$pred1 - y1m1.rmle)*(d$E)
                      U_y1m1 <- cbind(U_eta_y1m1,U_gamma_y1m1)
                      
                      # Constructing dU:
                      last_row1 <- cbind(col_Xym1,rep(-1,length(d$E)))
                      last_row0 <- cbind(col_Xym0,rep(-1,length(d$E)))
                      
                      EU_y1m0 <- rbind((cbind(((row_Xym1%*%(-d$E*col_Xym1))/length(d$E)),
                                 rep(0,nrow(row_Xym1)))),apply((1-d$E)*last_row1,2,mean))
                      EU_y0m0 <- rbind((cbind(((row_Xym0%*%(-(1-d$E)*col_Xym0))/length(d$E)),
                                 rep(0,nrow(row_Xym0)))),apply((1-d$E)*last_row0,2,mean))
                      EU_y1m1 <- rbind((cbind(((row_Xym1%*%(-d$E*col_Xym1))/length(d$E)),
                                 rep(0,nrow(row_Xym1)))),apply(d$E*last_row1,2,mean))
                      
                      U_direct <- cbind(U_y1m0,U_y0m0)
                      U_indirect <- cbind(U_y1m1,U_y1m0)
                      
                      mat0_y0m0 <- matrix(rep(0,nrow(EU_y0m0)*nrow(EU_y1m0)),ncol=nrow(EU_y1m0))
                      mat0_y1m0_1 <- matrix(rep(0,nrow(EU_y1m0)*nrow(EU_y0m0)),ncol=nrow(EU_y0m0))
                      mat0_y1m0_2 <- matrix(rep(0,nrow(EU_y1m0)*nrow(EU_y1m1)),ncol=nrow(EU_y1m1))
                      mat0_y1m1 <- matrix(rep(0,nrow(EU_y1m1)*nrow(EU_y1m0)),ncol=nrow(EU_y1m1))
                      EU_direct_p1 <- rbind(EU_y1m0,mat0_y0m0)
                      EU_direct_p2 <- rbind(mat0_y1m0_1,EU_y0m0)
                      EU_direct <- cbind(EU_direct_p1,EU_direct_p2)
                      EU_indirect_p1 <- rbind(EU_y1m1,mat0_y1m0_2)
                      EU_indirect_p2 <- rbind(mat0_y1m1,EU_y1m0)
                      EU_indirect <- cbind(EU_indirect_p1,EU_indirect_p2)
                      
                      # Variance U:
                      VAR_U_y1m0 <- var(U_y1m0)
                      VAR_U_y0m0 <- var(U_y0m0)
                      VAR_U_y1m1 <- var(U_y1m1)
                      
                      VAR_U_direct <- var(U_direct)
                      VAR_U_indirect <- var(U_indirect)
                      
                      # SE's:
                      SE.y1m0.rmle <- sqrt((1/length(d$E))*(solve(EU_y1m0)%*%VAR_U_y1m0%*%
                                      t(solve(EU_y1m0)))[nrow(EU_y1m0),nrow(EU_y1m0)])
                      SE.y0m0.rmle <- sqrt((1/length(d$E))*(solve(EU_y0m0)%*%VAR_U_y0m0%*%
                                      t(solve(EU_y0m0)))[nrow(EU_y0m0),nrow(EU_y0m0)])
                      SE.y1m1.rmle <- sqrt((1/length(d$E))*(solve(EU_y1m1)%*%VAR_U_y1m1%*%
                                      t(solve(EU_y1m1)))[nrow(EU_y1m1),nrow(EU_y1m1)])
                      
                      SE.direct.rmle <- sqrt((1/length(d$E))*((solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0),nrow(EU_y1m0)]) +
                                        ((1/length(d$E))*(solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0)+nrow(EU_y0m0),nrow(EU_y1m0)+nrow(EU_y0m0)])-
                                        2*((1/length(d$E))*(solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0),nrow(EU_y1m0)+nrow(EU_y0m0)]))
                      SE.indirect.rmle <- sqrt((1/length(d$E))*((solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1),nrow(EU_y1m1)]) +
                                          ((1/length(d$E))*(solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1)+nrow(EU_y1m0),nrow(EU_y1m1)+nrow(EU_y1m0)])-
                                          2*((1/length(d$E))*(solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1),nrow(EU_y1m1)+nrow(EU_y1m0)]))
                      
                      
                      ## Standard errors LE:
                      U_alpha <- d$E - d$P
                      U_gamma_eff_y1m0 <- ((d$pred1 - y1m0.le)*(1-d$E)) + (z1*(d$E-d$P))
                      U_eff_y1m0 <- cbind(U_alpha,U_eta_y1m0,U_gamma_eff_y1m0)
                      VAR_U_eff_y1m0 <- var(U_eff_y1m0)
                      
                      last_eff <- cbind(((1-d$E)%*%col_Xym1)/length(d$E),mean(-(1-d$E)))
                      last_3col_eff <- cbind(mean(-z1),last_eff)
                      
                      EU_eff_y1m0 <- rbind(c(-1,rep(0,ncol(col_Xym1)+1)),cbind(rep(0,ncol(col_Xym1)),
                                     ((row_Xym1%*%(-d$E*col_Xym1))/length(d$E)),rep(0,ncol(col_Xym1))),
                                     last_3col_eff)
                      SE.y1m0.le <- sqrt((1/length(d$E))*(solve(EU_eff_y1m0)%*%VAR_U_eff_y1m0%*%
                                    t(solve(EU_eff_y1m0)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)])
                      
                      U_gamma_eff_y0m0 <- (((1-d$E)/(1-d$P))*(d$Y-yfit0))+yfit0-y0m0.le
                      U_eff_y0m0 <- cbind(U_alpha,U_gamma_eff_y0m0)
                      VAR_U_eff_y0m0 <- var(U_eff_y0m0)
                      
                      last_eff_y0m0 <- c(mean(((1-d$E)/(1-d$P)^2)*(d$Y-yfit0)),-1)
                      EU_eff_y0m0 <- rbind(c(-1,0),last_eff_y0m0)
                      SE.y0m0.le <- sqrt((1/length(d$E))*(solve(EU_eff_y0m0)%*%VAR_U_eff_y0m0%*%
                                    t(solve(EU_eff_y0m0)))[nrow(EU_eff_y0m0),ncol(EU_eff_y0m0)])
                      
                      U_gamma_eff_y1m1 <- (((d$E)/(d$P))*(d$Y-yfit1))+yfit1-y1m1.le
                      U_eff_y1m1 <- cbind(U_alpha,U_gamma_eff_y1m1)
                      VAR_U_eff_y1m1 <- var(U_eff_y1m1)
                      
                      last_eff_y1m1 <- c(-mean(((d$E)/(d$P)^2)*(d$Y-yfit1)),-1)
                      EU_eff_y1m1 <- rbind(c(-1,0),last_eff_y1m1)
                      SE.y1m1.le <- sqrt((1/length(d$E))*(solve(EU_eff_y1m1)%*%VAR_U_eff_y1m1%*%
                                    t(solve(EU_eff_y1m1)))[nrow(EU_eff_y1m1),ncol(EU_eff_y1m1)])
                      
                      U_eff_direct <- cbind(U_eff_y1m0,U_gamma_eff_y0m0)
                      VAR_U_eff_direct <- var(U_eff_direct)
                      EU_eff_direct <- rbind(cbind(EU_eff_y1m0,(rep(0,nrow(EU_eff_y1m0)))),
                                       c(EU_eff_y0m0[-1,1],rep(0,(nrow(EU_eff_y1m0)-1)),
                                       EU_eff_y0m0[-1,-1]))
                      SE.direct.le <- sqrt(((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)]) +
                                      ((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[ncol(U_eff_direct),ncol(U_eff_direct)])-
                                      2*((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[nrow(EU_eff_y1m0),ncol(U_eff_direct)]))
                      
                      U_eff_indirect <- cbind(U_eff_y1m0,U_gamma_eff_y1m1)
                      VAR_U_eff_indirect <- var(U_eff_indirect)
                      EU_eff_indirect <- rbind(cbind(EU_eff_y1m0,(rep(0,nrow(EU_eff_y1m0)))),
                                         c(EU_eff_y1m1[-1,1],rep(0,(nrow(EU_eff_y1m0)-1)),
                                         EU_eff_y1m1[-1,-1]))
                      SE.indirect.le <- sqrt(((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)]) +
                                        ((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[ncol(U_eff_indirect),ncol(U_eff_indirect)])-
                                        2*((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[nrow(EU_eff_y1m0),ncol(U_eff_indirect)]))
                      
                      
                      ## Standard errors RE:
                      U_gamma_eff_y1m0_2 <- ((d$pred1 - y1m0.re)*(1-d$E)) + (z2*(d$E-d$P))
                      U_eff_y1m0_2 <- cbind(U_alpha,U_eta_y1m0,U_gamma_eff_y1m0_2)
                      VAR_U_eff_y1m0_2 <- var(U_eff_y1m0_2)
                      
                      last_eff_2 <- cbind(((1-d$E)%*%col_Xym1)/length(d$E),mean(-(1-d$E)))
                      last_3col_eff_2 <- c(mean(-z2),last_eff_2)
                      
                      EU_eff_y1m0_2 <- rbind(c(-1,rep(0,ncol(col_Xym1)+1)),cbind(rep(0,ncol(col_Xym1)),
                                       ((row_Xym1%*%(-d$E*col_Xym1))/length(d$E)),rep(0,ncol(col_Xym1))),
                                       last_3col_eff_2)
                      SE.y1m0.re <- sqrt((1/length(d$E))*(solve(EU_eff_y1m0_2)%*%VAR_U_eff_y1m0_2%*%
                                    t(solve(EU_eff_y1m0_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)])
                      
                      U_gamma_eff_y0m0_2 <- (((1-d$E)/(1-d$P))*(d$Y-yfit0a))+yfit0a-y0m0.re
                      U_eff_y0m0_2 <- cbind(U_alpha,U_gamma_eff_y0m0_2)
                      VAR_U_eff_y0m0_2 <- var(U_eff_y0m0_2)
                      last_eff_y0m0_2 <- c(mean(((1-d$E)/(1-d$P)^2)*(d$Y-yfit0a)),-1)
                      EU_eff_y0m0_2 <- rbind(c(-1,0),last_eff_y0m0_2)
                      SE.y0m0.re <- sqrt((1/length(d$E))*(solve(EU_eff_y0m0_2)%*%VAR_U_eff_y0m0_2%*%
                                    t(solve(EU_eff_y0m0_2)))[nrow(EU_eff_y0m0_2),ncol(EU_eff_y0m0_2)])
                      
                      U_gamma_eff_y1m1_2 <- (((d$E)/(d$P))*(d$Y-yfit1a))+yfit1a-y1m1.re
                      U_eff_y1m1_2 <- cbind(U_alpha,U_gamma_eff_y1m1_2)
                      VAR_U_eff_y1m1_2 <- var(U_eff_y1m1_2)
                      last_eff_y1m1_2 <- c(-mean(((d$E)/(d$P)^2)*(d$Y-yfit1a)),-1)
                      EU_eff_y1m1_2 <- rbind(c(-1,0),last_eff_y1m1_2)
                      SE.y1m1.re <- sqrt((1/length(d$E))*(solve(EU_eff_y1m1_2)%*%VAR_U_eff_y1m1_2%*%
                                    t(solve(EU_eff_y1m1_2)))[nrow(EU_eff_y1m1_2),ncol(EU_eff_y1m1_2)])
                      
                      U_eff_direct_2 <- cbind(U_eff_y1m0_2,U_gamma_eff_y0m0_2)
                      VAR_U_eff_direct_2 <- var(U_eff_direct_2)
                      EU_eff_direct_2 <- rbind(cbind(EU_eff_y1m0_2,(rep(0,nrow(EU_eff_y1m0_2)))),
                                         c(EU_eff_y0m0_2[-1,1],rep(0,(nrow(EU_eff_y1m0_2)-1)),
                                         EU_eff_y0m0_2[-1,-1]))
                      SE.direct.re <- sqrt(((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)]) +
                                      ((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[ncol(U_eff_direct_2),ncol(U_eff_direct_2)])-
                                      2*((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[nrow(EU_eff_y1m0_2),ncol(U_eff_direct_2)]))
                      
                      U_eff_indirect_2 <- cbind(U_eff_y1m0_2,U_gamma_eff_y1m1_2)
                      VAR_U_eff_indirect_2 <- var(U_eff_indirect_2)
                      EU_eff_indirect_2 <- rbind(cbind(EU_eff_y1m0_2,(rep(0,nrow(EU_eff_y1m0_2)))),
                                           c(EU_eff_y1m1_2[-1,1],rep(0,(nrow(EU_eff_y1m0_2)-1)),
                                           EU_eff_y1m1_2[-1,-1]))
                      SE.indirect.re <- sqrt(((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)]) +
                                        ((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[ncol(U_eff_indirect_2),ncol(U_eff_indirect_2)])-
                                        2*((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[nrow(EU_eff_y1m0_2),ncol(U_eff_indirect_2)]))
                      }
                      
                      else {
                        
                      d <- data.frame(Y=Y,M=M,E=E,Xoutc1=I(Xoutc1),Xoutc0=I(Xoutc0),Xsec1=I(Xsec1),
                           Xsec0=I(Xsec0),Xthird0=I(Xthird0))
                        
                      ## Restricted MLE:
                      ymfit1 <- glm(Y~Xoutc1+M,data=d[d$E==1,],family="binomial")
                      ymfit0 <- glm(Y~Xoutc0+M,data=d[d$E==0,],family="binomial")
                      d$pred1 <- predict(ymfit1,newdata=d,type="response")
                      d$pred0 <- predict(ymfit0,newdata=d,type="response")
                        
                      y1m0.rmle <- mean(predict(ymfit1,newdata=d[d$E==0,],type='response'))
                      y0m0.rmle <- mean(predict(ymfit0,newdata=d[d$E==0,],type='response'))
                      y1m1.rmle <- mean(predict(ymfit1,newdata=d[d$E==1,],type='response'))
                      direct.rmle <- mean(predict(ymfit1,newdata=d[d$E==0,],type='response')) - 
                                     mean(predict(ymfit0,newdata=d[d$E==0,],type='response'))
                      indirect.rmle <- mean(predict(ymfit1,newdata=d[d$E==1,],type='response')) - 
                                       mean(predict(ymfit1,newdata=d[d$E==0,],type='response'))
                        
                      ## Efficient estimators:
                      d$P <- mean(d$E)
                        
                        
                      ## Regression Model for Covariates on Outcome in the Control Arm:
                        
                      ## Approach 1:
                      yfit0 <- predict(glm(Y~Xsec0,data=d[d$E==0,],family=binomial,weights=(1/(1-d$P[d$E==0]))),
                      newdata=d,type="response")
                        
                      ## Approach 2:
                      d$outc_1 <- (Y*(1-d$E)/(1-d$P)) + ((d$E-d$P)/(1-d$P))*y0m0.rmle
                      d$x_vect <- t(t(Xsec0)-apply(Xsec0,2,mean))*((d$P-d$E)/(1-d$P))
                      lm0a <- as.vector(coef(lm(outc_1~x_vect,data=d)))
                      yfit0a <- apply((t(Xsec0)-apply(Xsec0,2,mean))*lm0a[-1],2,sum)
                        
                        
                      ## Regression Model for Covariates on Outcome in the Experimental Arm:
                        
                      ## Approach 1:
                      yfit1 <- predict(lm(Y~Xsec1,data=d[d$E==1,],family=binomial,weights=(1/(d$P[d$E==1]))),
                      newdata=d,type="response")
                        
                      ## Approach 2:
                      d$outc_2 <- Y*(d$E)/(d$P) - ((d$E-d$P)/(d$P))*y1m1.rmle
                      d$x_vect2 <- t(t(Xsec1)-apply(Xsec1,2,mean))*((d$E-d$P)/(d$P))
                      lm1a <- as.vector(coef(lm(outc_2~x_vect2,data=d)))
                      yfit1a <- apply((t(Xsec1)-apply(Xsec1,2,mean))*lm1a[-1],2,sum)
                        
                        
                      ## Regress Fitted Values from Outcome Model on Covariates in Control Arm:
                        
                      ## Approach 1:
                      z1 <- predict(glm(pred1~Xthird0,data=d[d$E==0,],family='binomial',weights=(1/(1-d$P[d$E==0]))),
                      newdata=d,type="response")
                        
                      ## Approach 2:
                      d$o1 <- d$pred1*(1-d$E)/(1-d$P) + ((d$E-d$P)/(1-d$P))*y1m0.rmle
                      d$xv <- t(t(Xthird0)-apply(Xthird0,2,mean))*((d$P-d$E)/(1-d$P))
                      B <- as.vector(coef(lm(o1~xv,data=d)))
                      z2 <- apply((t(Xthird0)-apply(Xthird0,2,mean))*B[-1],2,sum)
                        
                        
                      ## Proposal I: locally efficient estimators
                      y0m0.le <- mean(yfit0)
                      y1m1.le <- mean(yfit1)
                      y1m0.le <- mean(z1)
                      direct.le <- mean(z1) - mean(yfit0)
                      indirect.le <- mean(yfit1) - mean(z1)
                        
                        
                      ## Proposal II: restricted efficient estimators
                      y0m0.re <- mean((1/(1-d$P))*(((1-d$E)*Y)-(yfit0a*(d$P-d$E))))
                      y1m1.re <- mean((1/(d$P))*(((d$E)*Y)-(yfit1a*(d$E-d$P))))
                      y1m0.re <- mean((1/(1-d$P))*(((1-d$E)*d$pred1)-(z2*(d$P-d$E))))
                      direct.re <- y1m0.re - y0m0.re
                      indirect.re <- y1m1.re - y1m0.re
                        
                        
                      ## Standard errors restricted MLE:
                      col_Xym1 <- cbind(1,d$M,Xoutc1)
                      col_Xym0 <- cbind(1,d$M,Xoutc0)
                        
                      row_Xym1 <- rbind(1,d$M,t(Xoutc1))
                      row_Xym0 <- rbind(1,d$M,t(Xoutc0))
                        
                      # Constructing U:
                      r_y1m0 <- d$Y-d$pred1
                      U_eta_y1m0 <- (d$E*col_Xym1)*r_y1m0
                      U_gamma_y1m0 <- (d$pred1 - y1m0.rmle)*(1-d$E)
                      U_y1m0 <- cbind(U_eta_y1m0,U_gamma_y1m0)
                        
                      r_y0m0 <- d$Y-d$pred0
                      U_eta_y0m0 <- ((1-d$E)*col_Xym0)*r_y0m0
                      U_gamma_y0m0 <- (d$pred0 - y0m0.rmle)*(1-d$E)
                      U_y0m0 <- cbind(U_eta_y0m0,U_gamma_y0m0)
                        
                      r_y1m1 <- d$Y-d$pred1
                      U_eta_y1m1 <- (d$E*col_Xym1)*r_y1m1
                      U_gamma_y1m1 <- (d$pred1 - y1m1.rmle)*(d$E)
                      U_y1m1 <- cbind(U_eta_y1m1,U_gamma_y1m1)
                        
                      # Constructing dU:
                      e_y1m0 <- d$pred1/(1-d$pred1)
                      e_y1m1 <- d$pred1/(1-d$pred1)
                      e_y0m0 <- d$pred0/(1-d$pred0)
                      
                      exp_y1m0 <- e_y1m0/(1+e_y1m0)^2
                      exp_y1m1 <- e_y1m1/(1+e_y1m1)^2
                      exp_y0m0 <- e_y0m0/(1+e_y0m0)^2
                      
                      last_row10 <- cbind(col_Xym1*exp_y1m0,rep(-1,length(d$E)))
                      last_row00 <- cbind(col_Xym0*exp_y0m0,rep(-1,length(d$E)))
                      last_row11 <- cbind(col_Xym1*exp_y1m1,rep(-1,length(d$E)))
                      
                      EU_y1m0 <- rbind((cbind(((t(t(row_Xym1)*exp_y1m0)%*%(-d$E*col_Xym1))/length(d$E)),
                                 rep(0,nrow(row_Xym1)))),apply((1-d$E)*last_row10,2,mean))
                      EU_y0m0 <- rbind((cbind(((t(t(row_Xym0)*exp_y0m0)%*%(-(1-d$E)*col_Xym0))/length(d$E)),
                                 rep(0,nrow(row_Xym0)))),apply((1-d$E)*last_row00,2,mean))
                      EU_y1m1 <- rbind((cbind(((t(t(row_Xym1)*exp_y1m1)%*%(-d$E*col_Xym1))/length(d$E)),
                                 rep(0,nrow(row_Xym1)))),apply(d$E*last_row11,2,mean))
                        
                      U_direct <- cbind(U_y1m0,U_y0m0)
                      U_indirect <- cbind(U_y1m1,U_y1m0)
                        
                      mat0_y0m0 <- matrix(rep(0,nrow(EU_y0m0)*nrow(EU_y1m0)),ncol=nrow(EU_y1m0))
                      mat0_y1m0_1 <- matrix(rep(0,nrow(EU_y1m0)*nrow(EU_y0m0)),ncol=nrow(EU_y0m0))
                      mat0_y1m0_2 <- matrix(rep(0,nrow(EU_y1m0)*nrow(EU_y1m1)),ncol=nrow(EU_y1m1))
                      mat0_y1m1 <- matrix(rep(0,nrow(EU_y1m1)*nrow(EU_y1m0)),ncol=nrow(EU_y1m1))
                      EU_direct_p1 <- rbind(EU_y1m0,mat0_y0m0)
                      EU_direct_p2 <- rbind(mat0_y1m0_1,EU_y0m0)
                      EU_direct <- cbind(EU_direct_p1,EU_direct_p2)
                      EU_indirect_p1 <- rbind(EU_y1m1,mat0_y1m0_2)
                      EU_indirect_p2 <- rbind(mat0_y1m1,EU_y1m0)
                      EU_indirect <- cbind(EU_indirect_p1,EU_indirect_p2)
                        
                      # Variance U:
                      VAR_U_y1m0 <- var(U_y1m0)
                      VAR_U_y0m0 <- var(U_y0m0)
                      VAR_U_y1m1 <- var(U_y1m1)
                        
                      VAR_U_direct <- var(U_direct)
                      VAR_U_indirect <- var(U_indirect)
                        
                      # SE's:
                      SE.y1m0.rmle <- sqrt((1/length(d$E))*(solve(EU_y1m0)%*%VAR_U_y1m0%*%
                                      t(solve(EU_y1m0)))[nrow(EU_y1m0),nrow(EU_y1m0)])
                      SE.y0m0.rmle <- sqrt((1/length(d$E))*(solve(EU_y0m0)%*%VAR_U_y0m0%*%
                                      t(solve(EU_y0m0)))[nrow(EU_y0m0),nrow(EU_y0m0)])
                      SE.y1m1.rmle <- sqrt((1/length(d$E))*(solve(EU_y1m1)%*%VAR_U_y1m1%*%
                                      t(solve(EU_y1m1)))[nrow(EU_y1m1),nrow(EU_y1m1)])
                        
                      SE.direct.rmle <- sqrt((1/length(d$E))*((solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0),nrow(EU_y1m0)]) +
                                        ((1/length(d$E))*(solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0)+nrow(EU_y0m0),nrow(EU_y1m0)+nrow(EU_y0m0)])-
                                        2*((1/length(d$E))*(solve(EU_direct)%*%VAR_U_direct%*%
                                        t(solve(EU_direct)))[nrow(EU_y1m0),nrow(EU_y1m0)+nrow(EU_y0m0)]))
                      SE.indirect.rmle <- sqrt((1/length(d$E))*((solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1),nrow(EU_y1m1)]) +
                                          ((1/length(d$E))*(solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1)+nrow(EU_y1m0),nrow(EU_y1m1)+nrow(EU_y1m0)])-
                                          2*((1/length(d$E))*(solve(EU_indirect)%*%VAR_U_indirect%*%
                                          t(solve(EU_indirect)))[nrow(EU_y1m1),nrow(EU_y1m1)+nrow(EU_y1m0)]))
                        
                        
                      ## Standard errors LE:
                      U_alpha <- d$E - d$P
                      U_gamma_eff_y1m0 <- ((d$pred1 - y1m0.le)*(1-d$E)) + (z1*(d$E-d$P))
                      U_eff_y1m0 <- cbind(U_alpha,U_eta_y1m0,U_gamma_eff_y1m0)
                      VAR_U_eff_y1m0 <- var(U_eff_y1m0)
                        
                      last_eff <- cbind(((1-d$E)%*%(col_Xym1*exp_y1m0))/length(d$E),mean(-(1-d$E)))
                      last_3col_eff <- cbind(mean(-z1),last_eff)
                        
                      EU_eff_y1m0 <- rbind(c(-1,rep(0,ncol(col_Xym1)+1)),cbind(rep(0,ncol(col_Xym1)),
                                     ((row_Xym1%*%(-d$E*col_Xym1*exp_y1m0))/length(d$E)),rep(0,ncol(col_Xym1))),
                                     last_3col_eff)
                      SE.y1m0.le <- sqrt((1/length(d$E))*(solve(EU_eff_y1m0)%*%VAR_U_eff_y1m0%*%
                                    t(solve(EU_eff_y1m0)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)])
                        
                      U_gamma_eff_y0m0 <- (((1-d$E)/(1-d$P))*(d$Y-yfit0))+yfit0-y0m0.le
                      U_eff_y0m0 <- cbind(U_alpha,U_gamma_eff_y0m0)
                      VAR_U_eff_y0m0 <- var(U_eff_y0m0)
                        
                      last_eff_y0m0 <- c(mean(((1-d$E)/(1-d$P)^2)*(d$Y-yfit0)),-1)
                      EU_eff_y0m0 <- rbind(c(-1,0),last_eff_y0m0)
                      SE.y0m0.le <- sqrt((1/length(d$E))*(solve(EU_eff_y0m0)%*%VAR_U_eff_y0m0%*%
                                    t(solve(EU_eff_y0m0)))[nrow(EU_eff_y0m0),ncol(EU_eff_y0m0)])
                        
                      U_gamma_eff_y1m1 <- (((d$E)/(d$P))*(d$Y-yfit1))+yfit1-y1m1.le
                      U_eff_y1m1 <- cbind(U_alpha,U_gamma_eff_y1m1)
                      VAR_U_eff_y1m1 <- var(U_eff_y1m1)
                        
                      last_eff_y1m1 <- c(-mean(((d$E)/(d$P)^2)*(d$Y-yfit1)),-1)
                      EU_eff_y1m1 <- rbind(c(-1,0),last_eff_y1m1)
                      SE.y1m1.le <- sqrt((1/length(d$E))*(solve(EU_eff_y1m1)%*%VAR_U_eff_y1m1%*%
                                    t(solve(EU_eff_y1m1)))[nrow(EU_eff_y1m1),ncol(EU_eff_y1m1)])
                        
                      U_eff_direct <- cbind(U_eff_y1m0,U_gamma_eff_y0m0)
                      VAR_U_eff_direct <- var(U_eff_direct)
                      EU_eff_direct <- rbind(cbind(EU_eff_y1m0,(rep(0,nrow(EU_eff_y1m0)))),
                                       c(EU_eff_y0m0[-1,1],rep(0,(nrow(EU_eff_y1m0)-1)),
                                       EU_eff_y0m0[-1,-1]))
                      SE.direct.le <- sqrt(((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)]) +
                                      ((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[ncol(U_eff_direct),ncol(U_eff_direct)])-
                                      2*((1/length(d$E))*(solve(EU_eff_direct)%*%VAR_U_eff_direct%*%
                                      t(solve(EU_eff_direct)))[nrow(EU_eff_y1m0),ncol(U_eff_direct)]))
                        
                      U_eff_indirect <- cbind(U_eff_y1m0,U_gamma_eff_y1m1)
                      VAR_U_eff_indirect <- var(U_eff_indirect)
                      EU_eff_indirect <- rbind(cbind(EU_eff_y1m0,(rep(0,nrow(EU_eff_y1m0)))),
                                         c(EU_eff_y1m1[-1,1],rep(0,(nrow(EU_eff_y1m0)-1)),
                                         EU_eff_y1m1[-1,-1]))
                      SE.indirect.le <- sqrt(((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[nrow(EU_eff_y1m0),ncol(EU_eff_y1m0)]) +
                                        ((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[ncol(U_eff_indirect),ncol(U_eff_indirect)])-
                                        2*((1/length(d$E))*(solve(EU_eff_indirect)%*%VAR_U_eff_indirect%*%
                                        t(solve(EU_eff_indirect)))[nrow(EU_eff_y1m0),ncol(U_eff_indirect)]))
                        
                        
                      ## Standard errors RE:
                      U_gamma_eff_y1m0_2 <- ((d$pred1 - y1m0.re)*(1-d$E)) + (z2*(d$E-d$P))
                      U_eff_y1m0_2 <- cbind(U_alpha,U_eta_y1m0,U_gamma_eff_y1m0_2)
                      VAR_U_eff_y1m0_2 <- var(U_eff_y1m0_2)
                        
                      last_eff_2 <- cbind(((1-d$E)%*%(col_Xym1*exp_y1m0))/length(d$E),mean(-(1-d$E)))
                      last_3col_eff_2 <- c(mean(-z2),last_eff_2)
                        
                      EU_eff_y1m0_2 <- rbind(c(-1,rep(0,ncol(col_Xym1)+1)),cbind(rep(0,ncol(col_Xym1)),
                                       ((row_Xym1%*%(-d$E*col_Xym1*exp_y1m0))/length(d$E)),rep(0,ncol(col_Xym1))),
                                       last_3col_eff_2)
                      SE.y1m0.re <- sqrt((1/length(d$E))*(solve(EU_eff_y1m0_2)%*%VAR_U_eff_y1m0_2%*%
                                    t(solve(EU_eff_y1m0_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)])
                        
                      U_gamma_eff_y0m0_2 <- (((1-d$E)/(1-d$P))*(d$Y-yfit0a))+yfit0a-y0m0.re
                      U_eff_y0m0_2 <- cbind(U_alpha,U_gamma_eff_y0m0_2)
                      VAR_U_eff_y0m0_2 <- var(U_eff_y0m0_2)
                      last_eff_y0m0_2 <- c(mean(((1-d$E)/(1-d$P)^2)*(d$Y-yfit0a)),-1)
                      EU_eff_y0m0_2 <- rbind(c(-1,0),last_eff_y0m0_2)
                      SE.y0m0.re <- sqrt((1/length(d$E))*(solve(EU_eff_y0m0_2)%*%VAR_U_eff_y0m0_2%*%
                                    t(solve(EU_eff_y0m0_2)))[nrow(EU_eff_y0m0_2),ncol(EU_eff_y0m0_2)])
                        
                      U_gamma_eff_y1m1_2 <- (((d$E)/(d$P))*(d$Y-yfit1a))+yfit1a-y1m1.re
                      U_eff_y1m1_2 <- cbind(U_alpha,U_gamma_eff_y1m1_2)
                      VAR_U_eff_y1m1_2 <- var(U_eff_y1m1_2)
                      last_eff_y1m1_2 <- c(-mean(((d$E)/(d$P)^2)*(d$Y-yfit1a)),-1)
                      EU_eff_y1m1_2 <- rbind(c(-1,0),last_eff_y1m1_2)
                      SE.y1m1.re <- sqrt((1/length(d$E))*(solve(EU_eff_y1m1_2)%*%VAR_U_eff_y1m1_2%*%
                                    t(solve(EU_eff_y1m1_2)))[nrow(EU_eff_y1m1_2),ncol(EU_eff_y1m1_2)])
                        
                      U_eff_direct_2 <- cbind(U_eff_y1m0_2,U_gamma_eff_y0m0_2)
                      VAR_U_eff_direct_2 <- var(U_eff_direct_2)
                      EU_eff_direct_2 <- rbind(cbind(EU_eff_y1m0_2,(rep(0,nrow(EU_eff_y1m0_2)))),
                                         c(EU_eff_y0m0_2[-1,1],rep(0,(nrow(EU_eff_y1m0_2)-1)),
                                         EU_eff_y0m0_2[-1,-1]))
                      SE.direct.re <- sqrt(((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)]) +
                                      ((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[ncol(U_eff_direct_2),ncol(U_eff_direct_2)])-
                                      2*((1/length(d$E))*(solve(EU_eff_direct_2)%*%VAR_U_eff_direct_2%*%
                                      t(solve(EU_eff_direct_2)))[nrow(EU_eff_y1m0_2),ncol(U_eff_direct_2)]))
                        
                      U_eff_indirect_2 <- cbind(U_eff_y1m0_2,U_gamma_eff_y1m1_2)
                      VAR_U_eff_indirect_2 <- var(U_eff_indirect_2)
                      EU_eff_indirect_2 <- rbind(cbind(EU_eff_y1m0_2,(rep(0,nrow(EU_eff_y1m0_2)))),
                                           c(EU_eff_y1m1_2[-1,1],rep(0,(nrow(EU_eff_y1m0_2)-1)),
                                           EU_eff_y1m1_2[-1,-1]))
                      SE.indirect.re <- sqrt(((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[nrow(EU_eff_y1m0_2),ncol(EU_eff_y1m0_2)]) +
                                        ((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[ncol(U_eff_indirect_2),ncol(U_eff_indirect_2)])-
                                        2*((1/length(d$E))*(solve(EU_eff_indirect_2)%*%VAR_U_eff_indirect_2%*%
                                        t(solve(EU_eff_indirect_2)))[nrow(EU_eff_y1m0_2),ncol(U_eff_indirect_2)]))
                      }
                      
                      results <- data.frame(cbind(c(y1m0.rmle,y0m0.rmle,y1m1.rmle,direct.rmle,indirect.rmle,y1m0.le,
                                                    y0m0.le,y1m1.le,direct.le,indirect.le,y1m0.re,y0m0.re,y1m1.re,
                                                    direct.re,indirect.re),
                                                  c(SE.y1m0.rmle,SE.y0m0.rmle,SE.y1m1.rmle,SE.direct.rmle,
                                                    SE.indirect.rmle,SE.y1m0.le,SE.y0m0.le,SE.y1m1.le,SE.direct.le,
                                                    SE.indirect.le,SE.y1m0.re,SE.y0m0.re,SE.y1m1.re,SE.direct.re,
                                                    SE.indirect.re)))
                      rownames(results) <- c('RMLE Y1M0','RMLE Y0M0','RMLE Y1M1','RMLE Direct','RMLE Indirect',
                                             'LE Y1M0','LE Y0M0','LE Y1M1','LE Direct','LE Indirect',
                                             'RE Y1M0','RE Y0M0','RE Y1M1','RE Direct','RE Indirect')
                      colnames(results) <- c('Estimate','SE')
                      
                      results
               }


### Try function ###
eff.mediate(Y,M,E,Xoutc1,Xoutc0,Xsec1,Xsec0,Xthird0,typ='continuous')


### Binary Outcome ###
full <- read.csv("C:\\Users\\sjvdnber\\Documents\\PHD\\Simulations\\Linear-Logistic\\prevention.csv",head=T,sep=",")
data <- as.data.frame(full)

data <- data[data$SCR01!=8,]

data$ca1 <- ifelse(data$cad1==1,1,0)
data$ca2 <- ifelse(data$cad1==2,1,0)
data$ca3 <- ifelse(data$cad1==3,1,0)
data$ca4 <- ifelse(data$cad1==4,1,0)
data$site2 <- ifelse(data$site==2,1,0)
data$site3 <- ifelse(data$site==3,1,0)

data$binhamda <- ifelse(data$hamda>=20,1,0)

### Specifying variables and models ###

Y <- data$binhamda # outcome
M <- data$amedx # mediator
E <- data$interven # exposure

# Predictors outcome model in the treated (first model)
Xoutc1 <- cbind(data$ca1,data$ca2,data$ca3,data$ca4,data$hamda1,data$ssix01,data$SCR01,data$site2,data$site3) 

# Predictors outcome model in the untreated (first model)
Xoutc0 <- cbind(data$ca1,data$ca2,data$ca3,data$ca4,data$hamda1,data$ssix01,data$SCR01,data$site2,data$site3)

# Predictors third model in the untreated
Xthird0 <- cbind(data$ca1,data$ca2,data$ca3,data$ca4,data$hamda1,data$ssix01,data$SCR01,data$site2,data$site3)

# Predictors second model in the treated
Xsec1 <- cbind(data$ca1,data$ca2,data$ca3,data$ca4,data$hamda1,data$ssix01,data$SCR01,data$site2,data$site3)

# Predictors second model in the untreated
Xsec0 <- cbind(data$ca1,data$ca2,data$ca3,data$ca4,data$hamda1,data$ssix01,data$SCR01,data$site2,data$site3)


### Try function ###
eff.mediate(Y,M,E,Xoutc1,Xoutc0,Xsec1,Xsec0,Xthird0,typ='binary')
