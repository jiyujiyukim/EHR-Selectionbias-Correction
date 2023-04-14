library(MASS)
require(pbivnorm)
library("numDeriv")
library(betareg) 
library(sampleSelection)
library(dplyr)

b0<-10 
Py<- 0.1

set.seed(2021)
N<-100000; K<-2000
Ps<- 0.1; 

x1 <- rnorm(N, 0,1)
x2 <-rnorm(N, 0, 1)
z3<- rnorm(N, 0, 1)

#setting
rho <- 0.3
a0 <- -1.55
a1<-a2<-a3<- 0.5
b1<- b2 <- 0.5
#a0<- -2 # -1.55 for 10%; -2 for 5% 

area <- c(1:K)
area.pop<- rep(N/K, length(area))
id.area<- rep(c(1:K), each=N/K)

#read the file 
dat<- read.csv( ".../dat.csv")
dat <- dat[,-c(1)]
  
  ###### use different EHR S #######################################
  #aggregated 
  cen<- aggregate(dat, by=list(dat$id.area), FUN=mean)
  cen<- data.frame(cen, area.pop)
  names(cen)[names(cen)=="Freq"]<-"area.pop"
  sample.n<- as.vector(aggregate(dat, by=list(dat$id.area), FUN = sum)[8]) #S 
  colnames(sample.n)<-"sample.n"
  cen<- cbind(cen, sample.n)
  colnames(cen)[1]<-"zone"
  
  #EHR
  EHR<- subset(dat, S==1)
  
  #survey data 
  Ref1<- dat[which(dat$S==0),]
  Ref<- Ref1[which(Ref1$R==1),]
  
  tmp<- cen
  EHR<- merge(tmp[,c("id.area","x1", "x2")],EHR, by="id.area",all.y =T)
  head(EHR)
  colnames(EHR)<- c("zone", "x1.agg", "z2.agg","id", "x1", "x2", "Yp", "Sp", "Y", "S")
  
  ######## Estimate Two stage model ########
  ###Two stage model when individual inf available --- Ideal 
  #step 1. estimate alphas
  a0.hat.star<-coef(summary(glm(S~ x1+x2, family=binomial(link="probit"), dat)))[1,1]
  a1.hat.star<-coef(summary(glm(S~ x1+x2, family=binomial(link="probit"), dat)))[2,1]
  a2.hat.star<-coef(summary(glm(S~ x1+x2, family=binomial(link="probit"), dat)))[3,1]
  
  #step 2. calculate lambda for each individual in EHR 
  ld1<- a0.hat.star+a1.hat.star*EHR$x1+a2.hat.star*EHR$x2   # Best situation that we all observe individual level 
  nume<- dnorm(ld1,0,1); denume<- pnorm(ld1,0,1)
  EHR$lambda1.st<-nume/denume
  
  #step 3. estimate beta
  b0.hat.star1<-coef(summary(lm(Y~x1+lambda1.st,data=EHR)))[1]
  b1.hat.star1<- coef(summary(lm(Y~x1+lambda1.st,data=EHR)))[2]
  rho1.star<-coef(summary(lm(Y~x1+lambda1.st, data=EHR)))[3]
  
  ###or Using selectionSample package ###
  heckvan = heckit( S~ x1+x2, Y~x1, data=dat )
  a0.hat.star<-coef(heckvan)[1]
  a1.hat.star<-coef(heckvan)[2]
  a2.hat.star<-coef(heckvan)[3]
  b0.hat.star1<-coef(heckvan)[4]
  b1.hat.star1<-coef(heckvan)[5]
  rho.twost.ideal <- coef(heckvan)[8]

  #se
  g<-summary(glm(S~ x1+x2, family=binomial(link="probit"), dat))
  #############################Variance of beta ###############################
  #sigma.e.hat
  sv<- EHR$Y-b0.hat.star1-EHR$x1*b1.hat.star1-EHR$lambda1.st*rho1.star #n*1
  sv.t<- t(sv) #1*n
  sv<- t(sv.t)
  d<- EHR$lambda1.st*(EHR$lambda1.st+ (a0.hat.star+a1.hat.star*EHR$x1+a2.hat.star*EHR$x2)) #
  e.hat<-sv.t%*%sv/nrow(EHR)+ sum(d)/nrow(EHR)*rho1.star^2 #
  rho.hat<-sqrt(rho1.star^2/(e.hat)) ###rho hat! 
  rho.hat<-data.frame(rho.hat)[1,1]
  
  #varaince beta 
  V<- data.frame(vcov(g))
  V<- data.matrix(V)
  
  x.st<- cbind(1, EHR$x1, EHR$lambda1.st) #n*2 --> n*3
  x.st.t<- t(x.st) #2*n --> n*3
  d<- EHR$lambda1.st*(EHR$lambda1.st+ (a0.hat.star+a1.hat.star*EHR$x1+a2.hat.star*EHR$x2))
  del<- diag(d, nrow(EHR), nrow(EHR)) #n*n
  iden<- diag(nrow(EHR)) #n*n
  z<- cbind(1,EHR$x1, EHR$x2) # n*number of z 
  Q<- rho.hat^2*(x.st.t%*%del%*%z)%*%V%*%t((x.st.t%*%del%*%z))
  Varaince.beta<-solve(x.st.t%*%x.st)%*%(x.st.t%*%(iden-rho.hat^2*del)%*%x.st+Q)%*%solve(x.st.t%*%x.st)
  b1.hat.star1.se<-sqrt(data.frame(Varaince.beta)[2,2]*e.hat)
  ###############################################################################  
  b1.hat.star1.se <- b1.hat.star1.se[1,1]


  #################
  ###Two stage model when individual inf non available --- census tract 
  #step 1. Estimate alphas.--> Berkson's minimum chi-squre 
  cen$w <- (1/cen$area.pop)*cen$S*(1-cen$S)/dnorm(qnorm(cen$S))^2
  m <- lm(qnorm(cen$S)~ cen$x1 + cen$x2, weights=1/(cen$w))
  a0.hat<-coef(summary(m))[1,1]
  a1.hat<-coef(summary(m))[2,1]
  a2.hat<-coef(summary(m))[3,1]
  
  #step 2. Calculate lambda for each individual
  ld1<- a0.hat+a1.hat*EHR$x1+a2.hat*EHR$x2   # Best situation that we all observe individual level in EHR
  ld2<- a0.hat+a1.hat*EHR$x1+a2.hat*EHR$z2.agg  # Option 1. we use individual + agg
  
  nume<- dnorm(ld1,0,1); denume<- pnorm(ld1,0,1)
  EHR$lambda1<-nume/denume
  nume<-  dnorm(ld2,0,1); denume<- pnorm(ld2,0,1)
  EHR$lambda2<-nume/denume
  
  #step3. regress x and lamda "EHR data" 
  b0.hat1<-coef(summary(lm(Y~x1+lambda1, data=EHR)))[1]
  b0.hat1.agg<-coef(summary(lm(Y~x1+lambda2, data=EHR)))[1]
  
  b1.hat1<-coef(summary(lm(Y~x1+lambda1, data=EHR)))[2]
  b1.hat1.agg<-coef(summary(lm(Y~x1+lambda2, data=EHR)))[2]
  
  rho1<-coef(summary(lm(Y~x1+lambda1, data=EHR)))[3]
  rho1.agg<-coef(summary(lm(Y~x1+lambda2, data=EHR)))[3]
  
  #############################Variance of beta ###############################
  #sigma.e.hat
  sv<- EHR$Y-b0.hat1-EHR$x1*b1.hat1-EHR$lambda1*rho1 #
  sv.t<- t(sv)
  sv<- t(sv.t)
  d<- EHR$lambda1*(EHR$lambda1+ (a0.hat+a1.hat*EHR$x1+a2.hat*EHR$x2)) #
  e.hat<-sv.t%*%sv/nrow(EHR)+ sum(d)/nrow(EHR)*rho1^2 #  sigma^2 
  rho.hat<-rho1/sqrt(e.hat) ###rho hat! 
  rho.hat<-rho.hat[1,1]
  rho.twost.cen<- rho.hat
  
  #varaince beta 
  V<- data.frame(vcov(m))
  V<- data.matrix(V)
  
  x.st<- cbind(1, EHR$x1, EHR$lambda1) #n*2 --> n*3
  x.st.t<- t(x.st) #2*n --> n*3
  d<- EHR$lambda1*(EHR$lambda1+  (a0.hat+a1.hat*EHR$x1+a2.hat*EHR$x2))
  del<- diag(d, nrow(EHR), nrow(EHR)) #n*n
  iden<- diag(nrow(EHR)) #n*n
  z<- cbind(1,EHR$x1, EHR$x2) # n*number of z 
  Q<- rho.hat^2*(x.st.t%*%del%*%z)%*%V%*%t((x.st.t%*%del%*%z))
  Varaince.beta<-solve(x.st.t%*%x.st)%*%(x.st.t%*%(iden-rho.hat^2*del)%*%x.st+Q)%*%solve(x.st.t%*%x.st)
  b1.hat1.se<-sqrt(data.frame(Varaince.beta)[2,2]*e.hat)
  b1.hat1.se<- b1.hat1.se[1,1]
  ###############################################################################  
  
  #############################Variance of beta ###############################
  #sigma.e.hat
  sv<- EHR$Y-b0.hat1.agg-EHR$x1*b1.hat1.agg-EHR$lambda2*rho1.agg #
  sv.t<- t(sv)
  sv<- t(sv.t)
  d<- EHR$lambda2*(EHR$lambda2+ (a0.hat+a1.hat*EHR$x1+a2.hat*EHR$z2.agg)) #
  e.hat<-sv.t%*%sv/nrow(EHR)+ sum(d)/nrow(EHR)*rho1.agg^2 #  sigma^2 
  rho.hat<-rho1.agg/sqrt(e.hat) ###rho hat! 
  rho.hat<-rho.hat[1,1]
  rho.twost.cen.agg<- rho.hat
  
  #varaince beta 
  V<- data.frame(vcov(m))
  V<- data.matrix(V)
  
  x.st<- cbind(1, EHR$x1, EHR$lambda2) #n*2 --> n*3
  x.st.t<- t(x.st) #2*n --> n*3
  d<- EHR$lambda2*(EHR$lambda2+  (a0.hat+a1.hat*EHR$x1+a2.hat*EHR$z2.agg))
  del<- diag(d, nrow(EHR), nrow(EHR)) #n*n
  iden<- diag(nrow(EHR)) #n*n
  z<- cbind(1,EHR$x1, EHR$z2.agg) # n*number of z 
  Q<- rho.hat^2*(x.st.t%*%del%*%z)%*%V%*%t((x.st.t%*%del%*%z))
  Varaince.beta<-solve(x.st.t%*%x.st)%*%(x.st.t%*%(iden-rho.hat^2*del)%*%x.st+Q)%*%solve(x.st.t%*%x.st)
  b1.hat1.se.agg<-sqrt(data.frame(Varaince.beta)[2,2]*e.hat)
  b1.hat1.se.agg<- b1.hat1.se.agg[1,1]
  ###############################################################################  

 
  ####Use reference sample ############################################
  dat$p.ideal <-dat$delta  #delta or delta.p
  dat$p.ideal<- ifelse(dat$p.ideal==0.000, 0.00001, dat$p.ideal)
  dat$wt.ideal <-1/dat$p.ideal
  
  scale.w<-sum(dat$S==1)/sum(dat$wt.ideal[dat$S==1])  ##
  dat$wt.ideal<- dat$wt.ideal*scale.w
  head(dat)
  
  #EHR
  EHR<- dat[which(dat$S==1),]
  
  head(Ref); head(EHR)
  EHR$wt.r<- 1; EHR$comb<- 1; Ref$comb<-0
  EHR1<- EHR[, c("id","x1","x2","Yp","Sp","Y","S","id.area","delta","delta.p","R","p.ref","wt.r","comb")]
  head(Ref); head(EHR1)
  comb<- rbind(Ref, EHR1)
  
  head(comb)
  two.ref.a0<- coef(summary(glm(comb~ x1+x2, family=binomial(link="probit"), weights = wt.r, data=comb)))[1,1]
  two.ref.a1<- coef(summary(glm(comb~ x1+x2, family=binomial(link="probit"), weights = wt.r, data=comb)))[2,1]
  two.ref.a2<- coef(summary(glm(comb~ x1+x2, family=binomial(link="probit"), weights = wt.r, data=comb)))[3,1]
  
  ld1.ref<-  two.ref.a0+  two.ref.a1*EHR$x1+  two.ref.a2*EHR$x2 
  nume<- dnorm(ld1.ref,0,1); denume<- pnorm(ld1.ref,0,1)
  EHR$lambda.ref<-nume/denume
  
  #step 3. regress x and lamda 
  b0.hat1.ref<- coef(summary(lm(Y~x1+lambda.ref, data=EHR)))[1]
  b1.hat1.ref<- coef(summary(lm(Y~x1+lambda.ref, data=EHR)))[2] 
  rho.hat1.ref<- coef(summary(lm(Y~x1+lambda.ref, data=EHR)))[3] 
  
  #se 
  g<-summary(glm(comb~ x1+x2, family=binomial(link="probit"), weights = wt.r, data=comb))
  #############################Variance of beta ###############################
  #sigma.e.hat
  sv<- EHR$Y-b0.hat1.ref-EHR$x1*b1.hat1.ref-EHR$lambda.ref*rho.hat1.ref #
  sv.t<- t(sv)
  sv<- t(sv.t)
  d<- EHR$lambda.ref*(EHR$lambda.ref+ (two.ref.a0+two.ref.a1*EHR$x1+two.ref.a2*EHR$x2)) #
  e.hat<-sv.t%*%sv/nrow(EHR)+ sum(d)/nrow(EHR)*rho.hat1.ref^2 #  sigma^2 
  rho.hat<-rho.hat1.ref/sqrt(e.hat) ###rho hat! 
  rho.hat<-rho.hat[1,1]
  rho.ref.two<- rho.hat
  #varaince beta 
  V<- data.frame(vcov(g))
  V<- data.matrix(V)
  
  x.st<- cbind(1, EHR$x1, EHR$lambda.ref) #n*2 --> n*3
  x.st.t<- t(x.st) #2*n --> n*3
  d<- EHR$lambda.ref*(EHR$lambda.ref+  (two.ref.a0+two.ref.a1*EHR$x1+two.ref.a2*EHR$x2))
  del<- diag(d, nrow(EHR), nrow(EHR)) #n*n
  iden<- diag(nrow(EHR)) #n*n
  z<- cbind(1,EHR$x1, EHR$x2) # n*number of z 
  Q<- rho.hat^2*(x.st.t%*%del%*%z)%*%V%*%t((x.st.t%*%del%*%z))
  Varaince.beta<-solve(x.st.t%*%x.st)%*%(x.st.t%*%(iden-rho.hat^2*del)%*%x.st+Q)%*%solve(x.st.t%*%x.st)
  b1.hat1.se.ref<-sqrt(data.frame(Varaince.beta)[2,2]*e.hat)
  b1.hat1.se.ref<- b1.hat1.se.ref[1,1]
  ###############################################################################  

  
  #Prevalence estimates
  nonEHR<- subset(dat, S==0)
  tmp<- cen
  nonEHR<- merge(tmp[,c("id.area","x1", "x2")],nonEHR, by="id.area",all.y =T)
  colnames(nonEHR)<- c("zone", "x1.agg", "x2.agg","id", "x1", "x2", "Yp", "Sp", "Y", "S", "delta", "delta.p", "R", "p.ref", "wt.r")
  
  #ideal 
  ld1<- a0.hat.star+a1.hat.star*nonEHR$x1+a2.hat.star*nonEHR$x2   
  nume<- dnorm(ld1,0,1); denume<- 1-pnorm(ld1,0,1)
  nonEHR$ld<-nume/denume
  nonEHR$pp<- b0.hat.star1+ b1.hat.star1*nonEHR$x1 - data.frame(coef(heckvan)[7])[1,1]*data.frame(coef(heckvan)[8])[1,1]*nonEHR$ld

  #cen - two step, twostep z.agg
  cen.tmp<- cen
  EHR.summary<- aggregate(EHR, by=list(EHR$id.area), FUN=sum)
  colnames(EHR.summary)[1] <-"zone"
  meg<-merge(cen.tmp[,c("zone","id")], EHR.summary[,c("zone", "x1", "x2")], by="zone", all.x=T)
  head(meg); tail(meg)
  cen.tmp$x1.agg.notsample <- ((cen.tmp$area.pop)*(cen.tmp$x1) - meg$x1)/(cen.tmp$area.pop-cen.tmp$sample.n)
  cen.tmp$x2.agg.notsample <- ((cen.tmp$area.pop)*(cen.tmp$x2) - meg$x2)/(cen.tmp$area.pop-cen.tmp$sample.n)
  cen.tmp$x1.agg.notsample<- ifelse( is.na(cen.tmp$x1.agg.notsample), cen.tmp$x1, cen.tmp$x1.agg.notsample)
  cen.tmp$x2.agg.notsample<- ifelse( is.na(cen.tmp$x2.agg.notsample), cen.tmp$x2, cen.tmp$x2.agg.notsample)
  ld1<- a0.hat+a1.hat*cen.tmp$x1.agg.notsample+a2.hat*cen.tmp$x2.agg.notsample
  nume<- dnorm(ld1,0,1); denume<- 1-pnorm(ld1,0,1)
  cen.tmp$ld<-nume/denume
  census.sum <- (cen.tmp$area.pop-cen.tmp$sample.n)*(b0.hat1+ b1.hat1*cen.tmp$x1.agg.notsample - rho1* (cen.tmp$ld))
  pre3<-(sum(EHR$Y)+(sum(census.sum)))/nrow(dat)
  
  ld1<- a0.hat+a1.hat*cen.tmp$x1.agg.notsample+a2.hat*cen.tmp$x2.agg.notsample
  nume<- dnorm(ld1,0,1); denume<- 1-pnorm(ld1,0,1)
  cen.tmp$ld<-nume/denume
  census.sum <- (cen.tmp$area.pop-cen.tmp$sample.n)*(b0.hat1.agg+ b1.hat1.agg*cen.tmp$x1.agg.notsample - rho1.agg* (cen.tmp$ld))
  pre7<- (sum(EHR$Y)+(sum(census.sum)))/nrow(dat)
  

  #ref - two step and one step 
  Ref2<- Ref
  ld1.ref<-  two.ref.a0+  two.ref.a1*Ref2$x1+  two.ref.a2*Ref2$x2 
  nume<- dnorm(ld1.ref,0,1); denume<- 1-pnorm(ld1.ref,0,1)
  Ref2$ld<-nume/denume
  Ref2$pp5 <- b0.hat1.ref + b1.hat1.ref*Ref2$x1 - rho.hat1.ref*Ref2$ld

  ########### same   
  EHR$pp<- EHR$Y # 
  EHR$pp2<- EHR$Y 
  EHR$pp3<- EHR$Y; EHR$pp3.agg<- EHR$Y 
  EHR$pp4 <- EHR$Y 
  EHR$pp5 <- EHR$Y ;   EHR$pp6 <- EHR$Y 
  
  precomb<- rbind(nonEHR[,c("id", "pp")], EHR[,c("id", "pp")])
  pre1<- sum(precomb$pp)/N
  pre3<- pre3;   pre7 <- pre7
  pre5<-((nrow(EHR)*mean(EHR$pp5))+(nrow(dat)-nrow(EHR))*mean(Ref2$pp5))/N; 

  
  #########################
  #Naive
  naive.b1 <-coef(summary(lm(Y~x1, dat=EHR)))[2]
  prev.true<- sum(dat$Y)/nrow(dat)
  
  naive.prev<- sum(EHR$Y)/nrow(EHR) 
  naive.b1.se <- coef(summary(lm(Y~x1, dat=EHR)))[2,2]
  naive.prev.se <- sqrt(var(EHR$Y)/nrow(EHR))
  sele<- sum(dat$S)/N
  
  ##################IPW#######################
  #Ideal 
  head(dat)
  dat$p.ideal <-dat$delta  
  dat$wt.ideal <-1/dat$p.ideal
  
  scale.w<-sum(dat$S==1)/sum(dat$wt.ideal[dat$S==1])  ##
  dat$wt.ideal<- dat$wt.ideal*scale.w
  head(dat)
  
  #EHR
  EHR<- dat[which(dat$S==1),]
  head(EHR)
  
  #PS data analysis -- IDEAL 
  c<- suppressWarnings(lm(Y~  x1, weights = wt.ideal, data=EHR))
  ps.ideal.b1<- coef(summary(c))[2,1]
  ps.ideal.b1.se <- coef(summary(c))[2,2]
  
  #PS-- weighted mean 
  head(EHR)
  ps.prev.ideal=weighted.mean(EHR$Y, EHR$wt.ideal)  ## weired 
  
  #se-prev
  #ideal
  EHR$wt.ideal.sq<- EHR$wt.ideal^2
  ps.prev.ideal.se<- sqrt(var(EHR$Y)* (sum(EHR$wt.ideal.sq))/(sum(EHR$wt.ideal))^2)
  
  ###PS reference use ######### 1. IPW###
  head(Ref); head(EHR)
  EHR$wt.r<- 1
  EHR$comb<- 1
  Ref$comb<-0
  EHR1<- EHR[, c("id","x1","x2","Yp","Sp","Y","S","id.area","delta","delta.p","R","p.ref","wt.r","comb")]
  head(Ref); head(EHR1)
  comb<- rbind(Ref, EHR1)
  
  head(comb)
  
  a<-glm(comb~ x1+x2+Y, family="binomial", weights = wt.r, data=comb)
  comb$ipw<- predict(a, type="response")
  
  aa<-glm(comb~ x1+x2, family="binomial", weights = wt.r, data=comb)
  comb$ipw.ynot<- predict(aa, type="response")
  
  #prevalence 
  b<-glm(comb~ x1+x2+Y, family="binomial", weights = wt.r, data=comb)
  comb$ipw2<- predict(b, type="response")
  
  bb<-glm(comb~ x1+x2, family="binomial", weights = wt.r, data=comb)
  comb$ipw2.ynot<- predict(bb, type="response")
  
  #a-- beta 
  comb$wt.ipw <-1/comb$ipw  #--- y included 
  comb$wt.ipw.ynot <-1/comb$ipw.ynot  #---y not included
  
  #b-- prevalence
  comb$wt.ipw2 <-1/comb$ipw2  # --y included
  comb$wt.ipw2.ynot <-1/comb$ipw2.ynot # y not included 
  
  EHR2<- comb[which(comb$comb==1),]
  
  a<-lm(Y~ x1,  weights=wt.ipw, data=EHR2)
  ps.ipw.b1<- coef(summary(a))[2,1]

  aa<-lm(Y~ x1, weights=wt.ipw.ynot, data=EHR2)
  ps.ipw.b1.ynot<- coef(summary(aa))[2,1]
  
  ps.ipw.prev<- weighted.mean(EHR2$Y, EHR2$wt.ipw2) 
  ps.ipw.prev.ynot<- weighted.mean(EHR2$Y, EHR2$wt.ipw2.ynot)
  
  #se
  ps.ipw.b1.ynot.se <- coef(summary(aa))[2,2]
  
  #ipw
  EHR2$wt.ipw2.sq<- EHR2$wt.ipw2^2
  EHR2$wt.ipw2.ynot.sq<- EHR2$wt.ipw2.ynot^2
  
  ###PS reference use 2. PSAS ###
  ps2<- glm(comb~ x1+x2+Y, family=binomial(), data=comb)
  ps2.ynot<- glm(comb~ x1+x2, family=binomial(), data=comb)
  
  comb$psas <- predict(ps2, type="response")
  comb$psas.ynot<-predict(ps2.ynot, type="response")
  

  # y not include
  comb<-comb[order(comb$psas.ynot),]
  subgr1<- rep(1:100, each=round(nrow(comb)/100,0))
  subgr1<- subgr1[c(1: nrow(comb))]
  comb$subgr1<-subgr1
  #average psas
  ag<-aggregate(comb, by=list(comb$subgr1), FUN=mean)
  av.psas1<-rep(ag$psas.ynot, each=round(nrow(comb)/5,0))
  av.psas1<- av.psas1[c(1: nrow(comb))]
  comb$av.psas1<- av.psas1
  comb$wt.psas.ynot<- 1/comb$av.psas1
  scale.w<- sum(comb$comb==1)/sum(comb$wt.psas.ynot[comb$comb==1])
  comb$wt.psas.ynot<-  comb$wt.psas*scale.w
  
  #Probit prob
  ps3<- glm(comb~ x1+x2+Y, family="binomial", data=comb)
  ps3.ynot<- glm(comb~ x1+x2, family="binomial", data=comb)
  
  #1. yinclude
  comb<-comb[order(comb$psas),]
  subgr<- rep(1:100, each=round(nrow(comb)/100,0))
  subgr<- subgr[c(1: nrow(comb))]
  length(subgr) ;nrow(comb)
  comb$subgr<-subgr
  #average psas
  ag<-aggregate(comb, by=list(comb$subgr), FUN=mean)
  av.psas<-rep(ag$psas, each=round(nrow(comb)/5,0))
  av.psas<- av.psas[c(1: nrow(comb))]
  length(av.psas)
  comb$av.psas<- av.psas
  comb$wt.psas<- 1/comb$av.psas
  scale.w<- sum(comb$comb==1)/sum(comb$wt.psas[comb$comb==1])
  comb$wt.psas<-  comb$wt.psas*scale.w
  
  #y not included
  comb$psas2.ynot <- predict(ps3.ynot, type="response")
  comb<-comb[order(comb$psas2.ynot),]
  subgr1<- rep(1:100, each=round(nrow(comb)/100,0))
  subgr1<- subgr[c(1: nrow(comb))]
  comb$subgr1<-subgr1
  #average psas2
  ag2<-aggregate(comb, by=list(comb$subgr1), FUN=mean)
  av.psas2<-rep(ag2$psas2, each=round(nrow(comb)/5,0))
  av.psas2<- av.psas2[c(1: nrow(comb))]
  length(av.psas2)
  comb$av.psas2<- av.psas2
  comb$wt.psas2.ynot<- 1/comb$av.psas2
  scale.w<- sum(comb$comb==1)/sum(comb$wt.psas2.ynot[comb$comb==1])
  comb$wt.psas2.ynot<-  comb$wt.psas2.ynot*scale.w
  
  EHR3<- comb[which(comb$comb==1),]
  

  #y not include
  aa <- lm(Y~ x1, weights=wt.psas.ynot, data=EHR3)
  ps.psas.b1.ynot<-coef(summary(aa))[2,1]
  ps.psas.b1.ynot.se <- coef(summary(aa))[2,2]
  ps.psas.prev.ynot<- weighted.mean(EHR3$Y, EHR3$wt.psas2.ynot ) #-- y not include 
  
  
  #Obtain results 
  #############################beta###################
  #ideal
  b1.HII<- b1.hat.star1; b1.HII.se<- b1.hat.star1.se
  #census
  b1.HAI <- b1.hat1; b1.HAI.se<-  b1.hat1.se; 
  b1.HAA<- b1.hat1.agg; b1.HAA.se <- b1.hat1.se.agg
  #Survey
  b1.HSI<- b1.hat1.ref; b1.HSI.se <- b1.hat1.se.ref
  #PS
  b1.ipw.ynot<- ps.ipw.b1.ynot; ps.ipw.b1.ynot.se 
  b1.psas.ynot<- ps.psas.b1.ynot; ps.psas.b1.ynot.se
  #Naive
  naive.b1; naive.b1.se 
 
  ######################################prev  #####################
  #ideal
  prev.HII<- pre1
  #census
  prev.HAI<- pre3
  prev.HAA<- pre7 
  #Survey
  prev.HSI<- pre5
  #PS
  prev.ipw.ynot <-  ps.ipw.prev.ynot
  prev.psas.ynot<-  ps.psas.prev.ynot
  #Naive 
  naive.prev ;naive.prev.se
  

