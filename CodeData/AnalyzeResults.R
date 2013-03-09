####################################################################
#Can Start here now
##
setwd("U://Thesis//Simulation Code//Results")
Res<-read.csv("FullFinalResults.csv")[,-c(1)]



library(xtable)
library(ggplot2)
#ResFrame<-melt(Res,id=c("nu","n","Dist"),measure.var=c("HL1Error","ML2Error","ArithError","MedError"))
ResFrame<-melt(Res,id=c("nu","n","Dist","Sample"),measure.var=c("HL1Error","ML2Error","ArithError","MedError"))
colnames(ResFrame)[5:6]<-c("Estimator","Error")
levels(ResFrame$Estimator)<-c("Hartley L1","Manton L2","Projected Mean","Median")
x<-ddply(ResFrame,.(Dist,nu,n,Estimator),summarize,Median=round(median(Error),4),Bias=round(mean(Error),4),RMSE=round(sqrt(mean(Error^2)),4))


#Plot boxplots to see distribution of errors Do sample sizes seperately
ResFrame1<-ResFrame[ResFrame$nu==.25,]
m<-which.max(ResFrame1$Error)
ResFrame1<-ResFrame1[-m,]

max(ResFrame1$Error)

qplot(Estimator,Error,data=subset(ResFrame,nu==.5)[-13418,],geom="boxplot",main=expression(nu==0.50))+facet_grid(n~Dist,scales="free_y")
qplot(Estimator,Error,data=ResFrame[ResFrame$nu==.75,],geom="boxplot",main=expression(nu==0.75))+facet_grid(n~Dist,scales="free_y")
qplot(Estimator,Error,data=ResFrame1,geom="boxplot",main=expression(nu==0.25))+facet_grid(n~Dist,scales="free_y")


head(ResFrame[ResFrame$nu==.5,][-13418,])


#Get Summary stats of the error distributions
#Centered (around S_P) Median, Bias and RMSE for our four estimators
CMedian<-matrix(x$Median[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CMedian<-x$Median-kronecker(CMedian,matrix(rep(1,4),4))

CBias<-matrix(x$Bias[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CBias<-x$Bias-kronecker(CBias,matrix(rep(1,4),4))

CRMSE<-matrix(x$RMSE[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CRMSE<-x$RMSE-kronecker(CRMSE,matrix(rep(1,4),4))

#Line plots to compare estimators
qplot(as.factor(n),CRMSE,data=x,geom="path",colour=Estimator,size=I(1.5),group=Estimator,facets=nu~Dist,xlab="Sample Size")
qplot(as.factor(n),CMedian,data=x,geom="path",colour=Estimator,size=I(1.5),group=Estimator,facets=nu~Dist,xlab="Sample Size")
qplot(as.factor(n),CBias,data=x,geom="path",colour=Estimator,size=I(1.5),group=Estimator,facets=nu~Dist,xlab="Sample Size")

#Compare the Manton L2 and Projected Arithemetic 
qplot(ML2Error,ArithError,data=Res[Res$n==300,],xlab=expression(d[R](hat(bold(S))[L2],bold(S))),ylab=expression(d[R](hat(bold(S))[P],bold(S))))+
geom_abline(intercept=0,slope=1)+facet_grid(nu~Dist)

qplot(HL1Error,MedError,data=Res[Res$n==300,],xlab=expression(d[R](hat(bold(S))[L1],bold(S))),ylab=expression(d[R](hat(bold(S))[M],bold(S))))+
geom_abline(intercept=0,slope=1)+facet_grid(nu~Dist)


################################################################
###########Recreate plots in paper with right results###########
################################################################
################################################################

#Box plots for high variablility and n=300#
Largennu<-ResFrame[ResFrame$nu==.75,]
Largennu<-Largennu[Largennu$n==300,]
max<-max(Largennu$Error)
levels(Largennu$Estimator)<-c("SR1","SR2","SE2","SE1")

setwd("U://Thesis//PointEstimationPaper//images")
#pdf("CayleyBoxN300Nu75.pdf")
qplot(Estimator,Error,geom="boxplot",data=Largennu[Largennu$Dist=="Cayley",],ylim=c(0,max))
#dev.off()

#pdf("MisesBoxN300Nu75.pdf")
qplot(Estimator,Error,geom="boxplot",data=Largennu[Largennu$Dist=="Mises",],ylim=c(0,max))
#dev.off()

#pdf("FisherBoxN300Nu75.pdf")
qplot(Estimator,Error,geom="boxplot",data=Largennu[Largennu$Dist=="Fisher",],ylim=c(0,max))
#dev.off()

Centered (around S_P) Median, Bias and RMSE for our four estimators
CMedian<-matrix(x$Median[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CMedian<-x$Median-kronecker(CMedian,matrix(rep(1,4),4))

CBias<-matrix(x$Bias[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CBias<-x$Bias-kronecker(CBias,matrix(rep(1,4),4))

CRMSE<-matrix(x$RMSE[seq(3,nrow(x),by=4)],nrow(x)/4)
x$CRMSE<-x$RMSE-kronecker(CRMSE,matrix(rep(1,4),4))

levels(x$Dist)[c(2,3)]<-c("von Mises-Fisher","von Mises Circular")

x2<-x[x$Estimator!="Projected Mean",]
x2$Est<-factor(x2$Estimator[x2$Estimator!="Projected Mean"])
levels(x2$Est)<-c("SL1","SL2","SM")

#pdf("centeredMed.pdf",width=8,height=5)
qplot(as.factor(n),CMedian,geom="path",group=Est,colour=Est,data=x2,facets=nu~Dist,
	size=I(1.4),xlab="Sample Size",ylab="Centered Median")+geom_hline(yintercept=0)
#dev.off()

#pdf("centeredRMSE.pdf",width=8,height=5)
qplot(as.factor(n),CRMSE,geom="path",group=Est,lty=Est,data=x2,facets=nu~Dist,
	size=I(1.4),xlab="Sample Size",ylab="Centered RMSE")+geom_hline(yintercept=0)
#dev.off()

#pdf("centeredMean.pdf",width=8,height=5)
qplot(as.factor(n),CBias,geom="path",group=Est,lty=Est,data=x2,facets=nu~Dist,
	size=I(1.4),xlab="Sample Size",ylab="Centered Bias")+geom_hline(yintercept=0)
#dev.off()

#Tabular form of those plots when nu=0.75

x75<-x[x$nu==.75,]

xtable(x75)

####
##Direct comparisons of Euclid and Riemann for n=50 and the two extreme nu values
####

cResFrame<-cast(ResFrame,n+nu+Sample+Dist~Estimator)
names(cResFrame)[5:8]<-c("HL1","ML2","PM","Med")
Midn<-cResFrame[cResFrame$n==100,]

levels(Midn$Dist)[2:3]<-c("von Mises-Fisher","von Mises Circular")
xmax25<-max(Midn[Midn$nu==.25,5:8])
xmax75<-max(Midn[Midn$nu==.75,8])

#pdf("SMvsSL1Nu25.pdf",width=8,height=5)
qplot(Med,HL1,data=Midn[Midn$nu==.25,],facets=.~Dist,xlab=expression(hat(S)[M]),ylab=expression(hat(S)[L[1]]),
	ylim=c(0,xmax25),xlim=c(0,xmax25))+geom_abline(intercept=0,slope=1)
#dev.off()

#pdf("SMvsSL1Nu75.pdf",width=8,height=5)
qplot(Med,HL1,data=Midn[Midn$nu==.75,],facets=.~Dist,xlab=expression(hat(S)[M]),ylab=expression(hat(S)[L[1]]),
	ylim=c(0,.7),xlim=c(0,.7))+geom_abline(intercept=0,slope=1)
#dev.off()

#pdf("SPvsSL2Nu25.pdf",width=8,height=5)
qplot(PM,ML2,data=Midn[Midn$nu==.25,],facets=.~Dist,xlab=expression(hat(S)[P]),ylab=expression(hat(S)[L[2]]),
	ylim=c(0,xmax25),xlim=c(0,xmax25))+geom_abline(intercept=0,slope=1)
#dev.off()

#pdf("SPvsSL2Nu75.pdf",width=8,height=5)
qplot(PM,ML2,data=Midn[Midn$nu==.75,],facets=.~Dist,xlab=expression(hat(S)[P]),ylab=expression(hat(S)[L[2]]),
	ylim=c(0,.7),xlim=c(0,.7))+geom_abline(intercept=0,slope=1)
#dev.off()

#Tables that compute % above/below line and distance between
CaySumL1<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(Med-HL1),perc=sum(HL1<Med)/1000)
FisSumL1<-ddply(cResFrame[cResFrame$Dist=="Fisher",],.(nu,n),summarize,rbar=mean(Med-HL1),perc=sum(HL1<Med)/1000)
MisSumL1<-ddply(cResFrame[cResFrame$Dist=="Mises",],.(nu,n),summarize,rbar=mean(Med-HL1),perc=sum(HL1<Med)/1000)
SumL1<-cbind(CaySumL1,FisSumL1[,3:4],MisSumL1[,3:4])
xtable(SumL1,digits=3,label="tab:percL1")

CaySumL2<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(PM-ML2),perc=sum(ML2<PM)/1000)
FisSumL2<-ddply(cResFrame[cResFrame$Dist=="Fisher",],.(nu,n),summarize,rbar=mean(PM-ML2),perc=sum(ML2<PM)/1000)
MisSumL2<-ddply(cResFrame[cResFrame$Dist=="Mises",],.(nu,n),summarize,rbar=mean(PM-ML2),perc=sum(ML2<PM)/1000)
SumL2<-cbind(CaySumL2,FisSumL2[,3:4],MisSumL2[,3:4])
xtable(SumL2,digits=3,label="tab:percL2")

####
#Make tables to compare estimators for von Mises-Circular and nu=.75
####

vMRes<-x[x$nu==0.75 & x$Dist=="von Mises Circular",][,-c(1:2)]
vMTab<-cbind(vMRes[1:8,1:5],vMRes[9:16,c(1,2:5)])
xtable(vMTab,digits=3)




