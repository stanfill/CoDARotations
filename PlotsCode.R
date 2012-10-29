####################################################################
###########			All the density functions	##############
####################################################################
setwd("C:/Users/stanfill/My Dropbox/Rotation matrices/Paper/PointEstimationPaper/submissionv2/Supporting Materials/images")

library(ggplot2)
library(reshape2)
cayden<-function(r,kappa,Haar=T){
	den<-.5*gamma(kappa+2)/(sqrt(pi)*2^kappa*gamma(kappa+.5))*(1+cos(r))^kappa*(1-cos(r))

	if(Haar)
		return(den/(1-cos(r)))
	else
		return(den)
}
fishden<-function(r,kappa,Haar=T){
	den<-exp(2*kappa*cos(r))*(1-cos(r))/(2*pi*(besselI(2*kappa,0)-besselI(2*kappa,1)))

	if(Haar)
		return(den/(1-cos(r)))

	else
		return(den)
}

vonmden<-function(r,kappa,Haar=T){
	den<-1/(2*pi*besselI(kappa,0))*exp(kappa*cos(r))
	
	if(Haar)
		return(den/(1-cos(r)))

	else
		return(den)
}

r<-seq(-pi,pi,length=1000)

cayKap<-c(10,4,2,.065)
fishKap<-c(3.165,1.711,1.1564,.175)
misesKap<-c(2.40, 1.159, 0.5164)

denData<-data.frame(r=r,Fisher=fishden(r,fishKap[1],T),Mises=vonmden(r,misesKap[1],T),Cayley=cayden(r,cayKap[1],T),nu=0.25)
denData<-rbind(denData,data.frame(r=r,Fisher=fishden(r,fishKap[2],T),Mises=vonmden(r,misesKap[2],T),Cayley=cayden(r,cayKap[2],T),nu=0.5))
denData<-rbind(denData,data.frame(r=r,Fisher=fishden(r,fishKap[3],T),Mises=vonmden(r,misesKap[3],T),Cayley=cayden(r,cayKap[3],T),nu=0.75))

denData<-melt(denData,id=c("r","nu"))
colnames(denData)<-c("r","Variance","Density","Value")

levels(denData$Density)<-c("matrix Fisher","circ-von Mises","Cayley")
denData$Density<-factor(denData$Density,levels=levels(denData$Density)[c(3,1,2)])

denData2<-data.frame(r=r,Fisher=fishden(r,fishKap[1],F),Mises=vonmden(r,misesKap[1],F),Cayley=cayden(r,cayKap[1],F),nu=0.25)
denData2<-rbind(denData2,data.frame(r=r,Fisher=fishden(r,fishKap[2],F),Mises=vonmden(r,misesKap[2],F),Cayley=cayden(r,cayKap[2],F),nu=0.5))
denData2<-rbind(denData2,data.frame(r=r,Fisher=fishden(r,fishKap[3],F),Mises=vonmden(r,misesKap[3],F),Cayley=cayden(r,cayKap[3],F),nu=0.75))

denData2<-melt(denData2,id=c("r","nu"))
colnames(denData2)<-c("r","Variance","Density","Value")

levels(denData2$Density)<-c("matrix Fisher","circ-von Mises","Cayley")
denData2$Density<-factor(denData2$Density,levels=levels(denData2$Density)[c(3,1,2)])

#setwd("U:/Thesis/PointEstimationPaper/images")

##################################################
## Whole density pictures
##################################################

#With respect to Haar
theme_set(theme_black(16))
#pdf("Var25DensityHaar.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.25,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,12),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(0,1),legend.justification=c(0.1,1-.11),legend.background=element_rect(fill="white",linetype=0))+
  theme(legend.title=element_text(size=15,face="bold"),legend.text=element_text(size=14,face="bold"))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()

#pdf("Var5DensityHaar.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.5,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,4.5),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position="none")+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()

#pdf("Var75DensityHaar.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,2.25),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position="none")+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()


#Associate editor wants to see who distribution, do this for nu=0.75:
setwd("\\\\iastate.edu/cyfiles/stanfill/Desktop/GitHub/CoDARotations/images")

qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,25),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(0,1),legend.justification=c(0.1,1-.11),legend.background=element_rect(fill="white",linetype=0))+
  theme(legend.title=element_text(size=15,face="bold"),legend.text=element_text(size=14,face="bold"))+geom_abline(intercept=0,slope=0,colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))+geom_vline(xintercept=0)
ggsave("Var75DensityHaarFull.pdf",height=5,width=5)


qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(0,1),legend.justification=c(0.1,1-.11),legend.background=element_rect(fill="white",linetype=0))+
  theme(legend.title=element_text(size=15,face="bold"),legend.text=element_text(size=14,face="bold"))+geom_abline(intercept=0,slope=0,colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))+geom_vline(xintercept=0)
ggsave("Var75DensityHaarFullFull.pdf",height=5,width=5)

#With respect to lebesgue

#pdf("Var25Density.pdf",height=5,width=5)
qplot(r,Value,data=denData2[denData2$Variance==.25,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.75),ylab="f(r)")+
 scale_linetype_manual(values=c(1,12,3))+opts(legend.position=c(0,1),legend.justification=c(0.1,1-.11),legend.background=theme_rect(fill="white",linetype=0))+
 opts(legend.title=theme_text(size=10,face="bold"),legend.text=theme_text(size=10,face="bold"))+geom_abline(intercept=0,slope=c(0,100000))+
 opts(axis.text.x=theme_text(size=12))+opts(axis.text.y=theme_text(size=12))
#dev.off()

#pdf("Var5Density.pdf",height=5,width=5)
qplot(r,Value,data=denData2[denData2$Variance==.5,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.5),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+opts(legend.position="none")+geom_abline(intercept=0,slope=c(0,100000))+
  opts(axis.text.x=theme_text(size=12))+opts(axis.text.y=theme_text(size=12))
#dev.off()

#pdf("Var75Density.pdf",height=5,width=5)
qplot(r,Value,data=denData2[denData2$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.4),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+opts(legend.position="none")+geom_abline(intercept=0,slope=c(0,100000))+
  opts(axis.text.x=theme_text(size=12))+opts(axis.text.y=theme_text(size=12))
#dev.off()


##################################################
## Densities zoomed in
##################################################

Boxx<-c(1,1,pi,pi,1)
Boxy<-c(0,.8,.8,0,0)

#pdf("Var75DensityBox.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.75,],geom="blank",ylim=c(0,2.25),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position="none")+geom_path(aes(x=Boxx,y=Boxy),lwd=I(1.25),colour="gray")+
  geom_line(aes(x=r,y=Value,linetype=Density),lwd=I(1.25))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()


#pdf("Var75DensityZoom.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.8),xlim=c(1,pi),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(fill="white",linetype=0))+
  theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=16,face="bold"))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()

#pdf("Var75DensityZoomNoGuide.pdf",height=5,width=5)
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.8),xlim=c(1,pi),ylab="f(r)")+
  scale_linetype_manual(values=c(1,12,3))+theme(legend.position="NONE")+
  theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=16,face="bold"))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
  theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))
#dev.off()


###############################################################################
######                                                              ###########
######      Analyze simulation results                              ###########
######                                                              ###########
###############################################################################
library(ggplot2)
library(reshape2)
library(plyr)
library(grid)

#theme_set(theme_black(16))

setwd("U://Thesis//Simulation Code//Results")
Res<-read.csv("FullFinalResults.csv")[,-c(1)]
setwd("\\\\iastate.edu/cyfiles/stanfill/Desktop/GitHub/CoDARotations/images")

ResFrame<-melt(Res,id=c("nu","n","Dist","Sample"),measure.var=c("HL1Error","ML2Error","ArithError","MedError"))
colnames(ResFrame)[5:6]<-c("Estimator","Error")
levels(ResFrame$Estimator)<-c("R.Median","R.Mean","E.Mean","E.Median")
levels(ResFrame$Dist)<-c("Cayley","matrix Fisher","circular-von Mises")

x<-ddply(ResFrame,.(Dist,nu,n,Estimator),summarize,Median=round(median(Error),4),Mean=round(mean(Error),4),RMSE=round(sqrt(mean(Error^2)),4))

#Make previous Table 4 into a plot for Associate editor
my.labels <- list(bquote(widehat(S)[E]),bquote(widehat(S)[R]),bquote(widetilde(S)[E]),bquote(widetilde(S)[R]))

mx<-melt(x,id=c("Dist","Estimator","n","nu"),measure=c("Mean","RMSE"))
mx$n<-as.factor(mx$n)
mx75<-mx[(mx$nu==0.75&mx$Dist=="circular-von Mises"),]
mx75$Estimator<-factor(mx75$Estimator,levels=c("E.Mean","R.Mean","E.Median","R.Median"))

qplot(n,value,data=mx75,facets=.~variable,geom="path",group=Estimator,linetype=Estimator,lwd=I(1))+
  scale_linetype_manual(values=c(3,4,2,1),labels=my.labels)+coord_equal(ratio=4)+
  theme(legend.text=element_text(size=12),legend.key.width=unit(3,"line"),legend.title=element_text(size=12))+
  geom_hline(yintercept=0,colour="gray50")+
  theme(axis.text.x=element_text(size=12,color=1),axis.text.y=element_text(size=12,color=1))
ggsave("vonMisesnu75MeanRMSE.pdf",height=3,width=6)

#Plot boxplots as a function of nu for n=300
Largen<-ResFrame[ResFrame$n==100,]
Largen<-Largen[-which.max(Largen$Error),]
Largen<-Largen[-which.max(Largen$Error),]

Largen$nu<-as.factor(Largen$nu)
Largen$nu<-factor(Largen$nu,labels=c("nu == 0.25","nu == 0.50","nu == 0.75"))
levels(Largen$Dist)<-c("Cayley","matrix~~Fisher","circular-von~~Mises")

#THIS IS FOR github VERSION OF ggplot2
#setwd("U://Thesis//PointEstimationPaper//images")

qplot(Estimator,Error,geom="boxplot",data=Largen,xlab="",ylab=expression(d[R](bold(S),.)))+
  geom_hline(xintercept=0,colour="gray50")+
  facet_grid(nu~Dist,scales="free",labeller=label_parsed)+
  scale_x_discrete(limits=c("E.Mean","R.Mean","E.Median","R.Median"),breaks=c("E.Mean","R.Mean","E.Median","R.Median"),labels=c(expression(widehat(bold(S))[E]),expression(widehat(bold(S))[R]),expression(widetilde(bold(S))[E]),expression(widetilde(bold(S))[R])))+
  theme(axis.text.x=element_text(size=12,color=1,face='bold'),axis.text.y=element_text(size=12,color=1))
#ggsave("N100AllNuBoxes.pdf",width=9)

#Plot nu=.75 as a function of n.  Remove the bad observations that make the plot useless for n=50,100
Largenu<-ResFrame[ResFrame$nu==.75,]

Largenu50<-Largenu[Largenu$n==50,]
badeggs<-Largenu50[Largenu50$Error>2,]$Sample

Largenu100<-Largenu[Largenu$n==100,]
badeggs<-c(badeggs,Largenu100[Largenu100$Error>2,]$Sample)

beggs<-which(Largenu$Sample %in% badeggs)
Largenu<-Largenu[-beggs,]

Largenu$n<-as.factor(Largenu$n)
levels(Largenu$n)<-c("n = 10","n = 50","n = 100","n = 300")

qplot(Estimator,Error,geom="boxplot",data=Largenu,xlab="",ylab=expression(d[R](bold(S),.)))+
  facet_grid(n~Dist,scales="free")+
  geom_hline(xintercept=0,colour="gray50")+
  scale_x_discrete(limits=c("E.Mean","R.Mean","E.Median","R.Median"),breaks=c("E.Mean","R.Mean","E.Median","R.Median"),labels=c(expression(widehat(bold(S))[E]),expression(widehat(bold(S))[R]),expression(widetilde(bold(S))[E]),expression(widetilde(bold(S))[R])))+
  theme(axis.text.x=element_text(size=12,color=1,face='bold'),axis.text.y=element_text(size=12,color=1))
#ggsave("Nu75AllNBoxes.pdf",width=9)

####
##Direct comparisons of Euclid and Riemann for n=100 and the two extreme nu values
####
library(reshape2)
library(plyr)
cResFrame<-dcast(ResFrame,n+nu+Sample+Dist~Estimator,value.var="Error")
Midn<-cResFrame[cResFrame$n==100,]

#levels(Midn$Dist)[2:3]<-c("matrix Fisher","circular-von Mises")
xmax25<-max(Midn[Midn$nu==.25,5:8])
xmax75<-max(Midn[Midn$nu==.75,c(5,7:8)])

#setwd("U:/Thesis/PointEstimationPaper/images")

ggplot(Midn[Midn$nu==.25,],aes(E.Median,R.Median))+facet_grid(.~Dist)+geom_point()+geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray50")+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(d[R](bold(S),widetilde(bold(S))[E])),limits=c(0,xmax25))+
  scale_y_continuous(expression(d[R](bold(S),widetilde(bold(S))[R])),limits=c(0,xmax25))
#ggsave("SMvsSL1Nu25.pdf",width=8,height=4)

ggplot(Midn[Midn$nu==.75,],aes(E.Median,R.Median))+facet_grid(.~Dist)+geom_point()+geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray50")+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(d[R](bold(S),widetilde(bold(S))[E])),limits=c(0,xmax75))+
  scale_y_continuous(expression(d[R](bold(S),widetilde(bold(S))[R])),limits=c(0,xmax75))
#ggsave("SMvsSL1Nu75.pdf",width=8,height=4)


ggplot(Midn[Midn$nu==.25,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point()+geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray50")+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(d[R](bold(S),widehat(bold(S))[E])),limits=c(0,xmax25))+
  scale_y_continuous(expression(d[R](bold(S),widehat(bold(S))[R])),limits=c(0,xmax25))
#ggsave("SPvsSL2Nu25.pdf",width=8,height=4)


ggplot(Midn[Midn$nu==.75,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point()+geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray50")+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(d[R](bold(S),widehat(bold(S))[E])),limits=c(0,xmax75))+
  scale_y_continuous(expression(d[R](bold(S),widehat(bold(S))[R])),limits=c(0,xmax75))
#ggsave("SPvsSL2Nu75.pdf",width=8,height=4)

####
#Combine figures 5 and 6 ala reviewer 2

ggplot(Midn[Midn$nu==.25,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point(alpha=I(.75))+geom_point(aes(E.Median,R.Median),color="grey50",alpha=I(.75))+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(paste(d[E],"-based estimators")),limits=c(0,xmax25))+
  scale_y_continuous(expression(paste(d[R],"-based estimators")),limits=c(0,xmax25))+
  geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray70")+ggtitle(expression(nu==0.25))
#ggsave("EuclidRiemannNu25.pdf",width=9,height=4)

ggplot(Midn[Midn$nu==.75,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point(alpha=I(.75))+geom_point(aes(E.Median,R.Median),color="grey50",alpha=I(.75))+
  coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
  scale_x_continuous(expression(paste(d[E],"-based estimators")),limits=c(0,xmax75))+
  scale_y_continuous(expression(paste(d[R],"-based estimators")),limits=c(0,xmax75))+
  geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray70")+ggtitle(expression(nu==0.75))
#ggsave("EuclidRiemannNu75.pdf",width=9,height=4)

########
#Tables that compute % above/below line and distance between
#ORIGINAL VERSION WITH NO SEs
##

CaySumL1<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),perc=sum(R.Median<E.Median)/1000)
FisSumL1<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),perc=sum(R.Median<E.Median)/1000)
MisSumL1<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),perc=sum(R.Median<E.Median)/1000)
SumL1<-cbind(CaySumL1,FisSumL1[,3:4],MisSumL1[,3:4])
xtable(SumL1,digits=4,label="tab:percL1")

CaySumL2<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),perc=sum(R.Mean<E.Mean)/1000)
FisSumL2<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),perc=sum(R.Mean<E.Mean)/1000)
MisSumL2<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),perc=sum(R.Mean<E.Mean)/1000)
SumL2<-cbind(CaySumL2,FisSumL2[,3:4],MisSumL2[,3:4])
xtable(SumL2,digits=3,label="tab:percL2")

##########################################################################
##########################################################################
####   This section is dedicated to proving statistical significance  ####
####     between estimators                                           ####
##########################################################################
##########################################################################
####
#Non-parametric ANOVA Table to show differences in estimators
####
#i indicates nu, j indicates n and k indicates dist
nui<-unique(ResFrame$nu)
nj<-unique(ResFrame$n)
distk<-unique(ResFrame$Dist)
kres<-vector("list",length(nui)*length(nj)*length(distk))
ind<-1

for(i in nui){
  for(j in nj){
    for(k in distk){
      subResFrame<-subset(ResFrame,n==j & nu==i & Dist==k)
      kres[[ind]]<-kruskal.test(Error~Estimator,data=subResFrame)
      ind<-ind+1
    }
  }
}

kres

#Remove effect of n and nu then do one rank sum test for estimator and distribution

sumResFrame<-ddply(ResFrame,.(nu,n),summarize,yddd=mean(Error))
mergeResFrame<-merge(ResFrame,sumResFrame)
mergeResFrame$Error2<-mergeResFrame$Error-mergeResFrame$yddd

qplot(Estimator,Error2,facets=.~Dist,data=mergeResFrame,geom='boxplot')
kruskal.test(Error2~Dist+Estimator,data=mergeResFrame)

#Perhaps just get SEs for table in paper 
n100nu25ResFrame<-subset(ResFrame,n==100 & nu==0.25)
n100nu25ses<-ddply(n100nu25ResFrame,.(Dist,Estimator),summarize,mean=mean(Error),sd=sd(Error)/sqrt(length(Error)))
n100nu25ses

#Add ses to delta table
CaySumL1<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
FisSumL1<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
MisSumL1<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
SumL1<-cbind(CaySumL1[,2:5],FisSumL1[,3:5],MisSumL1[,3:5])
xtable(SumL1,digits=4,label="tab:percL1")

CaySumL2<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
FisSumL2<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
MisSumL2<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
SumL2<-cbind(CaySumL2[,2:5],FisSumL2[,3:5],MisSumL2[,3:5])
xtable(SumL2,digits=4,label="tab:percL2")

#Perform ANOVA on log transformed data
#Make one table treating nu/n/Dist combination as blocks
Error2<-log(ResFrame$Error)
qqnorm(Error2);qqline(Error2)

shapiro.test(Error2[5001:10000])
aovResFrame<-ResFrame
aovResFrame$logError<-log(aovResFrame$Error)
aovResFrame$Blocks<-as.factor(paste(aovResFrame$nu,aovResFrame$n,aovResFrame$Dist))

a1<-aov(logError~n*nu*Dist+Estimator,data=aovResFrame)
summary(a1)

a2<-aov(logError~Blocks+Estimator,data=aovResFrame)
summary(a2)
xtable(summary(a2))

#Do one ANOVA just for n=100, nu=0.25 to match table in Appendix
aovResFramen100nu25<-subset(aovResFrame,n==100&nu==0.25)
a3<-aov(logError~Estimator,data=aovResFramen100nu25)
summary(a3)
xtable(summary(a3))
##########################################################################
##########################################################################
##########################################################################
##########################################################################

####
#Make tables to compare estimators for von Mises-Circular and nu=.75
####

vMRes<-x[x$nu==0.75 & x$Dist=="von Mises Circular",][,-c(1:2)]
vMTab<-cbind(vMRes[1:8,1:5],vMRes[9:16,c(1,2:5)])
xtable(vMTab,digits=3)

####
#Make tables to compare estimators for n=100 across distributions
####

n100Res<-x[x$n==100 & x$nu==.25,][,-c(2:3)]
#n100Tab<-cbind(n100Res[1:8,1:5],n100Res[9:16,c(1,2:5)])
xtable(n100Res,digits=3)


###############################################################
####### Plot tail weight versus                 ###############
####### Estimator error                         ###############
###############################################################

setwd("\\\\iastate.edu/cyfiles/stanfill/Desktop/GitHub/CoDARotations")
library(rotations)
library(splines)

alldf<-read.csv("Nu75TailBehavior.csv") 
alldf$Dist<-factor(alldf$Dist,levels=levels(alldf$Dist)[c(1,3,2)])
tail75<-mean(c(2.109,2.185,1.975))
realProbs<-rep(0,3)
realProbs[1]<-integrate(dcayley,lower=tail75,upper=pi,Haar=F,kappa=2)$value*2
realProbs[2]<-integrate(dfisher,lower=tail75,upper=pi,Haar=F,kappa=1.15)$value*2
realProbs[3]<-integrate(dvmises,lower=tail75,upper=pi,Haar=F,kappa=0.52)$value*2

ggplot(alldf[alldf$n>99,],aes(Prop,Pdiff))+xlab("Proportion Observations in Tail")+ylab(expression(d[G](widehat(bold(S))[E],bold(S))-d[G](widetilde(bold(S))[E],bold(S))))+
  geom_hline(yintercept=0,colour="gray50")+
  geom_vline(xintercept=realProbs[1],colour=2,alpha=I(.5))+
  geom_vline(xintercept=realProbs[2],colour=3,alpha=I(.5))+
  geom_vline(xintercept=realProbs[3],colour=4,alpha=I(.5))+
  stat_smooth(method=lm,formula=y~ns(x,2),fullrange=T,colour=1)+
  geom_point(aes(shape=Dist),alpha=I(.6))+theme_bw()+scale_shape_manual(values=c(16,0,17))+
  facet_grid(n~.,scales="free_y",labeller=label_bquote(n==.(x)))

setwd("\\\\iastate.edu/cyfiles/stanfill/Desktop/GitHub/CoDARotations/images")

ggplot(alldf[alldf$n>299,],aes(Prop,Pdiff))+xlab("Proportion Observations in Tail")+
  ylab(expression(d[R](bold(S),widehat(bold(S))[E])-d[R](bold(S),widetilde(bold(S))[E])))+
  geom_hline(yintercept=0,colour="gray50")+
  stat_smooth(method=lm,formula=y~ns(x,2),fullrange=T,colour=1)+
  geom_point(aes(shape=Dist),alpha=I(.6))+theme_bw()+
  scale_shape_manual(values=c(16,0,17),name="Distribution")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))
#ggsave("Nu75N300TailBehavior.pdf",width=8,height=4)

alldf$ScalePdiff<-alldf$Pdiff/(alldf$PMean+alldf$PMedian)

ggplot(alldf[alldf$n>299,],aes(Prop,ScalePdiff))+xlab("Proportion Observations in Tail")+
  ylab(expression(frac(d[R](bold(S),widehat(bold(S))[E])-d[R](bold(S),widetilde(bold(S))[E]),d[R](bold(S),widehat(bold(S))[E])+d[R](bold(S),widetilde(bold(S))[E]))))+
  geom_hline(yintercept=0,colour="gray50")+
  stat_smooth(method=lm,formula=y~ns(x,2),fullrange=T,colour=1)+
  geom_point(aes(shape=Dist),alpha=I(.6))+theme_bw()+
  scale_shape_manual(values=c(16,0,17),name="Distribution")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(legend.justification=c(0,1),legend.position=c(0,1))
#ggsave("Nu75N300TailBehaviorStandard.pdf",width=7,height=4)


###############################################################
####### Find the cases where the                ###############
####### Riemannian estimate really sucks        ###############
###############################################################

ExtDat<-Midn[Midn$nu==.75,]
qplot(E.Mean,R.Mean,data=ExtDat,facets=.~Dist,xlab=expression(hat(S)[P]),ylab=expression(hat(S)[L[2]]))+geom_abline(intercept=0,slope=1)

ext<-which.max(ExtDat$R.Mean)
ExtSampMax<-ExtDat[ext,]$Sample
ExtDatFish<-ExtDat[ExtDat$Dist=="von Mises-Fisher",]
SampMax<-which(ExtDatFish$Sample==ExtSampMax)

setwd("U:/Thesis/R Package")
source("RoxygenFun.R")
setwd("U:/Thesis/Simulation Code/DataFiles")
fisherSamples<-read.csv("AltFisherVar-75.csv")[,-1]

rowMax<-10*1000+50*1000+SampMax*100

ExtSampleMax<-as.matrix(fisherSamples[(rowMax-100+1):(rowMax),],100,9)
SPMax<-riedist(arith.mean(ExtSampleMax))
Confirm<-Res[Res$Sample==ExtSampMax,]
#This will be very near zero if sample is correct
SPMax-Confirm$ArithError

riedist(MantonL2(ExtSampleMax)$R)

#Given the sample, whats the sampling dist look like?
rMax<-rep(0,100)
for(i in 1:100){
  rMax[i]<-riedist(ExtSampleMax[i,])
}
hist(rMax,breaks=100)

#Make a matrix of all the pairwise distances
MaxDistMat<-matrix(0,100,100)
i=j=1

for(i in 1:100){
  j=1
  while(j<i){
    MaxDistMat[i,j]<-riedist(matrix(ExtSampleMax[i,],3,3),matrix(ExtSampleMax[j,],3,3))
    j=j+1
  }
}

max(MaxDistMat)

#Rearrange the rows
ExtSampleMax<-ExtSampleMax[c(4,1:3,5:100),]
riedist(MantonL2(ExtSampleMax)$R)

#Now find sample where Rieman works well
MinExtSample<-ExtDatFish[which.min(ExtDatFish$R.Mean-ExtDatFish$E.Mean),]$Sample
SampMin<-which(ExtDatFish$Sample==MinExtSample)


rowMin<-10*1000+50*1000+SampMin*100

ExtSampleMin<-as.matrix(fisherSamples[(rowMin-100+1):(rowMin),],100,9)
SPMin<-riedist(arith.mean(ExtSampleMin))
Confirm<-Res[Res$Sample==MinExtSample,]
#This will be very near zero if sample is correct
SPMin-Confirm$ArithError

riedist(MantonL2(ExtSampleMin)$R)

#Given the sample, whats the sampling dist look like?
rMin<-rep(0,100)
for(i in 1:100){
  rMin[i]<-riedist(ExtSampleMin[i,])
}


rs<-data.frame(r=c(rMin,rMax),Type=c(rep("Min",100),rep("Max",100)))
qplot(r,data=rs,colour=Type,group=Type,geom="density")

#Make a matrix of all the pairwise distances
MinDistMat<-matrix(0,100,100)
i=j=1

for(i in 1:100){
  j=1
  while(j<i){
    MinDistMat[i,j]<-riedist(matrix(ExtSampleMin[i,],3,3),matrix(ExtSampleMin[j,],3,3))
    j=j+1
  }
}

max(MinDistMat)





