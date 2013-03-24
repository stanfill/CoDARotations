#Load the necessary libraries and create the 'tests' function used to create 
#Table 2 of the Supplementary Materials
library(rotations)
library(reshape2)
library(plyr)
library(splines)
library(xtable)

tests<-function(x){
	num<-ncol(x)
	pval<-rep(0,num)
	for(i in 1:num){
		pval[i]<-t.test(x[,i])$p.value
	} 
	return(pval)
}

######################
##In this section the plots in the Section 4 (Simulation Study) 
## and some of the Supplementary Materials plots are made
######################

r<-seq(-pi,pi,length=1000)

cayKap<-c(10,4,2,.065)
fishKap<-c(3.165,1.711,1.1564,.175)
misesKap<-c(2.40, 1.159, 0.5164)

denData<-data.frame(r=r,Fisher=dfisher(r,fishKap[1],T),Mises=dvmises(r,misesKap[1],T),Cayley=dcayley(r,cayKap[1],T),nu=0.25)
denData<-rbind(denData,data.frame(r=r,Fisher=dfisher(r,fishKap[2],T),Mises=dvmises(r,misesKap[2],T),Cayley=dcayley(r,cayKap[2],T),nu=0.5))
denData<-rbind(denData,data.frame(r=r,Fisher=dfisher(r,fishKap[3],T),Mises=dvmises(r,misesKap[3],T),Cayley=dcayley(r,cayKap[3],T),nu=0.75))

denData<-melt(denData,id=c("r","nu"))
colnames(denData)<-c("r","Variance","Density","Value")

levels(denData$Density)<-c("matrix Fisher","circ-von Mises","Cayley")
denData$Density<-factor(denData$Density,levels=levels(denData$Density)[c(3,1,2)])

#This is Figure 1(a) of the manuscript
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,25),ylab="f(r)")+geom_abline(intercept=0,slope=0,colour="gray50")+
	theme(legend.position='NONE')+scale_linetype_manual(values=c(1,12,3))+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))+geom_vline(xintercept=0)


#This is Figure 1(b) of the manuscript
Boxx<-c(1.75,1.75,2.5,2.5,1.75)
Boxy<-c(0,.25,.25,0,0)
qplot(r,Value,data=denData[denData$Variance==.75,],geom="blank",ylim=c(0,2.25),ylab="f(r)")+
	scale_linetype_manual(values=c(1,12,3))+theme(legend.position="none")+geom_path(aes(x=Boxx,y=Boxy),lwd=I(1.25),colour="gray")+
	geom_line(aes(x=r,y=Value,linetype=Density),lwd=I(1.25))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))

#This is Figure 1(c) of the manuscript
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.25),xlim=c(1.75,2.5),ylab="f(r)")+
	scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(fill="white",linetype=0))+
	theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=16,face="bold"),legend.key.width=unit(3,"line"))+
	geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))


#This will produce a plot similar to Figure 3(a) and Figures 2(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rcayley,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
ggsave("eye-cayley.pdf",  height=5, width=5)
plot(Rs,center=id.SO3, col=2) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=3) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")






#This will produce a plot similar to Figure 3(b) and Figures 3(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rfisher,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
ggsave("eye-fisher.pdf",  height=5, width=5)
plot(Rs,center=id.SO3, col=2) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=3) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")



#This will produce a plot similar to Figure 3(c) and Figures 4(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rvmises,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
ggsave("eye-vmises.pdf",  height=5, width=5)
plot(Rs,center=id.SO3, col=2) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=3) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")

######################
##In this section the plots in the Section 5 (Simulation Results)
######################

#Read in the data and manipulate it to use ggplot2 easily
Res<-read.csv("FullFinalResults.csv")[,-1]
alldf<-read.csv("Nu75TailBehavior.csv")[,-1]

ResFrame<-melt(Res,id=c("nu","n","Dist","Sample"),measure.var=c("HL1Error","ML2Error","ArithError","MedError"))
colnames(ResFrame)[5:6]<-c("Estimator","Error")
levels(ResFrame$Estimator)<-c("R.Median","R.Mean","E.Mean","E.Median")
levels(ResFrame$Dist)<-c("Cayley","matrix Fisher","circular-von Mises")

Largen<-ResFrame[ResFrame$n==100,]
Largen<-Largen[-which.max(Largen$Error),]
Largen<-Largen[-which.max(Largen$Error),]

Largen$nu<-as.factor(Largen$nu)
Largen$nu<-factor(Largen$nu,labels=c("nu == 0.25","nu == 0.50","nu == 0.75"))
levels(Largen$Dist)<-c("Cayley","matrix~~Fisher","circular-von~~Mises")

x<-ddply(ResFrame,.(Dist,nu,n,Estimator),summarize,Median=round(median(Error),4),Mean=round(mean(Error),4),RMSE=round(sqrt(mean(Error^2)),4))
my.labels <- list(bquote(widehat(bold(S))[E]),bquote(widehat(bold(S))[R]),bquote(widetilde(bold(S))[E]),bquote(widetilde(bold(S))[R]))
mx<-melt(x,id=c("Dist","Estimator","n","nu"),measure=c("Mean","RMSE"))
mx$n<-as.factor(mx$n)
mx75<-mx[(mx$nu==0.75&mx$Dist=="circular-von Mises"),]
mx75$Estimator<-factor(mx75$Estimator,levels=c("E.Mean","R.Mean","E.Median","R.Median"))


alldf$Dist<-factor(alldf$Dist,levels=levels(alldf$Dist)[c(1,3,2)])
tail75<-mean(c(2.109,2.185,1.975))
realProbs<-rep(0,3)
realProbs[1]<-integrate(dcayley,lower=tail75,upper=pi,Haar=F,kappa=2)$value*2
realProbs[2]<-integrate(dfisher,lower=tail75,upper=pi,Haar=F,kappa=1.15)$value*2
realProbs[3]<-integrate(dvmises,lower=tail75,upper=pi,Haar=F,kappa=0.52)$value*2
alldf$ScalePdiff<-alldf$Pdiff/(alldf$PMean+alldf$PMedian)

cResFrame<-dcast(ResFrame,n+nu+Sample+Dist~Estimator,value.var="Error")
Midn<-cResFrame[cResFrame$n==100,]
xmax25<-max(Midn[Midn$nu==.25,5:8])
xmax75<-max(Midn[Midn$nu==.75,c(5,7:8)])
Largenu<-ResFrame[ResFrame$nu==.75,]
Largenu50<-Largenu[Largenu$n==50,]
badeggs<-Largenu50[Largenu50$Error>2,]$Sample
Largenu100<-Largenu[Largenu$n==100,]
badeggs<-c(badeggs,Largenu100[Largenu100$Error>2,]$Sample)
beggs<-which(Largenu$Sample %in% badeggs)
Largenu<-Largenu[-beggs,]
Largenu$n<-as.factor(Largenu$n)
levels(Largenu$n)<-c("n = 10","n = 50","n = 100","n = 300")

#Use the data frames to make plots

#This will make Figure 4
qplot(Estimator,Error,geom="boxplot",data=Largen,xlab="",ylab=expression(d[R](bold(S),.)))+
	geom_hline(xintercept=0,colour="gray50")+
	facet_grid(nu~Dist,scales="free",labeller=label_parsed)+
	scale_x_discrete(limits=c("E.Mean","R.Mean","E.Median","R.Median"),breaks=c("E.Mean","R.Mean","E.Median","R.Median"),labels=c(expression(widehat(bold(S))[E]),expression(widehat(bold(S))[R]),expression(widetilde(bold(S))[E]),expression(widetilde(bold(S))[R])))+
	theme(axis.text.x=element_text(size=12,color=1,face='bold'),axis.text.y=element_text(size=12,color=1))

#This will make Figure 5
qplot(n,value,data=mx75,facets=.~variable,geom="path",group=Estimator,linetype=Estimator,lwd=I(1))+
	scale_linetype_manual(values=c(3,4,2,1),labels=my.labels)+coord_equal(ratio=4)+
	theme(legend.text=element_text(size=12),legend.key.width=unit(3,"line"),legend.title=element_text(size=12))+
	geom_hline(yintercept=0,colour="gray50")+
	theme(axis.text.x=element_text(size=12,color=1),axis.text.y=element_text(size=12,color=1))

#This will make Figure 6
ggplot(alldf[alldf$n>299,],aes(Prop,ScalePdiff))+xlab("Proportion Observations in Tail")+
	ylab(expression(frac(d[R](bold(S),widehat(bold(S))[E])-d[R](bold(S),widetilde(bold(S))[E]),d[R](bold(S),widehat(bold(S))[E])+d[R](bold(S),widetilde(bold(S))[E]))))+
	geom_hline(yintercept=0,colour="gray50")+
	stat_smooth(method=lm,formula=y~ns(x,2),fullrange=T,colour=1)+
	geom_point(aes(shape=Dist),alpha=I(.6))+theme_bw()+
	scale_shape_manual(values=c(16,0,17),name="Distribution")+
	theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
	theme(legend.justification=c(0,1),legend.position=c(0,1))

#This will make Figure 7(a)
ggplot(Midn[Midn$nu==.25,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point(alpha=I(.75))+geom_point(aes(E.Median,R.Median),color="grey50",alpha=I(.75))+
	coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
	scale_x_continuous(expression(paste(d[E],"-based estimators")),limits=c(0,xmax25))+
	scale_y_continuous(expression(paste(d[R],"-based estimators")),limits=c(0,xmax25))+
	geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray70")

#This will make Figure 7(b)
ggplot(Midn[Midn$nu==.75,],aes(E.Mean,R.Mean))+facet_grid(.~Dist)+geom_point(alpha=I(.75))+geom_point(aes(E.Median,R.Median),color="grey50",alpha=I(.75))+
	coord_equal(ratio=1)+theme(axis.text.x=element_text(size=14,colour=1),axis.text.y=element_text(size=14,color=1))+
	scale_x_continuous(expression(paste(d[E],"-based estimators")),limits=c(0,xmax75))+
	scale_y_continuous(expression(paste(d[R],"-based estimators")),limits=c(0,xmax75))+
	geom_abline(intercept=0,slope=c(1,0,100000000),colour="gray70")


#This will make Figure 1 of the Supplementary Materials
qplot(Estimator,Error,geom="boxplot",data=Largenu,xlab="",ylab=expression(d[R](bold(S),.)))+
	facet_grid(n~Dist,scales="free")+
	geom_hline(xintercept=0,colour="gray50")+
	scale_x_discrete(limits=c("E.Mean","R.Mean","E.Median","R.Median"),breaks=c("E.Mean","R.Mean","E.Median","R.Median"),labels=c(expression(widehat(bold(S))[E]),expression(widehat(bold(S))[R]),expression(widetilde(bold(S))[E]),expression(widetilde(bold(S))[R])))+
	theme(axis.text.x=element_text(size=12,color=1,face='bold'),axis.text.y=element_text(size=12,color=1))

#This will make Table 1 of the Supplementary Materials
n100ResFrame<-subset(ResFrame,n==100)
n100ses<-ddply(n100ResFrame,.(nu,Dist,Estimator),summarize,mean=mean(Error),SE=sd(Error)/sqrt(length(Error)),RMSE=sqrt(mean((Error)^2)))
levels(n100ses$Estimator)<-c("GeomMedian","GeomMean","ProjMean","ProjMedian")
n100nu25Tab<-cbind(n100ses[c(2,3,1,4),c(1,3:6)],n100ses[c(2,3,1,4)+4,4:6],n100ses[c(2,3,1,4)+8,4:6])
n100nu5Tab<-cbind(n100ses[c(2,3,1,4)+12,c(1,3:6)],n100ses[c(2,3,1,4)+16,4:6],n100ses[c(2,3,1,4)+20,4:6])
n100nu75Tab<-cbind(n100ses[c(2,3,1,4)+24,c(1,3:6)],n100ses[c(2,3,1,4)+28,4:6],n100ses[c(2,3,1,4)+32,4:6])
xtable(rbind(n100nu25Tab,n100nu5Tab,n100nu75Tab),digits=4)

#This will make Table 2
Res$Medians<-Res$HL1Error-Res$MedError
Res$Means<-Res$ML2Error-Res$ArithError
Res$Euclid<-Res$ArithError-Res$MedError
Res$Rieman<-Res$ML2Error-Res$HL1Error
Res$HL1Arith<-Res$HL1Error-Res$ArithError
Res$MedML2<-Res$MedError-Res$ML2Error

ResforTab<-Res[Res$n%in%c(10,100) & Res$nu%in%c(0.25,0.75) & Res$Dist%in%c("Cayley","Fisher"),]
Omits<-which(ResforTab$Means>.22 & ResforTab$ML2Error>=1)
ResforTabOmit<-ResforTab[-Omits,]
results<-matrix(NA,24,10)
dist<-c("Cayley","Fisher")
B<-c(50,100,1000)
n<-c(10,100)
nu<-c(.25,0.75)
rowind<-1
ps<-0
for(i in 1:2){
	for(j in 1:2){
		for(k in 1:2){
			for(l in 1:3){
				ResSub<-ResforTabOmit[ResforTabOmit$nu==nu[j] & ResforTabOmit$n==n[k] & ResforTabOmit$Dist==dist[i],12:17]
				convN<-nrow(ResSub)
				if(B[l]!=1000){
					for(m in 1:100){
						samp<-sample(1:convN,B[l],replace=F)
						ps<-ps+tests(ResSub[samp,])/100
					}
				}else{
					addTo<-1000-convN
					addSamp<-sample(1:convN,addTo,replace=F)
					ResSubtemp<-rbind(ResSub,ResSub[addSamp,])
					ps<-tests(ResSubtemp)
				}
				ps<-round(ps,5)
				results[rowind,]<-c(dist[i],nu[j],n[k],B[l],ps)
				rowind<-rowind+1
				ps<-0
			}
		}
	}
}
results<-data.frame(results)
colnames(results)<-c("Distribution","nu","n","B","Medians","Means","Euclid","Reiman","HL1Arith","MedML2")

for(k in 5:10)
	results[,k]<-as.numeric(as.character(results[,k]))

results[5:10]<-round(results[5:10],4)
xtable(results,digits=4)

#This will make Table 3 of the Supplementary Materials
CaySumL1<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
FisSumL1<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
MisSumL1<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Median-R.Median),sdrbar=sd(E.Median-R.Median)/sqrt(1000),perc=sum(R.Median<E.Median)/1000)
SumL1<-cbind(CaySumL1[,2:5],FisSumL1[,3:5],MisSumL1[,3:5])
xtable(SumL1,digits=4)

#This will make Table 4 of the Supplementary Materials
CaySumL2<-ddply(cResFrame[cResFrame$Dist=="Cayley",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
FisSumL2<-ddply(cResFrame[cResFrame$Dist=="matrix Fisher",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
MisSumL2<-ddply(cResFrame[cResFrame$Dist=="circular-von Mises",],.(nu,n),summarize,rbar=mean(E.Mean-R.Mean),sdrbar=sd(E.Mean-R.Mean)/sqrt(1000),perc=sum(R.Mean<E.Mean)/1000)
SumL2<-cbind(CaySumL2[,2:5],FisSumL2[,3:5],MisSumL2[,3:5])
xtable(SumL2,digits=4)