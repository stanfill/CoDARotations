#To install the 'rotations' package, which is not currently available on CRAN, use the following code.
#In the near future, this package will include the dataset used in the 
#Data Applications section of this paper

library(devtools)
install_github("rotations","heike")
library(rotations)

#Load the necessary libraries and create the 'tests' function used to create 
#Table 2 of the Supplementary Materials
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
##In this section we demonstrate how to perform a 
##simulation study similar to the one performed in
##the manuscript. Repeat this for a variety of sample sizes (n)
##and circular variances (nu)
######################

#Generate a random sample from the UARS distribution of interest
n<-10
nu<-0.25
Rs<-ruars(n,rcayley,S=id.SO3,nu=nu)   	 #Cayley-UARS distribution
#Rs<-ruars(n,rfisher,S=id.SO3,nu=nu)		 #Fisher-UARS distribution
#Rs<-ruars(n,rvmises,S=id.SO3,nu=nu)		 #circular von Mises-UARS distribution

ShatE<-mean(Rs,type='projected')     	#Projected Mean
ShatR<-mean(Rs,type='intrinsic')			#Intrinsic Mean
StildeE<-median(Rs,type='projected')	#Projected Median
StildeR<-median(Rs,type='intrinsic')	#Intrinsic Median

#Find the geodesic distance between the true mean (id.SO3) and each estimate
errorShatE<-dist(ShatE,id.SO3,method='intrinsic')
errorShatR<-dist(ShatR,id.SO3,method='intrinsic')
errorStildeE<-dist(StildeE,id.SO3,method='intrinsic')
errorStildeR<-dist(StildeR,id.SO3,method='intrinsic')

######################
##In this section the plots in the Section 4 (Simulation Study) 
## and some of the Supplementary Materials plots are made
######################

r<-seq(-pi,pi,length=1000)

cayKap<-c(10,4,2,.065)
fishKap<-c(3.165,1.711,1.1564,.175)
misesKap<-c(2.40, 1.159, 0.5164)

denData<-data.frame(r=r,Fisher=dfisher(r,fishKap[1]),Mises=dvmises(r,misesKap[1]),Cayley=dcayley(r,cayKap[1]),nu=0.25)
denData<-rbind(denData,data.frame(r=r,Fisher=dfisher(r,fishKap[2]),Mises=dvmises(r,misesKap[2]),Cayley=dcayley(r,cayKap[2]),nu=0.5))
denData<-rbind(denData,data.frame(r=r,Fisher=dfisher(r,fishKap[3]),Mises=dvmises(r,misesKap[3]),Cayley=dcayley(r,cayKap[3]),nu=0.75))

denData<-melt(denData,id=c("r","nu"))
colnames(denData)<-c("r","Variance","Density","Value")

levels(denData$Density)<-c("matrix Fisher","circ-von Mises","Cayley")
denData$Density<-factor(denData$Density,levels=levels(denData$Density)[c(3,1,2)])

#This is Figure 2(a) of the manuscript
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,25),ylab="f(r)")+geom_abline(intercept=0,slope=0,colour="gray50")+
	theme(legend.position='NONE')+scale_linetype_manual(values=c(1,12,3))+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))+geom_vline(xintercept=0)


#This is Figure 2(b) of the manuscript
Boxx<-c(1.75,1.75,2.5,2.5,1.75)
Boxy<-c(0,.25,.25,0,0)
qplot(r,Value,data=denData[denData$Variance==.75,],geom="blank",ylim=c(0,2.25),ylab="f(r)")+
	scale_linetype_manual(values=c(1,12,3))+theme(legend.position="none")+geom_path(aes(x=Boxx,y=Boxy),lwd=I(1.25),colour="gray")+
	geom_line(aes(x=r,y=Value,linetype=Density),lwd=I(1.25))+geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))

#This is Figure 2(c) of the manuscript
qplot(r,Value,data=denData[denData$Variance==.75,],linetype=Density,geom="line",lwd=I(1.25),ylim=c(0,.25),xlim=c(1.75,2.5),ylab="f(r)")+
	scale_linetype_manual(values=c(1,12,3))+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_rect(fill="white",linetype=0))+
	theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=16,face="bold"),legend.key.width=unit(3,"line"))+
	geom_abline(intercept=0,slope=c(0,100000),colour="gray50")+
	theme(axis.text.x=element_text(size=14,colour=1))+theme(axis.text.y=element_text(size=14,colour=1))


#This will produce a plot similar to Figure 3(a) and Figures 2(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rcayley,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=2) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=3) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")


#This will produce a plot similar to Figure 3(b) and Figures 3(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rfisher,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=2) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
plot(Rs,center=id.SO3, col=3) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")



#This will produce a plot similar to Figure 3(c) and Figures 4(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rvmises,nu=0.25)
plot(Rs,center=id.SO3, col=1) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5)) + theme(legend.position="none")
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
mx75$variable<-factor(mx75$variable,labels=c("Mean Error","RMSE"))
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
	geom_hline(yintercept=0,colour="gray50")+ylab("")+
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

#This will make Table 2 of the Supplementary Materials
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



######################
##In this section we demonstrate how to reproduce the results in the Data
##Application section of the paper.  To gain access to the data, please email
##Dr. Melissa Bingham, based in University of Wisconsin - La Crosse as of this
##writing.  In the near future this data will be included in the rotations
##package.
######################

load("datasetnickel.RData")
require(plyr)
dat.out <- adply(data, .margins= c(1,3), function(x) {
	as.vector(x)
})
names(dat.out)[c(1,2)] <- c("rep", "location")
dat.out$xpos <- xpos[dat.out$location]
dat.out$ypos <- ypos[dat.out$location]

checks <- adply(dat.out, .margins=1, function(x) {
	is.SO3(unlist(x[3:11]))
})
dat.out$check <- checks$V1

dat.ests <- dlply(dat.out, .(location), function(x) {
	res <- na.omit(x)
	res <- subset(res, check==TRUE)
	
	n <- nrow(res) 
	SE2 <- SR2 <- SE1 <- SR1 <- NULL
	if (n == 1) {
		R <- as.SO3(matrix(unlist(res[1,3:11]), 3, 3))
		SE2 <- SR2 <- SE1 <- SR1 <- R
	} else if (n > 0) {
		rots <- as.SO3(as.matrix(res[,3:11]))
		SE2 <- mean(rots)
		SR2 <- mean(rots, type='intrinsic')
		SE1 <- median(rots)
		SR1 <- median(rots, type='intrinsic')
	}
	location <- as.numeric(as.character(unique(x$location)))
	return(list(location=location, n=n, SE2=SE2, SR2=SR2, SE1=SE1, SR1=SR1))
})

## find distances between estimators and angles to identity for each
loc.stats <- ldply(dat.ests, function(x) {  
	location <- as.numeric(as.character(unique(x$location)))
	if (x$n > 0)
		data.frame(location=x$location, n=x$n, 
							 dE1=angle(x$SE1), dE2=angle(x$SE2),
							 dR1=angle(x$SR1), dR2=angle(x$SR2),
							 dE=dist(x$SE1, x$SE2, method="intrinsic", p=1),
							 dR=dist(x$SR1, x$SR2, method="intrinsic", p=1),
							 d1=dist(x$SE1, x$SR1, method="intrinsic", p=1),
							 d2=dist(x$SE2, x$SR2, method="intrinsic", p=1)
		)
})

loc.stats$xpos <- xpos[loc.stats$location]
loc.stats$ypos <- ypos[loc.stats$location]
loc.stats <- adply(loc.stats, .margin=1, transform, sd=sd(c(dE1, dE2, dR1, dR2)))
idx <- which.max(loc.stats$dEdegree)

require(ggplot2)
#This will make Figure 8(a)
d <- ggplot(loc.stats, aes(xpos, ypos, color=dE1))
d2 <- d + geom_point(size=2.25) + scale_colour_gradient(expression(d[R](tilde(S)[E], I["3x3"])), low="grey99", high="grey10", limits=c(0, pi), 
					breaks=c( pi/4, pi/2, 3*pi/4), labels=expression( pi/4, pi/2, 3*pi/4)) + theme_bw() + xlab("") + ylab("") + coord_equal() + 
					scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), 
					labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
					scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
					geom_point(shape="o", colour="yellow", size=5, data=loc.stats[idx,])  + theme(plot.margin=unit(rep(0,4), "lines"))
d2

label <- "in degrees"
mains <- bquote(.(parse(text=paste("d[R](tilde(S)[E], hat(S)[E]", quote("in degrees"), sep="\n"))) )
mains <- expression(paste(d[R](tilde(S)[E], hat(S)[E]), "  (in ",degree,")  "))
loc.stats$dEdegree <- loc.stats$dE*180/pi
d <- ggplot(loc.stats, aes(xpos, ypos, color=dEdegree))

#This will make Figure 8(b)
d + geom_point(size=2.2) + scale_colour_gradient(mains, low="white", high="grey10", trans = "sqrt", breaks=c( 0, 0.5, 5, 20))+ theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m))  + geom_point(shape="o", colour="yellow", size=5, data=loc.stats[idx,]) + theme(plot.margin=unit(rep(0,4), "lines"))


#Select the sample we would like to illustrate more cleanly
idx <- which.max(loc.stats$dEdegree)
datdf <- subset(dat.out, location %in% loc.stats[idx,]$location)

#This will make Table 4
xtable(datdf[,c(1,3:12)],digits=3)

subdf <- as.matrix(subset(datdf, check)[,3:11])
pmed <- median(as.SO3(subdf), type="projected")
gmed <- median(as.SO3(subdf), type="intrinsic")
pmean <- mean(as.SO3(subdf), type="projected")
gmean <- mean(as.SO3(subdf), type="intrinsic")
ests <- data.frame(rbind(as.vector(pmean), as.vector(pmed), as.vector(gmean), as.vector(gmed)))
ests$ID <- 1:4
require(reshape2)
em <- melt(ests, id.var="ID")
em$variable <- as.numeric(gsub("X", "", as.character(em$variable)))
labels <- c(expression(hat(S)[E]), expression(tilde(S)[E]), expression(hat(S)[R]), expression(tilde(S)[R]))
require(RColorBrewer)
cols <- brewer.pal(4, "Paired")
n <- nrow(datdf)
cay <- rcayley(n, kappa=4000)
jitdf <- genR(cay)
for (i in 1:nrow(datdf)) {
	R <- matrix(unlist(datdf[i,3:11]),3,3)
	S <- matrix(unlist(jitdf[i,]),3,3)
	jitdf[i,] <- unlist(R %*% S)
}
class(jitdf) <- "matrix"
jitdf <- data.frame(jitdf)
names(jitdf) <- expression(x[11],x[12],x[13],x[21],x[22],x[23],x[31],x[32],x[33])

#This will create Figure 9(a)
ggpcp(jitdf) + geom_line() + ylim(c(-1,1)) + 
	scale_x_discrete("", labels = expression(x[11],x[12],x[13],x[21],x[22],x[23],x[31],x[32],x[33])) +
	geom_line(aes(x=variable, y=value, colour=factor(ID), group=ID), data=subset(em, ID %in% c(1,2)), size=1, inherit.aes=F) + scale_colour_manual("Estimates", values=cols[4:3], labels=labels[1:2]) + theme_bw()

#This will create Figure 9(b) and the other axes 
mid <- median(as.SO3(subdf), type="projected") 
plot.SO3(subdf, center=mid, col=1, show_estimates=c("proj.mean", "proj.median")) + theme(legend.position="none")+ scale_colour_manual("Estimates", values=cols[4:3], labels=labels[1:2])
plot.SO3(as.matrix(jitdf), center=mid, col=2, show_estimates=c("proj.mean", "proj.median")) + theme(legend.position="none")+ scale_colour_manual("Estimates", values=cols[4:3], labels=labels[1:2])
plot.SO3(subdf, center=mid, col=3, show_estimates=c("proj.mean", "proj.median"))+ theme(legend.position="none")+ scale_colour_manual("Estimates", values=cols[4:3], labels=labels[1:2])

