library(rotations)
library(reshape2)
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
plot(Rs,center=id.SO3,col=1)
plot(Rs,center=id.SO3,col=2)
plot(Rs,center=id.SO3,col=3)


#This will produce a plot similar to Figure 3(b) and Figures 3(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rfisher,nu=0.25)
plot(Rs,center=id.SO3,col=1)
plot(Rs,center=id.SO3,col=2)
plot(Rs,center=id.SO3,col=3)



#This will produce a plot similar to Figure 3(c) and Figures 4(a)-(c) in Supplementary Matrials
Rs<-ruars(100,rvmises,nu=0.25)
plot(Rs,center=id.SO3,col=1)
plot(Rs,center=id.SO3,col=2)
plot(Rs,center=id.SO3,col=3)
