"plot.TwoWaySurvfit" <-
function(x,...)
{
factor.names<-x$factor.names
grid.frame<-x$grid.frame
varying.frame<-x$varying.frame
deviation.frame<-x$deviation.frame
p<-x$p

attach(varying.frame);attach(deviation.frame)


las<-1;cex.main<-0.9;tcl<- -0.1;cex.lab<-0.7;cex.axis<-0.7;lwd<-1


#plot t
y.range.t.baseline<-range(c(varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline,varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline))

plot(grid.frame$grid.t,varying.frame$alpha.t.Baseline,xlab="Duration time (t)",ylab="",cex=0.1,main="Baseline",ylim=y.range.t.baseline,...)
vector.minus<-varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline
vector.plus<-varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline
polygon(cbind(c(grid.frame$grid.t,grid.frame$grid.t[length(grid.frame$grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame$grid.t,varying.frame$alpha.t.Baseline,lwd=lwd)
lines(grid.frame$grid.t,varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline,cex=0.08,col=3)
lines(grid.frame$grid.t,varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame$grid.t,varying.frame$alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)


#plot b
y.range.b.baseline<-range(c(varying.frame$alpha.b.Baseline-deviation.frame$deviation.b.Baseline,varying.frame$alpha.b.Baseline+deviation.frame$deviation.b.Baseline))

plot(grid.frame$grid.b,varying.frame$alpha.b.Baseline,xlab="Entry time (b)",ylab="",cex=0.1,main="Baseline",axes=TRUE,ylim=y.range.b.baseline,...)
vector.minus<-varying.frame$alpha.b.Baseline-deviation.frame$deviation.b.Baseline
vector.plus<-varying.frame$alpha.b.Baseline+deviation.frame$deviation.b.Baseline
polygon(cbind(c(grid.frame$grid.b,grid.frame$grid.b[length(grid.frame$grid.b):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame$grid.b,varying.frame$alpha.b.Baseline,lwd=lwd)
lines(grid.frame$grid.b,varying.frame$alpha.b.Baseline-deviation.frame$deviation.b.Baseline,cex=0.08,col=3)
lines(grid.frame$grid.b,varying.frame$alpha.b.Baseline+deviation.frame$deviation.b.Baseline,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame$grid.b,varying.frame$alpha.b.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)



if (dim(x$varying.frame)[2] > 2)
  {
#ranges for y in plots
y.range.t<-range(c(unlist(mget(paste("alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("varying.frame")))-unlist(mget(paste("deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("deviation.frame"))),unlist(mget(paste("alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("varying.frame")))+unlist(mget(paste("deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("deviation.frame")))))

y.range.b<-range(c(unlist(mget(paste("alpha.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("varying.frame")))-unlist(mget(paste("deviation.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("deviation.frame"))),unlist(mget(paste("alpha.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("varying.frame")))+unlist(mget(paste("deviation.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment("deviation.frame")))))

y.range<-range(y.range.t,y.range.b)  
    
  
for (k in 1:p)
{
#plot.t
plot(grid.frame$grid.t,get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame"),xlab="Duration time (t)",ylab="",cex=0.1,main=factor.names[k],ylim=y.range,...)
vector.minus<-get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame")-get(paste("deviation.t.",factor.names[k],sep=""),pos="deviation.frame")
vector.plus<-get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame")+get(paste("deviation.t.",factor.names[k],sep=""),pos="deviation.frame")
polygon(cbind(c(grid.frame$grid.t,grid.frame$grid.t[length(grid.frame$grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame$grid.t,get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame"),lwd=lwd)
lines(grid.frame$grid.t,get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame")-get(paste("deviation.t.",factor.names[k],sep=""),pos="deviation.frame"),cex=0.08,col=3)
lines(grid.frame$grid.t,get(paste("alpha.t.",factor.names[k],sep=""),pos="varying.frame")+get(paste("deviation.t.",factor.names[k],sep=""),pos="deviation.frame"),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame$grid.t,alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)


#plot.b
plot(grid.frame$grid.b,get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame"),xlab="Entry time (b)",ylab="",cex=0.1,main=factor.names[k],axes=TRUE,ylim=y.range,...)
vector.minus<-get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame")-get(paste("deviation.b.",factor.names[k],sep=""),pos="deviation.frame")
vector.plus<-get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame")+get(paste("deviation.b.",factor.names[k],sep=""),pos="deviation.frame")
polygon(cbind(c(grid.frame$grid.b,grid.frame$grid.b[length(grid.frame$grid.b):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame$grid.b,get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame"),lwd=lwd)
lines(grid.frame$grid.b,get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame")-get(paste("deviation.b.",factor.names[k],sep=""),pos="deviation.frame"),cex=0.08,col=3)
lines(grid.frame$grid.b,get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame")+get(paste("deviation.b.",factor.names[k],sep=""),pos="deviation.frame"),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame$grid.b,get(paste("alpha.b.",factor.names[k],sep=""),pos="varying.frame"),type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)
}
}

detach(varying.frame);detach(deviation.frame)

}

