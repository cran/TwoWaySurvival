"plot.OneWaySurvfit" <-
function(x,...)
{
factor.names<-x$factor.names
grid.frame<-x$grid.frame
varying.frame<-x$varying.frame
deviation.frame<-x$deviation.frame
p<-x$p

assign("grid.t",grid.frame$grid.t)
for (i in 1:dim(varying.frame)[2]) assign(names(varying.frame)[i],varying.frame[,i])
for (i in 1:dim(deviation.frame)[2]) assign(names(deviation.frame)[i],deviation.frame[,i])


las<-1;cex.main<-0.9;tcl<- -0.1;cex.lab<-0.7;cex.axis<-0.7;lwd<-1


y.range.t.baseline<-range(c(varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline,varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline))

plot(grid.t,varying.frame$alpha.t.Baseline,xlab="Duration time (t)",ylab="",cex=0.1,main="Baseline",ylim=y.range.t.baseline,...)
vector.minus<-varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline
vector.plus<-varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline
polygon(cbind(c(grid.t,grid.t[length(grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.t,varying.frame$alpha.t.Baseline,lwd=lwd)
lines(grid.t,varying.frame$alpha.t.Baseline-deviation.frame$deviation.t.Baseline,cex=0.08,col=3)
lines(grid.t,varying.frame$alpha.t.Baseline+deviation.frame$deviation.t.Baseline,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.t,varying.frame$alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)



if (dim(x$varying.frame)[2] > 1)
  {
#ranges for y in plots
y.range.t<-range(c(unlist(mget(paste("varying.frame$alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))-unlist(mget(paste("deviation.frame$deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1))),unlist(mget(paste("varying.frame$alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))+unlist(mget(paste("deviation.frame$deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))))

y.range<-range(y.range.t)  
    
  
for (k in 1:p)
{
plot(grid.t,get(paste("varying.frame$alpha.t.",factor.names[k],sep="")),xlab="Duration time (t)",ylab="",cex=0.1,main=factor.names[k],ylim=y.range,...)
vector.minus<-get(paste("varying.frame$alpha.t.",factor.names[k],sep=""))-get(paste("deviation.frame$deviation.t.",factor.names[k],sep=""))
vector.plus<-get(paste("varying.frame$alpha.t.",factor.names[k],sep=""))+get(paste("deviation.frame$deviation.t.",factor.names[k],sep=""))
polygon(cbind(c(grid.t,grid.t[length(grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.t,get(paste("varying.frame$alpha.t.",factor.names[k],sep="")),lwd=lwd)
lines(grid.t,get(paste("varying.frame$alpha.t.",factor.names[k],sep=""))-get(paste("deviation.frame$deviation.t.",factor.names[k],sep="")),cex=0.08,col=3)
lines(grid.t,get(paste("varying.frame$alpha.t.",factor.names[k],sep=""))+get(paste("deviation.frame$deviation.t.",factor.names[k],sep="")),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.t,alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)
}
}


}

