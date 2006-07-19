"plot.TwoWaySurvfit" <-
function(x,...)
{
factor.names<-x$factor.names
grid.frame<-x$grid.frame
varying.frame<-x$varying.frame
deviation.frame<-x$deviation.frame


assign("grid.t",grid.frame$grid.t)
assign("grid.b",grid.frame$grid.b)
for (i in 1:dim(varying.frame)[2]) assign(names(varying.frame)[i],varying.frame[,i])
for (i in 1:dim(deviation.frame)[2]) assign(names(deviation.frame)[i],deviation.frame[,i])


plot(grid.t,alpha.t.Baseline,xlab="Duration time (t)",ylab="",cex=0.5,main="Baseline",ylim=range(c(alpha.t.Baseline-deviation.t.Baseline,alpha.t.Baseline+deviation.t.Baseline)),...)
vector.minus<-alpha.t.Baseline-deviation.t.Baseline
vector.plus<-alpha.t.Baseline+deviation.t.Baseline
polygon(cbind(c(grid.t,grid.t[length(grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.t,alpha.t.Baseline,lwd=3)
lines(grid.t,alpha.t.Baseline-deviation.t.Baseline,cex=0.08,col=3)
lines(grid.t,alpha.t.Baseline+deviation.t.Baseline,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.t,alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)


plot(grid.b,alpha.b.Baseline,xlab="Entry time (b)",ylab="",cex=0.5,main="Baseline",axes=TRUE,ylim=range(c(alpha.b.Baseline-deviation.b.Baseline,alpha.b.Baseline+deviation.b.Baseline)),...)
vector.minus<-alpha.b.Baseline-deviation.b.Baseline
vector.plus<-alpha.b.Baseline+deviation.b.Baseline
polygon(cbind(c(grid.b,grid.b[length(grid.b):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.b,alpha.b.Baseline,lwd=3)
lines(grid.b,alpha.b.Baseline-deviation.b.Baseline,cex=0.08,col=3)
lines(grid.b,alpha.b.Baseline+deviation.b.Baseline,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.b,alpha.b.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)


if (dim(x$varying.frame)[2] > 2)
  {
#ranges for y in plots
y.range.t<-range(c(unlist(mget(paste("alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))-unlist(mget(paste("deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1))),unlist(mget(paste("alpha.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))+unlist(mget(paste("deviation.t.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))))

y.range.b<-range(c(unlist(mget(paste("alpha.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))-unlist(mget(paste("deviation.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1))),unlist(mget(paste("alpha.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))+unlist(mget(paste("deviation.b.",factor.names[1:length(factor.names)],sep=""),envir=as.environment(-1)))))

y.range<-range(y.range.t,y.range.b)  
    
  
for (k in 1:x$p)
{
#plot.t
plot(grid.t,get(paste("alpha.t.",factor.names[k],sep="")),xlab="Duration time (t)",ylab="",cex=0.5,main=factor.names[k],ylim=y.range,...)
vector.minus<-get(paste("alpha.t.",factor.names[k],sep=""))-get(paste("deviation.t.",factor.names[k],sep=""))
vector.plus<-get(paste("alpha.t.",factor.names[k],sep=""))+get(paste("deviation.t.",factor.names[k],sep=""))
polygon(cbind(c(grid.t,grid.t[length(grid.t):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.t,get(paste("alpha.t.",factor.names[k],sep="")),lwd=3)
lines(grid.t,get(paste("alpha.t.",factor.names[k],sep=""))-get(paste("deviation.t.",factor.names[k],sep="")),cex=0.08,col=3)
lines(grid.t,get(paste("alpha.t.",factor.names[k],sep=""))+get(paste("deviation.t.",factor.names[k],sep="")),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.t,alpha.t.Baseline,type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)


#plot.b
plot(grid.b,get(paste("alpha.b.",factor.names[k],sep="")),xlab="Entry time (b)",ylab="",cex=0.5,main=factor.names[k],axes=TRUE,ylim=y.range,...)
vector.minus<-get(paste("alpha.b.",factor.names[k],sep=""))-get(paste("deviation.b.",factor.names[k],sep=""))
vector.plus<-get(paste("alpha.b.",factor.names[k],sep=""))+get(paste("deviation.b.",factor.names[k],sep=""))
polygon(cbind(c(grid.b,grid.b[length(grid.b):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.b,get(paste("alpha.b.",factor.names[k],sep="")),lwd=3)
lines(grid.b,get(paste("alpha.b.",factor.names[k],sep=""))-get(paste("deviation.b.",factor.names[k],sep="")),cex=0.08,col=3)
lines(grid.b,get(paste("alpha.b.",factor.names[k],sep=""))+get(paste("deviation.b.",factor.names[k],sep="")),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.b,get(paste("alpha.b.",factor.names[k],sep="")),type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)
}
}


}

