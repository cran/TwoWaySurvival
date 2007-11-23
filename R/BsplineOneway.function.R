BsplineOneway<-function(data.set=data.set,control=control)
{

#############################################################
##One Way Hazard Model With Varying Coefficients#############
#############################################################

#############################################################  
#modelling with B-splines (all coefficients are penalized)###
#############################################################  

library(splines)


if (names(data.set)[3]=="Intercept") data.set<-data.set[,-3]

N<-length(data.set[,1])   #total number of observations in data.set



#################################################
#knots and order of B-spline Bases###############
#################################################
default.knots<-function(x,num.knots)
{
bs.order<-4   
bs.degree<-bs.order-1   

if (missing(num.knots)) num.knots<-20

knots.inner<-seq(min(x),max(x),le=num.knots)
intervall<-knots.inner[length(knots.inner)]-knots.inner[1]

#aquidistant periodic boundary knots
knots.left<-knots.inner[(length(knots.inner)-bs.degree):(length(knots.inner)-1)]-intervall 
knots.right<-knots.inner[2:(bs.degree+1)]+intervall
knots<-c(knots.left,knots.inner,knots.right)

#coincident boundary.knots
#knots.left<-rep(knots.inner[1],bs.degree)  
#knots.right<-rep(knots.inner[length(knots.inner)],bs.degree)
#knots<-c(knots.left,knots.inner,knots.right)

knots

} #end default.knots



#################################################
#create artificial poisson data##################
#################################################
survival.to.poisson <- function (time=time, status=status, x=NULL)
{
N <- length(time)
event.time <- time[status==1]

#create the grid
#grid <- c(unique(sort(event.time)),max(time+1))  #all event.time points
#grid<-c(sort(sample(unique(time[status==1]),control$number.int)),max(time+1)) #just a sample of event.time points
grid<-round(quantile(unique(time),prob=seq(0,1,le=control$number.int)),max(time+1)) #quantile of specified length

m <- length(grid)
grid.minus1 <- c(0,grid[1:(m-1)])
grid.plus1 <- c(grid[2:m],max(time)+2)

#create poisson data
Yt.list<-sapply(seq(time),FUN=function(i) c(rep(0,sum(grid<time[i])),status[i]))
Mt.list<-sapply(seq(time),FUN=function(i) sum(grid<time[i])+1)
Mt<-unlist(Mt.list)

#inflate the matrix with variables accordingly 
Xt<-sapply(seq(Mt),FUN=function(i) kronecker(matrix(1,Mt[i],1),t(matrix(as.matrix(x[i,])))))
Xt<-eval(parse(text=paste("rbind(",paste("Xt[[",1:length(Xt),"]]",sep="",collapse=","),")",sep="")))

#offset parameter
Ot.list<-sapply(seq(Mt),FUN=function(i) log(0.5*(apply(cbind(grid.plus1[1:Mt[i]],rep(time[i],Mt[i])),MARGIN=1,FUN=min)- apply(cbind(grid.minus1[1:Mt[i]],rep(time[i],Mt[i])),MARGIN=1,FUN=min))))
#beachte: diese berechnung von offset weicht am rechten rand von der angegebnen formel ab  

return(list(grid=grid,y.list=Yt.list,m.list=Mt.list,o.list=Ot.list,x=Xt))
} #end survival.to.poisson



#################################################
#create MP-inverse of a matrix###################
#################################################
ginverse<-function(X,tol=1e-100)
{
Xsvd<-svd(X,LINPACK=TRUE)
if (is.complex(X)) Xsvd$u<-Conj(Xsvd$u)
Positive<-Xsvd$d > max(tol * Xsvd$d[1], 0)
if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
else if (!any(Positive)) array(0, dim(X)[2:1])
else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}




#create a model.matrix for covariables
if (dim(data.set)[2] > 2) model.matrix.x<-model.matrix(as.formula(paste("~",paste(names(data.set)[3:length(names(data.set))],collapse="+"))),data=data.set) else model.matrix.x<-model.matrix(~-1+rep(1,N))



help.surv<-survival.to.poisson(data.set$time,data.set$status,model.matrix.x)  
grid.points<-help.surv$grid   #t-points, at which the integral is being approximated
#for each observation there is an element of the list with accordingly:
y.poisson.list<-help.surv$y.list  #poisson.data
y.poisson<-unlist(y.poisson.list)
offset.list<-help.surv$o.list     #offset.parameter
offset.model<-unlist(offset.list)
m.list<-help.surv$m.list          #number of poisson data




Design.variables<-help.surv$x
if (dim(data.set)[2] > 2) colnames(Design.variables)<-c("Intercept",colnames(model.matrix.x)[2:length(colnames(model.matrix.x))]) else colnames(Design.variables)<-c("Intercept")


#number of varaibles (including factor levels by factor variables) in the defined design matrix
p<-ncol(Design.variables)-1 



#knots and order of B-splines
bs.order<-4   
bs.degree<-bs.order-1 
if (length(grep("num.knots.t",names(control)))==0) knots.t<-default.knots(data.set$time) else knots.t<-default.knots(data.set$time,num.knots.t)



#define survival time for expanded data set, and the b-spline bases
Z.time.list<-sapply(data.set$time,FUN=function(x) c(grid.points[grid.points<x],x))
Z.time<-unlist(Z.time.list)
Basis.t<-splineDesign(knots=knots.t,ord=bs.order,x=Z.time)
K.t<-dim(Basis.t)[2]



#build Difference matrix
Diff.2.t<-matrix(0,nrow=K.t-2,ncol=K.t)    #diff.penalty of 2.order
for (i in 1:(K.t-2)) Diff.2.t[i,i:(i+2)]<-c(1,-2,1)
Penalty.matrix.t<-crossprod(Diff.2.t)



#initialize penalty parameters through variances of random effects
variance.penalty.t.list<-list()

if (length(grep("start.penalty.t",names(control)))==0)
{
if (p > 0) 
  {
for (i in 1:(p+1)) variance.penalty.t.list[[i]]<-10e-1
names(variance.penalty.t.list)<-c("variance.t.Baseline",paste("variance.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
  } else
  {
variance.penalty.t.list[[1]]<-10e-1
names(variance.penalty.t.list)<-c("variance.t.Baseline")
  }

} else

{
if (p > 0) 
  {
for (i in 1:(p+1)) variance.penalty.t.list[[i]]<-1/control$start.penalty.t
names(variance.penalty.t.list)<-c("variance.t.Baseline",paste("variance.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
  } else
  {
variance.penalty.t.list[[1]]<-1/control$start.penalty.t
names(variance.penalty.t.list)<-c("variance.t.Baseline")
  }
}

variance.penalty.t<-unlist(variance.penalty.t.list)                  
variance.penalty<-c(variance.penalty.t)




#initialize random effects
u.t.list<-list()
if (p > 0)
  {
for (i in 1:(p+1)) u.t.list[[i]]<-rep(1,K.t)
names(u.t.list)<-c("u.t.Baseline.",paste("u.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))                         
} else
  {
u.t.list[[1]]<-rep(1,K.t);names(u.t.list)<-c("u.t.Baseline.")
}

u.t<-unlist(u.t.list)



#combine in theta vector
teta.t.Baseline<-u.t[1:length(u.t.list[[1]])]
teta.t.list<-list()
teta.t.list[[1]]<-teta.t.Baseline
if (p > 0) for (i in 1:p) teta.t.list[[i+1]]<-u.t[(i*length(u.t.list[[i]])+1):((i+1)*length(u.t.list[[i]]))]
teta.t<-unlist(teta.t.list)



#create penalty matrix
Penalty.matrix.teta.t<-matrix(0,nrow=length(teta.t),ncol=length(teta.t))
Penalty.matrix.teta.t[1:K.t,1:K.t]<-(1/variance.penalty.t[1])*Penalty.matrix.t
if (p > 0)
for (k in 1:p) Penalty.matrix.teta.t[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)]<-(1/variance.penalty.t[k+1])*Penalty.matrix.t
                


#construct design matrix for the model
if (p > 0)
{  
variables.t<-Design.variables[,2:ncol(Design.variables),drop=FALSE]
Design.matrix.t<-cbind(Basis.t,eval(parse(text=paste("cbind(",paste("variables.t[,",1:p,"]*Basis.t",sep="",collapse=","),")"))))
} else
Design.matrix.t<-Basis.t




#remove not more used objects
rm(Z.time.list,Z.time,help.surv,Basis.t,model.matrix.x)
for (i in 1:30) gc()



###############################################################################
#initialize list object for penalty variances from following iterations########
###############################################################################
variance.penalty.list<-list()  
variance.penalty.list[[1]]<-variance.penalty   #first element of the list is the initial variance vector




status<-data.set$status


log.lik.margin.t.list<-list()   #saves the values of marginal likelihood for all iterations


for (iter.epoch in 1:control$niter.epoch)
 {
if (control$print.epoch) cat("epoch= ", iter.epoch, "\n\n")


 
#likelihood for alpha.t
offset.tij.list<-lapply(1:N,FUN=function(i) offset.list[[i]])
unlist.offset.tij.list<-unlist(offset.tij.list)
log.lik.t.pen<-function(x)
{
log.lik.t<-sum(unlist(y.poisson.list)*Design.matrix.t%*%x-exp(Design.matrix.t%*%x+unlist.offset.tij.list))
log.lik.t-0.5*t(x)%*%Penalty.matrix.teta.t%*%x
}


#score fot theta.t
S.t<-crossprod(Design.matrix.t,as.vector(y.poisson-exp(Design.matrix.t%*%teta.t+unlist.offset.tij.list)))
S.t.pen<-S.t-Penalty.matrix.teta.t%*%teta.t
#Fisher.type.Matrix(observed) 
I.t<-crossprod(Design.matrix.t,Design.matrix.t*as.vector(exp(Design.matrix.t%*%teta.t+unlist.offset.tij.list)))
I.t.pen<-I.t+Penalty.matrix.teta.t
#add a small positiv value of the identity matrix (for the case not being invertaible)
delta<-1e-05
tau<-max(0,delta-min(Re(eigen(I.t.pen)$values)))
I.t.pen<-I.t.pen+tau*diag(1,ncol(I.t.pen))
inverse.t <- qr.solve(I.t.pen,tol=1e-100)

teta.t.old<-teta.t
delta.teta.t<-inverse.t%*%S.t.pen


###########################################################
#optimization goes by backtracking.line.search procedure###
###########################################################
#fix parameters
alpha.backtracking<-0.01
beta.backtracking<-0.5
#optimal value of the step s
s<-1
while (-log.lik.t.pen(teta.t.old+s*delta.teta.t) > -log.lik.t.pen(teta.t.old)+alpha.backtracking*s*crossprod(-S.t.pen,delta.teta.t)) {s<-s*beta.backtracking}

teta.t.new<-teta.t.old+s*delta.teta.t
names(teta.t.new)<-names(teta.t)    
teta.t<-teta.t.new
names(teta.t)<-names(teta.t.new)    

#discent.vector.t<-crossprod(-S.t.pen,inverse.t%*%S.t.pen)


#update of random effects
u.t<-teta.t.new[grep("u.t",names(teta.t),fixed=TRUE)]


#print estimates
if (control$print.epoch)
{       
print(u.t)
cat("","\n\n")
}





###################################################################
#optimization procedure for penalties##############################
#####goes by fix.point iteration###################################
###################################################################



#define first the submatrices of the (penalized) fisher.information, corresponding to the random coefficients
index.vector.t<-1:K.t  #for u.t.Baseline
if (p > 0) for (k in 1:p) index.vector.t<-c(index.vector.t,(K.t+(k-1)*K.t+1):(K.t+k*K.t)) #for covariables



#browser()


################################################################################
#chose optimization routine for smoothing.parameters; only fix.point goes#######
################################################################################
if (control$method !="fix") stop("only method='fix' allowed") else
{

for (iter.penalty in 1:control$niter.penalty)
{  
if (control$print.penalty) cat("iteration.penalty= ",iter.penalty,"\n\n") 

  
I.t.u.t<-I.t.pen[index.vector.t,index.vector.t]    #penalized fisher.type matrix for t
inverse.I.t.u.t<-qr.solve(I.t.u.t,tol=1e-100)


variance.penalty.t.old<-variance.penalty.t
variance.penalty.old<-c(variance.penalty.t.old)


variance.t.Baseline<-as.vector((crossprod(u.t[grep("Baseline",names(u.t),fixed=TRUE)],Penalty.matrix.t%*%u.t[grep("Baseline",names(u.t),fixed=TRUE)]) + sum(diag(inverse.I.t.u.t[1:K.t,1:K.t]%*%Penalty.matrix.t[1:K.t,1:K.t])))/(K.t-2))
variance.penalty.t[1]<-variance.t.Baseline
if ( p > 0)  #number of covariables
for (k in 1:p) variance.penalty.t[k+1]<-as.vector((crossprod(u.t[(k*K.t+1):((k+1)*K.t)],Penalty.matrix.t%*%u.t[(k*K.t+1):((k+1)*K.t)]) + sum(diag(inverse.I.t.u.t[(k*K.t+1):((k+1)*K.t),(k*K.t+1):((k+1)*K.t)]%*%Penalty.matrix.t)))/(K.t-2))
#update of penalty matrix
Penalty.matrix.teta.t[1:K.t,1:K.t]<-(1/variance.penalty.t[1])*Penalty.matrix.t
if (p > 0) for (k in 1:p) Penalty.matrix.teta.t[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)]<-(1/variance.penalty.t[k+1])*Penalty.matrix.t

#update for penalized Fisher within a penalty.update-loop
I.t.pen<-I.t+Penalty.matrix.teta.t


variance.penalty<-c(variance.penalty.t)


penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))



if (control$print.penalty) {print(penalty);cat("","\n\n")}
if (sum((variance.penalty-variance.penalty.old)^2) < control$tol.penalty) break   #stop for fixpoint iteration


}   #end for (iter.penalty in 1:control$niter.penalty)


} #end if (control$method !="fix") else


#stoping criterion for penalized.backfitting
variance.penalty.list[[iter.epoch+1]]<-variance.penalty
if (sum((variance.penalty.list[[iter.epoch+1]]-variance.penalty.list[[iter.epoch]])^2) < control$tol.epoch.variance & sum((teta.t.new-teta.t.old)^2) < control$tol.epoch.theta) break #stop for outer loop (epoch=update of theta and penalties)




#marginal likelihood
g.inverse.t<-ginverse(Penalty.matrix.teta.t[index.vector.t,index.vector.t],tol=1e-100)
#make.sym.t<-(t(g.inverse.t)+g.inverse.t)/2
#for t-direction
#log.lik.margin.t<-log.lik.t.pen(teta.t)-0.5*(sum(log(eigen(I.t.pen[index.vector.t,index.vector.t])$values)) + sum(log(eigen(make.sym.t)$values)))
log.lik.margin.t<-log.lik.t.pen(teta.t)-0.5*Re(sum(log(eigen(I.t.pen[index.vector.t,index.vector.t]%*%g.inverse.t)$values))) 


if (control$print.log.lik)
  {
cat("","\n\n")    
cat("log.lik.margin.t= ",log.lik.margin.t,"\n")
  }


log.lik.margin<-c(log.lik.margin.t=log.lik.margin.t)
log.lik.margin.t.list[[iter.epoch]]<-log.lik.margin.t
  

} #end for (iter.epoch in 1:control$niter.epoch) 





cat("","\n\n")

#penalties
penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))


#print resulting estimates and penalty values
if (control$print.estimates)
  {
cat("resulting estimates: random effects and penalties","\n\n\n")  #last updates
print(u.t)
cat("","\n\n")
print(penalty)
  }


#likelihood from the 1.iteration
log.lik.margin.start<-c(log.lik.margin.t.list[[1]])
names(log.lik.margin.start)<-c("log.lik.margin.t")


factor.names<-colnames(Design.variables)


rm(offset.list,offset.tij.list,unlist.offset.tij.list,y.poisson.list,y.poisson,I.t,I.t.pen,Design.variables,inverse.t)
for (i in 1:30) gc()




##################################################################
#calculate variances of the resulting estimates for theta#########
##################################################################


#penalized and unpenalized Fiher.type matrix for the full design
I.t.0<-crossprod(Design.matrix.t,Design.matrix.t*exp(as.vector(Design.matrix.t%*%teta.t)+offset.model))
I.t<-I.t.0+Penalty.matrix.teta.t


#estimate covariances through sandwich estimator
inverse.I.t<-qr.solve(I.t,tol=1e-100)
co.variance.teta<-inverse.I.t%*%I.t.0%*%inverse.I.t
variance.teta<-diag(co.variance.teta)
co.variance.teta.t<-co.variance.teta[1:(K.t+p*K.t),1:(K.t+p*K.t)]
names(variance.teta)<-c(names(teta.t))
variance.teta.t<-variance.teta[1:(K.t+p*K.t)]


#degrees.of.freedom for components
#total
df.total<-sum(diag(inverse.I.t%*%I.t.0))
#components specific
df.t.Baseline<-sum(diag((inverse.I.t%*%I.t.0)[1:K.t,1:K.t]))
if (p > 0) for (k in 1:p)
  {
assign(paste("df.t.",factor.names[k+1],sep=""),sum(diag((inverse.I.t%*%I.t.0)[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)])))
  }
#write out d.f.
df<-matrix(0,nrow=p+1,ncol=1)
df[1,]<-c(df.t.Baseline)
if (p > 0) for (k in 1:p) df[k+1,]<-c(get(paste("df.t.",factor.names[k+1],sep="")))
dimnames(df)<-list(c("Baseline",factor.names[-1]),c("t"))




#Grid for survival.time and saison.time
grid.t<-seq(min(data.set$time),max(data.set$time),le=1000)
names(grid.t)<-1:length(grid.t)
#specify knots
if (length(grep("num.knots.t",names(control)))==0) knots.t<-default.knots(grid.t) else knots.t<-default.knots(grid.t,num.knots.t)
#specify spline.bases
B.grid.t<-splineDesign(knots=knots.t,ord=bs.order,x=grid.t)



#Confidence Bands
C.t.grid<-B.grid.t
#for Baseline.t
co.variance.teta.t.Baseline<-co.variance.teta.t[1:K.t,1:K.t]
variance.t.Baseline<-apply(C.t.grid,1,FUN=function(help.row) t(help.row)%*%co.variance.teta.t.Baseline%*%help.row)
deviation.t.Baseline<-qnorm(0.975)*sqrt(variance.t.Baseline)
#for covariables
if (p > 0)
for (k in 1:p)
  {
assign(paste("co.variance.teta.t.",factor.names[k+1],sep=""),co.variance.teta.t[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)])
assign(paste("variance.t.",factor.names[k+1],sep=""),apply(C.t.grid,1,FUN=function(help.row) t(help.row)%*%get(paste("co.variance.teta.t.",factor.names[k+1],sep=""))%*%help.row))
assign(paste("deviation.t.",factor.names[k+1],sep=""),qnorm(0.975)*sqrt(get(paste("variance.t.",factor.names[k+1],sep=""))))
}


#penalty values
penalties<-1/variance.penalty.list[[length(variance.penalty.list)]]
names(penalties)<-sub("variance","penalty",names(penalties))
random.coef<-c(u.t)
var.random<-c(variance.teta.t[grep("u",names(variance.teta.t))])
names(var.random)<-sub("u","variance.u",names(var.random))




#varying coefficients for Baseline and covariables
alpha.t.Baseline<-B.grid.t%*%u.t[grep("Baseline",names(u.t),fixed=TRUE)]
if (p > 0)
{
for (k in 1:p)
  {
assign(paste("alpha.t.",factor.names[k+1],sep=""),B.grid.t%*%u.t[grep(paste("u.t.",factor.names[k+1],sep=""),names(u.t),fixed=TRUE)])  
}
}




#write out varying coefficient, deviation and grid vectors in frames

#varying coefficients
list.t.frame<-list()
t.frame<-NULL;b.frame<-NULL
if (p > 0)
{  
for (k in 1:p)
  {list.t.frame[[k]]<-get(paste("alpha.t.",factor.names[k+1],sep=""))
   t.frame<-cbind(t.frame,list.t.frame[[k]])
  }
varying.frame<-data.frame(alpha.t.Baseline,t.frame)
names(varying.frame)<-c("alpha.t.Baseline",paste("alpha.t.",factor.names[2:(p+1)],sep=""))
} else
{varying.frame<-data.frame(alpha.t.Baseline)
 names(varying.frame)<-c("alpha.t.Baseline")
}


#deviations
list.t.frame<-list()
t.frame<-NULL
if (p > 0)
{  
for (k in 1:p)
  {list.t.frame[[k]]<-get(paste("deviation.t.",factor.names[k+1],sep=""))
   t.frame<-cbind(t.frame,list.t.frame[[k]])
  } 
deviation.frame<-data.frame(deviation.t.Baseline,t.frame)
names(deviation.frame)<-c("deviation.t.Baseline",paste("deviation.t.",factor.names[2:(p+1)],sep=""))
} else
{deviation.frame<-data.frame(deviation.t.Baseline)
names(deviation.frame)<-c("deviation.t.Baseline")
}


#grids
grid.frame<-data.frame(grid.t=grid.t)


#output elements
list(random.coef=random.coef,penalty=penalty,var.random=var.random,log.lik.margin.start=log.lik.margin.start,log.lik.margin=log.lik.margin,df=df,df.total=sum(df),niter.epoch=iter.epoch,varying.frame=varying.frame,deviation.frame=deviation.frame,grid.frame=grid.frame,p=p,factor.names=factor.names[-1])



}#end of function

 
