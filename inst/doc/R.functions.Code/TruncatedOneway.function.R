TruncatedOneway<-function(data.set=data.set,control=control)  
{  

#############################################################
##One Way Hazard Model With Varying Coefficients#############
#############################################################

#############################################################  
#modelling with truncated polynomilas########################
#############################################################  
  
if (names(data.set)[3]=="Intercept") data.set<-data.set[,-3]

N<-length(data.set[,1])   #total number of observations in data.set




############################################################################################
#knots selection (s. Ngo, Wand: Smoothing with mixed model software, 2003, Part 3)##########
############################################################################################
default.knots<-function(x,num.knots)
{
   if (missing(num.knots))
      num.knots <- max(5,min(floor(length(unique(x))/4),35))
   return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))]))
}


if (length(grep("num.knots.t",names(control)))==0)  knots.t<-default.knots(data.set$time[data.set$status==1]) else knots.t<-default.knots(data.set$time[data.set$status==1],as.vector(unlist(control[grep("num.knots.t",names(control))])))
K.t<-length(knots.t)




#################################################
#create artificial poisson data##################
#################################################
survival.to.poisson <- function(time=time, status=status, x=NULL)
{
N <- length(time)
event.time <- time[status==1]

#create the grid
#grid<-c(unique(sort(event.time)),max(time+1))
#grid<-c(sort(sample(unique(time[status==1]),control$number.int)),max(time+1))
grid<-round(quantile(unique(time),prob=seq(0,1,le=control$number.int)),max(time+1))

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
}

#create a model.matrix for covariables
if (dim(data.set)[2] > 2) model.matrix.x<-model.matrix(as.formula(paste("~",paste(names(data.set)[3:length(names(data.set))],collapse="+"))),data=data.set) else model.matrix.x<-model.matrix(~-1+rep(1,N))


help.surv<-survival.to.poisson(data.set$time,data.set$status,model.matrix.x)  
grid.points<-help.surv$grid   #t-points, at which the integral is being approximated
#for each observation there is an element of the list with accordingly:
y.poisson.list<-help.surv$y.list  #poisson.data
y.poisson<-unlist(y.poisson.list)
offset.list<-help.surv$o.list     #offset.parameter
m.list<-help.surv$m.list          #number of poisson data


Design.variables<-help.surv$x
if (dim(data.set)[2] > 2) colnames(Design.variables)<-c("Intercept",colnames(model.matrix.x)[2:length(colnames(model.matrix.x))]) else colnames(Design.variables)<-c("Intercept")


#number of varaibles (including factor levels by factor variables) in the defined design matrix
p<-ncol(Design.variables)-1 


#define survival time and birth time for expanded data set, and the (truncated polynomial) spline bases
Z.time.list<-sapply(data.set$time,FUN=function(x) c(grid.points[x>grid.points],x))
Z.time<-unlist(Z.time.list)  #Z.time=t




#define the start values for beta parameter by means of glm regression for y.poisson#########
if (p > 0)
  {
variables.time<-Design.variables[,2:ncol(Design.variables)]*Z.time
variables.time<-as.matrix(variables.time)
colnames(variables.time)<- paste(colnames(Design.variables)[2:ncol(Design.variables)],".t",sep="")
rownames(variables.time)<-NULL

beta.start<-coef(glm(as.formula(paste("y.poisson~ ",paste(c("Z.time",colnames(Design.variables)[2:ncol(Design.variables)],colnames(variables.time),"offset(unlist(offset.list))"),collapse="+"))),data=data.frame(Z.time,Design.variables,variables.time),family=poisson))
names(beta.start)<-c("beta.t.Baseline.intercept","beta.t.Baseline.slope",paste("beta.t.",colnames(Design.variables)[2:ncol(Design.variables)],".intercept",sep=""),paste("beta.t.",colnames(Design.variables)[2:ncol(Design.variables)],".slope",sep=""))
} else 
{
beta.start<-coef(glm(as.formula(paste("y.poisson~ ",paste(c("Z.time","offset(unlist(offset.list))"),collapse="+"))),data=data.frame(Z.time),family=poisson))
names(beta.start)<-c("beta.t.Baseline.intercept","beta.t.Baseline.slope") 
}




#coefficients at Baseline                      
beta.t.Baseline.intercept<-beta.start[1]    #for Intercept
beta.t.Baseline.slope<- beta.start[2]  #for t-Trend
names(beta.t.Baseline.intercept)<-"beta.t.Baseline.intercept"
names(beta.t.Baseline.slope)<-"beta.t.Baseline.slope"

#coefficients at covariables                      
beta.t.intercept.list<-list()
beta.t.slope.list<-list()                   
beta.t.intercept.list[[1]]<-beta.t.Baseline.intercept
beta.t.slope.list[[1]]<-beta.t.Baseline.slope
if (p > 0)
  {
for (i in 3:(3+p-1)) beta.t.intercept.list[[i-1]]<-beta.start[i]
for (i in (3+p-1+1):(3+2*p-1)) beta.t.slope.list[[i-(3+p-1-1)]]<-beta.start[i]
}

beta.t.intercept<-unlist(beta.t.intercept.list)
beta.t.slope<-unlist(beta.t.slope.list)
                     



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



#create penalty matrix
Lambda.t<-eval(parse(text=paste("diag(c(",paste(paste("0,0,rep(1/variance.penalty.t[",1:length(variance.penalty.t),sep=""),"],K.t)",sep="",collapse=","),"))",sep="")))
D.t <-diag(rep(c(0,0,rep(1,K.t)),p+1))


#initialize random effects
u.t.list<-list()
if (p > 0)
  {
for (i in 1:(p+1)) u.t.list[[i]]<-rep(0,K.t)
names(u.t.list)<-c("u.t.Baseline.",paste("u.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))                         
} else
  {
u.t.list[[1]]<-rep(0,K.t);names(u.t.list)<-c("u.t.Baseline.")
}

u.t<-unlist(u.t.list)                                                 
                               


#combine in theta vector
teta.t.Baseline<-c(beta.t.intercept.list[[1]],beta.t.slope.list[[1]],u.t[1:length(u.t.list[[1]])])
teta.t.list<-list()
teta.t.list[[1]]<-teta.t.Baseline
if (p > 0) for (i in 1:p) teta.t.list[[i+1]]<-c(beta.t.intercept.list[[i+1]],beta.t.slope.list[[i+1]],u.t[(i*length(u.t.list[[i]])+1):((i+1)*length(u.t.list[[i]]))]) 
teta.t<-unlist(teta.t.list)





#construct design matrix for the model
Basis.t<-outer(Z.time,knots.t,FUN="-");Basis.t<-Basis.t*(Basis.t>0)
if (p > 0)
{
variables.t<-Design.variables[,2:ncol(Design.variables),drop=FALSE]
Design.matrix.t<-cbind(1,Z.time,Basis.t,eval(parse(text=paste("cbind(",paste("variables.t[,",1:p,"]*cbind(1,Z.time,Basis.t)",sep="",collapse=","),")"))))
} else
{
Design.matrix.t<-cbind(1,Z.time,Basis.t)
}


if (control$print.epoch)
  {
cat("start values: fixed parameters of the varying coefficients","\n\n\n")
print(beta.t.intercept)
print(beta.t.slope)
cat("","\n\n")

penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))
print(penalty)
cat("","\n\n")
}




#initialize list object for penalty variances from following iterations
variance.penalty.list<-list()
variance.penalty.list[[1]]<-variance.penalty   #first element of the list is the initial variance vector



status<-data.set$status


#saves the values of likelihood for all iterations
log.lik.margin.t.list<-list()  



###############################
#optimization routine##########
###############################
for (iter.epoch in 1:control$niter.epoch)
 {
if (control$print.epoch) cat("epoch= ", iter.epoch, "\n\n")



#likelihood for alpha.t
offset.tij.list<-lapply(1:N,FUN=function(i) offset.list[[i]])


#score for theta.t
S.t<-crossprod(Design.matrix.t,as.vector(unlist(y.poisson.list)-exp(Design.matrix.t%*%teta.t+unlist(offset.tij.list))))
S.t.pen<-S.t-Lambda.t%*%D.t%*%teta.t
#Fisher.type.Matrix(observed)
I.t<-crossprod(Design.matrix.t,Design.matrix.t*as.vector(exp(Design.matrix.t%*%teta.t+unlist(offset.tij.list))))
I.t.pen<-I.t+Lambda.t%*%D.t 
#update of theta.t
inverse.t <- qr.solve(I.t.pen,tol=1e-100)
teta.t.old<-teta.t
teta.t.new<-teta.t.old+inverse.t%*%S.t.pen
names(teta.t.new)<-names(teta.t)    
teta.t<-teta.t.new
names(teta.t)<-names(teta.t.new)    




#update of beta's:
beta.t<-teta.t.new[grep("beta.t",names(teta.t),fixed=TRUE)]
#update of random effects
u.t<-teta.t.new[grep("u.t",names(teta.t),fixed=TRUE)]


if (control$print.epoch)
{                         
print(beta.t)
cat("","\n\n")
}



###################################################################
#optimization procedure for penalties##############################
###################################################################


#define first the submatrices of the (penalized) fisher.information, correspondig to the random coefficients
index.vector.t<-(2+1):(2+K.t)  #for u.t.Baseline
if (p > 0) for (k in 1:p) index.vector.t<-c(index.vector.t,(2+K.t+2*k+1+(k-1)*K.t):(2+K.t+2*k+k*K.t))



#chose optimization routine for smoothing.parameters; one of both fix.point or Newton.Raphson
if (control$method=="fix")
  {
    
####################
#fix.point iteration
####################
for (iter.penalty in 1:control$niter.penalty)
  {  
if (control$print.penalty) cat("iteration.penalty= ",iter.penalty,"\n\n")  


I.t.u.t<-I.t.pen[index.vector.t,index.vector.t]    #penalizied fisher.type matrix for t
inverse.I.t.u.t<-qr.solve(I.t.u.t,tol=1e-100)


variance.penalty.t.old<-variance.penalty.t
variance.penalty.old<-c(variance.penalty.t.old)   


variance.t.Baseline<-(sum(u.t[grep("Baseline",names(u.t),fixed=TRUE)]^2) + sum(diag(inverse.I.t.u.t[1:K.t,1:K.t])))/K.t
variance.penalty.t[1]<-variance.t.Baseline
if ( p > 0) 
for (k in 1:p) variance.penalty.t[k+1]<-(sum(u.t[(k*K.t+1):((k+1)*K.t)]^2) + sum(diag(inverse.I.t.u.t[(k*K.t+1):((k+1)*K.t),(k*K.t+1):((k+1)*K.t)])))/K.t

Lambda.t<-eval(parse(text=paste("diag(c(",paste(paste("0,0,rep(1/variance.penalty.t[",1:length(variance.penalty.t),sep=""),"],K.t)",sep="",collapse=","),"))",sep="")))
I.t.pen<-I.t+Lambda.t%*%D.t 



variance.penalty<-c(variance.penalty.t)


penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))



if (control$print.penalty) {print(penalty);cat("","\n\n")}
if (sum((variance.penalty-variance.penalty.old)^2) < control$tol.penalty) break   #stop for fixpoint iteration


}   #end for (iter.penalty in 1:control$niter.penalty)


} else    #one of both via fix.point or Newton.Raphson


{
#########################
#Newton.Raphson iteration
#########################
I.t.u.t<-I.t[index.vector.t,index.vector.t]  #not penalizied fisher.type matrix for t

for (iter.penalty in 1:control$niter.penalty)
  {  
if (control$print.penalty) cat("iteration.penalty= ",iter.penalty,"\n\n")  


variance.penalty.t.old<-variance.penalty.t
variance.penalty.old<-c(variance.penalty.t.old)   

lambda.t<-1/variance.penalty.t



##################################
#I.derivative w.r.t lambda.t######
##################################
D.lambda.t<-diag(1,(p+1)*K.t)*rep(lambda.t,rep(K.t,p+1))

deriv.D.lambda.t<-list()
deriv.D.lambda.t[[1]]<-diag(c(rep(1,K.t),rep(0,p*K.t)))  #for Baseline
if (p > 0)
for (k in 1:p) deriv.D.lambda.t[[k+1]]<-diag(c(rep(0,K.t),rep(0,(k-1)*K.t),rep(1,K.t),rep(0,(p-k)*K.t)))  #for covariables

#Score for log.lik.margin.t  
S.lambda.t<-numeric(p+1)
#for Baseline
S.lambda.t[1]<- -0.5*sum(u.t[grep("Baseline",names(u.t),fixed=TRUE)]^2)+0.5*K.t/lambda.t[1]-0.5*sum(diag(qr.solve(I.t.u.t+D.lambda.t,tol=1e-100)%*%deriv.D.lambda.t[[1]]))
#for covariables
if (p > 0)
for (k in 1:p) S.lambda.t[k+1]<- -0.5*sum(u.t[(k*K.t+1):((k+1)*K.t)]^2)+0.5*K.t/lambda.t[k+1]-0.5*sum(diag(qr.solve(I.t.u.t+D.lambda.t,tol=1e-100)%*%deriv.D.lambda.t[[k+1]]))


##################################
#II.derivative w.r.t lambda.t#####
##################################
help.inverse<-qr.solve(I.t.u.t+D.lambda.t,tol=1e-100)
help.apply<-apply(expand.grid(1:(p+1),1:(p+1)),1,FUN=function(index) -0.5*sum(diag(help.inverse%*%deriv.D.lambda.t[[index[1]]]%*%help.inverse%*%deriv.D.lambda.t[[index[2]]])))
I.lambda.t<- -(matrix(help.apply,nrow=p+1,ncol=p+1,byrow=FALSE)-0.5*diag(K.t/as.vector(lambda.t)^2,length(lambda.t)))   #-II.derivative

##################################################
#Newton.Raphson iteration (update of lambda.t)####
##################################################
lambda.t.old<-lambda.t
lambda.t.new<-as.vector(lambda.t.old+qr.solve(I.lambda.t,tol=1e-100)%*%S.lambda.t)    
lambda.t<-lambda.t.new


variance.penalty.t<-1/lambda.t
names(variance.penalty.t)<-names(variance.penalty.t.list)

Lambda.t<-eval(parse(text=paste("diag(c(",paste(paste("0,0,rep(1/variance.penalty.t[",1:length(variance.penalty.t),sep=""),"],K.t)",sep="",collapse=","),"))",sep="")))



variance.penalty<-c(variance.penalty.t)


penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))


if (control$print.penalty) {print(penalty);cat("","\n\n")}
if (sum((variance.penalty-variance.penalty.old)^2) < control$tol.penalty) break   #stop for Newton.Raphson iteration




}   #end for (iter.penalty in 1:control$niter.penalty)


}  #end one of both: fix.point or Netwon.Raphson



#stoping criterion for penalized.backfitting
variance.penalty.list[[iter.epoch+1]]<-variance.penalty
if (sum((variance.penalty.list[[iter.epoch+1]]-variance.penalty.list[[iter.epoch]])^2) < control$tol.epoch.variance & sum((teta.t.new-teta.t.old)^2) < control$tol.epoch.theta) break    #stop for outer loop (epoch=update of theta and penalties)



#marginal likelihood
log.lik.pen.t<-sum(as.vector(unlist(y.poisson.list))*(Design.matrix.t%*%teta.t)-exp(Design.matrix.t%*%teta.t+unlist(offset.tij.list)))-0.5*t(teta.t)%*%Lambda.t%*%D.t%*%teta.t 
log.lik.margin.t<-log.lik.pen.t-0.5*sum(log(eigen(I.t.pen[index.vector.t,index.vector.t]%*%qr.solve(Lambda.t[index.vector.t,index.vector.t],tol=1e-100))$values))



if (control$print.log.lik)
  {
cat("","\n\n")    
cat("log.lik.margin.t= ",log.lik.margin.t,"\n")
  }



log.lik.margin<-c(log.lik.margin.t=log.lik.margin.t)

log.lik.margin.t.list[[iter.epoch]]<-log.lik.margin.t


}  #end (for iter.epoch in 1:control$niter.epoch)



cat("","\n\n\n")

penalty.t<-1/variance.penalty.t
penalty<-c(penalty.t)
names(penalty)<-sub("variance","penalty",names(variance.penalty))


if (control$print.estimates)
  {
cat("resulting estimates: beta and penalties","\n\n\n")  #last updates
print(beta.t)
cat("","\n\n")
print(penalty)
  }


#likelihood from the 1.iteration
log.lik.margin.start<-c(log.lik.margin.t.list[[1]])
names(log.lik.margin.start)<-c("log.lik.margin.t")



factor.names<-colnames(Design.variables)

rm(offset.tij.list,y.poisson.list,y.poisson,I.t,I.t.pen,Design.variables,Z.time,inverse.t,m.list)
for (i in 1:30) gc()




##################################################################
##################################################################
#calculate variances of the resulting estimates for theta#########
##################################################################
##################################################################

#fisher.type.matrix for (complete) model
I.t.0<-crossprod(Design.matrix.t,Design.matrix.t*exp(as.vector(Design.matrix.t%*%teta.t)+unlist(offset.list)))
I.t<-I.t.0+diag(diag(Lambda.t))


#estimate covariances through sandwich estimator
inverse.I.t<-qr.solve(I.t,tol=1e-100)
co.variance.teta<-inverse.I.t%*%I.t.0%*%inverse.I.t
variance.teta<-diag(co.variance.teta)
co.variance.teta.t<-co.variance.teta[1:(2+K.t+p+p+p*K.t),1:(2+K.t+p+p+p*K.t)]
names(variance.teta)<-c(names(teta.t))
variance.teta.t<-variance.teta[1:(2+K.t+p+p+p*K.t)]


##########################################################
#define degrees.of.freedom for components#################
##########################################################
#total
#df.total<-sum(diag(inverse.I.t%*%I.t.0))
#components specific
df.t.Baseline<-sum(diag((inverse.I.t%*%I.t.0)[1:(2+K.t),1:(2+K.t)]))
if (p > 0) for (k in 1:p)
  {
assign(paste("df.t.",factor.names[k+1],sep=""),sum(diag((inverse.I.t%*%I.t.0)[(2+K.t+2*(k-1)+(k-1)*K.t+1):(2+K.t+2*k+k*K.t),(2+K.t+2*(k-1)+(k-1)*K.t+1):(2+K.t+2*k+k*K.t)])))
  }
#write out d.f.
df<-matrix(0,nrow=p+1,ncol=1)
df[1,]<-c(df.t.Baseline)
if (p > 0) for (k in 1:p) df[k+1,]<-c(get(paste("df.t.",factor.names[k+1],sep="")))
dimnames(df)<-list(c("Baseline",factor.names[-1]),c("t"))





#Grid for time and birth
grid.t<-seq(min(data.set$time),max(data.set$time),le=1000)
B.grid.t<-outer(grid.t,knots.t,FUN="-")
B.grid.t<-B.grid.t*(B.grid.t>0)



#############################
#Confidence Bands############
#############################
C.t.grid<-cbind(1,grid.t,B.grid.t)
#for Baseline.t
co.variance.teta.t.Baseline<-co.variance.teta.t[1:(2+K.t),1:(2+K.t)]
variance.t.Baseline<-apply(C.t.grid,1,FUN=function(help.row) t(help.row)%*%co.variance.teta.t.Baseline%*%help.row)  
deviation.t.Baseline<-qnorm(0.975)*sqrt(variance.t.Baseline)
#for covariables
if (p > 0)
for (k in 1:p)
  {
assign(paste("co.variance.teta.t.",factor.names[k+1],sep=""),co.variance.teta.t[(2+K.t+2*k-1+(k-1)*K.t):(2+K.t+2*k+k*K.t),(2+K.t+2*k-1+(k-1)*K.t):(2+K.t+2*k+k*K.t)])
assign(paste("variance.t.",factor.names[k+1],sep=""),apply(C.t.grid,1,FUN=function(help.row) t(help.row)%*%get(paste("co.variance.teta.t.",factor.names[k+1],sep=""))%*%help.row))
assign(paste("deviation.t.",factor.names[k+1],sep=""),qnorm(0.975)*sqrt(get(paste("variance.t.",factor.names[k+1],sep=""))))
}


fix.coef<-c(beta.t)
penalties<-1/variance.penalty.list[[length(variance.penalty.list)]]
names(penalties)<-sub("variance","penalty",names(penalties))
random.coef<-c(u.t)
var.fix<-c(variance.teta.t[grep("beta",names(variance.teta.t),fixed=TRUE)])
names(var.fix)<-sub("beta","variance.beta",names(var.fix))
var.random<-c(variance.teta.t[grep("u",names(variance.teta.t),fixed=TRUE)])
names(var.random)<-sub("u","variance.u",names(var.random))




#varying coefficients for Baseline and covariables
alpha.t.Baseline <-cbind(1,grid.t)%*%beta.t[grep("Baseline",names(beta.t),fixed=TRUE)]+B.grid.t%*%u.t[grep("Baseline",names(u.t),fixed=TRUE)] 
if (p > 0)
{
for (k in 1:p)
  {
assign(paste("alpha.t.",factor.names[k+1],sep=""),cbind(1,grid.t)%*%beta.t[grep(paste("beta.t.",factor.names[k+1],sep=""),names(beta.t),fixed=TRUE)]+B.grid.t%*%u.t[grep(paste("u.t.",factor.names[k+1],sep=""),names(u.t),fixed=TRUE)])
  }
}



#write out varying coefficient, deviation and grid vectors in frames


#varying coefficients
list.t.frame<-list()
t.frame<-NULL
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



list(fix.coef=fix.coef,random.coef=random.coef,penalty=penalty,var.fix=var.fix,var.random=var.random,log.lik.margin.start=log.lik.margin.start,log.lik.margin=log.lik.margin,df=df,df.total=sum(df),niter.epoch=iter.epoch,varying.frame=varying.frame,deviation.frame=deviation.frame,grid.frame=grid.frame,p=p,factor.names=factor.names[-1])



}  #end of function
