BackfitBspline<-function(data.set=data.set,control=control)
{

#############################################################
##impelmentation accordingly to: Kauermann, Khomski:#########
##Additiv Two Way Hazard Model With Variing Coefficients#####
##with periodic 2nd time-scale component#####################
#############################################################

#############################################################  
#modelling with B-splines (all coefficients are penalized)###
#############################################################  

library(splines)

if (names(data.set)[4]=="Intercept") data.set<-data.set[,-4]

N<-length(data.set[,1])   #total number of observations in data.set



#################################################
#knots and order of B-spline Bases###############
#################################################
default.knots<-function(x,num.knots)
{
bs.order<-4   
bs.degree<-bs.order-1   

if (missing(num.knots)) num.knots<-20  #10

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




ginverse<-function(X,tol=1e-100)
{
Xsvd<-svd(X,LINPACK=TRUE)
if (is.complex(X)) Xsvd$u<-Conj(Xsvd$u)
Positive<-Xsvd$d > max(tol * Xsvd$d[1], 0)
if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
else if (!any(Positive)) array(0, dim(X)[2:1])
else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}








##############################################################################################################
#create a model.matrix for (co)variables######################################################################
#N.B.!!! the  reference categories for factors must be defined in the analizied data set bevore applied hier##
##############################################################################################################
if (dim(data.set)[2] > 3) model.matrix.x<-model.matrix(as.formula(paste("~",paste(names(data.set)[4:length(names(data.set))],collapse="+"))),data=data.set) else model.matrix.x<-model.matrix(~-1+rep(1,N))



help.surv<-survival.to.poisson(data.set$time,data.set$status,model.matrix.x)  
grid.points<-help.surv$grid   #t-points, at which the integral is being approximated
#for each observation there is an element of the list with accordingly:
y.poisson.list<-help.surv$y.list  #poisson.data
y.poisson<-unlist(y.poisson.list)
offset.list<-help.surv$o.list     #offset.parameter
offset.model<-unlist(offset.list)
m.list<-help.surv$m.list          #number of poisson data




Design.variables<-help.surv$x
if (dim(data.set)[2] > 3) colnames(Design.variables)<-c("Intercept",colnames(model.matrix.x)[2:length(colnames(model.matrix.x))]) else colnames(Design.variables)<-c("Intercept")

p<-ncol(Design.variables)-1 #number of varaibles (including factor levels by factor variables) in the defined design matrix



#########################################
#knots and order of B-splines############
#########################################
bs.order<-4   
bs.degree<-bs.order-1 

if (length(grep("num.knots.t",names(control)))==0) knots.t<-default.knots(data.set$time) else knots.t<-default.knots(data.set$time,num.knots)
if (length(grep("num.knots.b",names(control)))==0) knots.b<-default.knots(data.set$birth) else knots.b<-default.knots(data.set$birth,num.knots)



##############################################################################################################
#define survival time and saison.component time for expanded data set, and the b-spline bases#################
##############################################################################################################
Z.time.list<-sapply(data.set$time,FUN=function(x) c(grid.points[grid.points<x],x))
Z.time<-unlist(Z.time.list)
Basis.t<-splineDesign(knots=knots.t,ord=bs.order,x=Z.time)
K.t<-dim(Basis.t)[2]
Z.birth.list<-sapply(seq(data.set$birth),FUN=function(i) rep(data.set$birth[i],length(Z.time.list[[i]])))
Z.birth<-unlist(Z.birth.list)
names(Z.birth)<-1:length(Z.birth)
Basis.b<-splineDesign(knots=knots.b,ord=bs.order,x=Z.birth)
####################################################
#periodic continuation (wrapping) in b-direction####
####################################################
Basis.b[,1:bs.degree]<-Basis.b[,1:bs.degree]+Basis.b[,(ncol(Basis.b)-(bs.degree-1)):ncol(Basis.b)]
Basis.b<-Basis.b[,1:(ncol(Basis.b)-bs.degree)]

K.b<-dim(Basis.b)[2]


############################################################################################
#N.B.: at the start point "0"(=begin of the period) set varying.coeff in b-direction=0######
############################################################################################
#transform the wrapped B-spline matrix in b-direction appropriately#########################
############################################################################################
c.vector<- -Basis.b[as.numeric(names(Z.birth)[Z.birth==min(Z.birth)][1]),-1]/unique(Basis.b[Z.birth==min(Z.birth),1])

Basis.b<-outer(Basis.b[,1],c.vector)+Basis.b[,-1]




########################################
#build Difference matrices##############
########################################

#for t-direction
################
Diff.2.t<-matrix(0,nrow=K.t-2,ncol=K.t)    #diff.penalty of 2.order
for (i in 1:(K.t-2)) Diff.2.t[i,i:(i+2)]<-c(1,-2,1)
Penalty.matrix.t<-crossprod(Diff.2.t)

#for b-direction: cyclic(=wrapped)
##################################
Diff.2.b<-matrix(0,nrow=K.b-2+bs.degree,ncol=K.b+bs.degree)
for (i in 1:(K.b+bs.degree-2)) Diff.2.b[i,i:(i+2)]<-c(1,-2,1)
Diff.2.b[,1:bs.degree]<-Diff.2.b[,1:bs.degree]+Diff.2.b[,(ncol(Diff.2.b)-(bs.degree-1)):ncol(Diff.2.b)]
Diff.2.b<-Diff.2.b[,1:(ncol(Diff.2.b)-bs.degree)]
Penalty.matrix.b<-crossprod(Diff.2.b)

#transform the penalty matrix appropriately
Penalty.matrix.b<-outer(c.vector*Penalty.matrix.b[1,1],c.vector)+outer(c.vector,Penalty.matrix.b[1,-1])+outer(Penalty.matrix.b[-1,1],c.vector)+Penalty.matrix.b[-1,-1]

#new setting becouse of transformation
K.b<-dim(Basis.b)[2]   



#remove not more used objects
rm(Z.time.list,Z.birth.list)
for (i in 1:30) gc()



############################################################################
#initialize penalty parameters through variances of random effects##########
############################################################################
variance.penalty.t.list<-list()
variance.penalty.b.list<-list()

if (length(grep("start.penalty.t",names(control)))==0 | length(grep("start.penalty.b",names(control)))==0)
{
if (p > 0) 
  {
for (i in 1:(p+1)) variance.penalty.t.list[[i]]<-10e-1
names(variance.penalty.t.list)<-c("variance.t.Baseline",paste("variance.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
for (i in 1:(p+1)) variance.penalty.b.list[[i]]<-10e-1
names(variance.penalty.b.list)<-c("variance.b.Baseline",paste("variance.b.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
  } else
  {
variance.penalty.t.list[[1]]<-10e-1
names(variance.penalty.t.list)<-c("variance.t.Baseline")
variance.penalty.b.list[[1]]<-10e-1
names(variance.penalty.b.list)<-c("variance.b.Baseline")
  }

} else

{
if (p > 0) 
  {
for (i in 1:(p+1)) variance.penalty.t.list[[i]]<-1/control$start.penalty.t
names(variance.penalty.t.list)<-c("variance.t.Baseline",paste("variance.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
for (i in 1:(p+1)) variance.penalty.b.list[[i]]<-1/control$start.penalty.b
names(variance.penalty.b.list)<-c("variance.b.Baseline",paste("variance.b.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
  } else
  {
variance.penalty.t.list[[1]]<-1/control$start.penalty.t
names(variance.penalty.t.list)<-c("variance.t.Baseline")
variance.penalty.b.list[[1]]<-1/control$start.penalty.b
names(variance.penalty.b.list)<-c("variance.b.Baseline")
  }
}

variance.penalty.t<-unlist(variance.penalty.t.list)
variance.penalty.b<-unlist(variance.penalty.b.list)                               
variance.penalty<-c(variance.penalty.t,variance.penalty.b)



############################################
#initialize random effects##################
############################################
u.t.list<-list()
u.b.list<-list()
if (p > 0)
  {
for (i in 1:(p+1)) u.t.list[[i]]<-rep(1,K.t)
names(u.t.list)<-c("u.t.Baseline.",paste("u.t.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))                         
for (i in 1:(p+1)) u.b.list[[i]]<-rep(1,K.b)
names(u.b.list)<-c("u.b.Baseline",paste("u.b.",colnames(Design.variables)[2:ncol(Design.variables)],sep=""))
} else
  {
u.t.list[[1]]<-rep(1,K.t);names(u.t.list)<-c("u.t.Baseline.")
u.b.list[[1]]<-rep(1,K.b);names(u.b.list)<-c("u.b.Baseline.")
}

u.t<-unlist(u.t.list)
u.b<-unlist(u.b.list)   




####################################
#combine in theta vector############
####################################
teta.t.Baseline<-u.t[1:length(u.t.list[[1]])]
teta.t.list<-list()
teta.t.list[[1]]<-teta.t.Baseline
if (p > 0) for (i in 1:p) teta.t.list[[i+1]]<-u.t[(i*length(u.t.list[[i]])+1):((i+1)*length(u.t.list[[i]]))]
teta.t<-unlist(teta.t.list)

teta.b.Baseline<-c(u.b[1:length(u.b.list[[1]])])
teta.b.list<-list()
teta.b.list[[1]]<-teta.b.Baseline  
if (p > 0) for (i in 1:p) teta.b.list[[i+1]]<-u.b[(i*length(u.b.list[[i]])+1):((i+1)*length(u.b.list[[i]]))] 
teta.b<-unlist(teta.b.list)
                      

teta<-c(teta.t,teta.b)



###############################################
#create penalty matrices#######################
###############################################
#for t-
Penalty.matrix.teta.t<-matrix(0,nrow=length(teta.t),ncol=length(teta.t))
Penalty.matrix.teta.t[1:K.t,1:K.t]<-(1/variance.penalty.t[1])*Penalty.matrix.t
if (p > 0)
for (k in 1:p) Penalty.matrix.teta.t[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)]<-(1/variance.penalty.t[k+1])*Penalty.matrix.t

#for b-
Penalty.matrix.teta.b<-matrix(0,nrow=length(teta.b),ncol=length(teta.b))
Penalty.matrix.teta.b[1:K.b,1:K.b]<-(1/variance.penalty.b[1])*Penalty.matrix.b
if (p > 0)
for (k in 1:p) Penalty.matrix.teta.b[(K.b+(k-1)*K.b+1):(K.b+k*K.b),(K.b+(k-1)*K.b+1):(K.b+k*K.b)]<-(1/variance.penalty.b[k+1])*Penalty.matrix.b


Penalty.matrix.tb<-matrix(0,nrow=length(teta.t)+length(teta.b),ncol=length(teta.t)+length(teta.b))
Penalty.matrix.tb[1:length(teta.t),1:length(teta.t)]<-Penalty.matrix.teta.t
Penalty.matrix.tb[(length(teta.t)+1):(length(teta.t)+length(teta.b)),(length(teta.t)+1):(length(teta.t)+length(teta.b))]<-Penalty.matrix.teta.b
                      


#remove not more used objects
rm(help.surv,Basis.b)
for (i in 1:30) gc()



########################################################################
########################################################################
#construct design matrices for the model################################
######################################################################## 
########################################################################
if (p > 0)
{  
variables.t<-Design.variables[,2:ncol(Design.variables),drop=FALSE]
Design.matrix.t<-cbind(Basis.t,eval(parse(text=paste("cbind(",paste("variables.t[,",1:p,"]*Basis.t",sep="",collapse=","),")"))))
} else
Design.matrix.t<-Basis.t

#need renewed calculation of the basis
Z.birth<-data.set$birth
Basis.b<-splineDesign(knots=knots.b,ord=bs.order,x=Z.birth)
Basis.b[,1:bs.degree]<-Basis.b[,1:bs.degree]+Basis.b[,(ncol(Basis.b)-(bs.degree-1)):ncol(Basis.b)]
Basis.b<-Basis.b[,1:(ncol(Basis.b)-bs.degree)]

#transform Basis.b
Basis.b<-outer(Basis.b[,1],c.vector)+Basis.b[,-1]

if (p > 0)
{  
variables.b<-as.matrix(model.matrix.x[,2:ncol(model.matrix.x)])
Design.matrix.b<-cbind(Basis.b,eval(parse(text=paste("cbind(",paste("variables.b[,",1:p,"]*Basis.b",sep="",collapse=","),")"))))
} else
Design.matrix.b<-Basis.b



#remove not more used objects
rm(Basis.b,Basis.t,model.matrix.x)
for (i in 1:30) gc()



###############################################################################
#initialize list object for penalty variances from following iterations########
###############################################################################
variance.penalty.list<-list()  
variance.penalty.list[[1]]<-variance.penalty   #first element of the list is the initial variance vector



status<-data.set$status


log.lik.margin.t.list<-list()   #saves the values of marginal likelihood for all iterations
log.lik.margin.b.list<-list()
#log.lik.margin.tb.list<-list()


for (iter.epoch in 1:control$niter.epoch)
 {
if (control$print.epoch) cat("epoch= ", iter.epoch, "\n\n")


###############################################
#likelihood for alpha.b, alpha.t be given######
##############################################
help.cut<-cumsum(unlist(m.list))   
help.matrix<-cbind(c(0,help.cut[-length(help.cut)])+1,help.cut)  
indizes<-lapply(1:N,FUN=function(i) help.matrix[i,1]:help.matrix[i,2])  
offset.bi<-unlist(lapply(indizes, FUN=function(ind)  log(sum(exp(Design.matrix.t[ind,]%*%teta.t+offset.model[ind])))))


##############################
#log.likelihood for b#########
##############################
log.lik.b.pen<-function(x)
{
log.lik.b<-sum(status*Design.matrix.b%*%x-exp(Design.matrix.b%*%x+offset.bi))
log.lik.b-0.5*t(x)%*%Penalty.matrix.teta.b%*%x
}

#score for theta.b
S.b<-crossprod(Design.matrix.b,status-exp(Design.matrix.b%*%teta.b+offset.bi))
S.b.pen<-S.b-Penalty.matrix.teta.b%*%teta.b
#Fisher.type.Matrix(observed)
I.b<-crossprod(Design.matrix.b,Design.matrix.b*as.vector(exp(Design.matrix.b%*%teta.b+offset.bi)))
I.b.pen<-I.b+Penalty.matrix.teta.b
#add a small postive value of the identity matrix (for the case not being invertaible)
delta<-1e-05
tau<-max(0,delta-min(Re(eigen(I.b.pen)$values)))
I.b.pen<-I.b.pen+tau*diag(1,ncol(I.b.pen))
inverse.b <- qr.solve(I.b.pen,tol=1e-100)

teta.b.old<-teta.b
delta.teta.b<-inverse.b%*%S.b.pen


###########################################################
#optimization goes by backtracking.line.search procedure###
###########################################################
#fix parameters
alpha.backtracking<-0.01
beta.backtracking<-0.5
#optimal value of the step s
s<-1
while (-log.lik.b.pen(teta.b.old+s*delta.teta.b) > -log.lik.b.pen(teta.b.old)+alpha.backtracking*s*crossprod(-S.b.pen,delta.teta.b)) s<-s*beta.backtracking

#update of theta.b
teta.b.new<-teta.b.old+s*delta.teta.b
names(teta.b.new)<-names(teta.b)    
teta.b<-teta.b.new
names(teta.b)<-names(teta.b.new)     

#discent.vector.b<-crossprod(-S.b.pen,inverse.b%*%S.b.pen)



##############################################    
#likelihood for alpha.t, alpha.b be given#####
##############################################
offset.tij.list<-lapply(1:N,FUN=function(i) matrix(rep(Design.matrix.b[i,],m.list[i]),nrow=m.list[i],byrow=TRUE)%*%teta.b + offset.list[[i]])
unlist.offset.tij.list<-unlist(offset.tij.list)


##############################
#log.likelihood for t#########
##############################
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
#add a small postive value of the identity matrix (for the case not being invertaible)
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
while (-log.lik.t.pen(teta.t.old+s*delta.teta.t) > -log.lik.t.pen(teta.t.old)+alpha.backtracking*s*crossprod(-S.t.pen,delta.teta.t)) s<-s*beta.backtracking

teta.t.new<-teta.t.old+s*delta.teta.t
names(teta.t.new)<-names(teta.t)    
teta.t<-teta.t.new
names(teta.t)<-names(teta.t.new)    

#discent.vector.t<-crossprod(-S.t.pen,inverse.t%*%S.t.pen)


################################
#update of random effects#######
################################
u.t<-teta.t.new[grep("u.t",names(teta.t),fixed=TRUE)]
u.b<-teta.b.new[grep("u.b",names(teta.b),fixed=TRUE)]


if (control$print.epoch)
{       
print(u.t)
print(u.b)
cat("","\n\n")
}




###################################################################    
###################################################################
#optimization procedure for penalties##############################
#####goes by fix.point iteration###################################
###################################################################
###################################################################


########################################################################################################################
#define first the submatrices of the (penalized) fisher.information, corresponding to the random coefficients###########
########################################################################################################################
index.vector.t<-1:K.t  #for u.t.Baseline
if (p > 0) for (k in 1:p) index.vector.t<-c(index.vector.t,(K.t+(k-1)*K.t+1):(K.t+k*K.t))
index.vector.b<-1:K.b  #for u.b.Baseline
if (p > 0) for (k in 1:p) index.vector.b<-c(index.vector.b,(K.b+(k-1)*K.b+1):(K.b+k*K.b))



################################################################################
#chose optimization routine for smoothing.parameters; only fix.point goes#######
################################################################################
if (control$method !="fix") stop("only method='fix' allowed") else
{

for (iter.penalty in 1:control$niter.penalty)
{  
if (control$print.penalty) cat("iteration.penalty= ",iter.penalty,"\n\n") 

  
I.t.u.t<-I.t.pen[index.vector.t,index.vector.t]    #penalizied fisher.type matrix for t
inverse.I.t.u.t<-qr.solve(I.t.u.t,tol=1e-100)
I.b.u.b<-I.b.pen[index.vector.b,index.vector.b]    #penalizied fisher.type matrix for b
inverse.I.b.u.b<-qr.solve(I.b.u.b,tol=1e-100)



variance.penalty.t.old<-variance.penalty.t
variance.penalty.b.old<-variance.penalty.b
variance.penalty.old<-c(variance.penalty.t.old,variance.penalty.b.old)


#for t-direction
variance.t.Baseline<-as.vector((crossprod(u.t[grep("Baseline",names(u.t),fixed=TRUE)],Penalty.matrix.t%*%u.t[grep("Baseline",names(u.t),fixed=TRUE)]) + sum(diag(inverse.I.t.u.t[1:K.t,1:K.t]%*%Penalty.matrix.t[1:K.t,1:K.t])))/(K.t-2))
variance.penalty.t[1]<-variance.t.Baseline
#############################################################
#for-loop (in t) for the factor levels of covariables########
#############################################################
if ( p > 0)  #number of covariables
for (k in 1:p) variance.penalty.t[k+1]<-as.vector((crossprod(u.t[(k*K.t+1):((k+1)*K.t)],Penalty.matrix.t%*%u.t[(k*K.t+1):((k+1)*K.t)]) + sum(diag(inverse.I.t.u.t[(k*K.t+1):((k+1)*K.t),(k*K.t+1):((k+1)*K.t)]%*%Penalty.matrix.t)))/(K.t-2))
#update of penalty matrix in t
Penalty.matrix.teta.t[1:K.t,1:K.t]<-(1/variance.penalty.t[1])*Penalty.matrix.t
if (p > 0) for (k in 1:p) Penalty.matrix.teta.t[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)]<-(1/variance.penalty.t[k+1])*Penalty.matrix.t


#for b-direction
variance.b.Baseline<-as.vector((crossprod(u.b[grep("Baseline",names(u.b),fixed=TRUE)],Penalty.matrix.b%*%u.b[grep("Baseline",names(u.b),fixed=TRUE)]) + sum(diag(inverse.I.b.u.b[1:K.b,1:K.b]%*%Penalty.matrix.b[1:K.b,1:K.b])))/(K.b-2))
variance.penalty.b[1]<-variance.b.Baseline
#############################################################
#for-loop (in b) for the factor levels of covariables########
#############################################################
if (p > 0)
for (k in 1:p) variance.penalty.b[k+1]<-as.vector((crossprod(u.b[(k*K.b+1):((k+1)*K.b)],Penalty.matrix.b%*%u.b[(k*K.b+1):((k+1)*K.b)]) + sum(diag(inverse.I.b.u.b[(k*K.b+1):((k+1)*K.b),(k*K.b+1):((k+1)*K.b)]%*%Penalty.matrix.b)))/(K.b-2))
#update of penalty matrix in b
Penalty.matrix.teta.b[1:K.b,1:K.b]<-(1/variance.penalty.b[1])*Penalty.matrix.b
if (p > 0) for (k in 1:p) Penalty.matrix.teta.b[(K.b+(k-1)*K.b+1):(K.b+k*K.b),(K.b+(k-1)*K.b+1):(K.b+k*K.b)]<-(1/variance.penalty.b[k+1])*Penalty.matrix.b



variance.penalty<-c(variance.penalty.t,variance.penalty.b)


penalty.t<-1/variance.penalty.t
penalty.b<-1/variance.penalty.b
penalty<-c(penalty.t,penalty.b)
names(penalty)<-sub("variance","penalty",names(variance.penalty))



if (control$print.penalty) {print(penalty);cat("","\n\n")}
if (sum((variance.penalty-variance.penalty.old)^2) < control$tol.penalty) break   #stop for fixpoint iteration


}   #end for (iter.penalty in 1:control$niter.penalty)


} #end if (control$method !="fix") else


#stoping criterion for penalized.backfitting
variance.penalty.list[[iter.epoch+1]]<-variance.penalty
if (sum((variance.penalty.list[[iter.epoch+1]]-variance.penalty.list[[iter.epoch]])^2) < control$tol.epoch.variance & sum((teta.t.new-teta.t.old)^2)+sum((teta.b.new-teta.b.old)^2) < control$tol.epoch.theta) break #stop for outer loop (epoch=update of theta and penalties)



###################################################################
#marginal likelihood###############################################
###################################################################
g.inverse.t<-ginverse(Penalty.matrix.teta.t[index.vector.t,index.vector.t],tol=1e-100)
g.inverse.b<-ginverse(Penalty.matrix.teta.b[index.vector.b,index.vector.b],tol=1e-100)
#make.sym.t<-(t(g.inverse.t)+g.inverse.t)/2
#make.sym.b<-(t(g.inverse.b)+g.inverse.b)/2
#for t-direction
#log.lik.margin.t<-log.lik.t.pen(teta.t)-0.5*(sum(log(eigen(I.t.pen[index.vector.t,index.vector.t])$values)) + sum(log(eigen(make.sym.t)$values)))
log.lik.margin.t<-log.lik.t.pen(teta.t)-0.5*Re(sum(log(eigen(I.t.pen[index.vector.t,index.vector.t]%*%g.inverse.t)$values))) 
#for b-direction
#log.lik.margin.b<-log.lik.b.pen(teta.b)-0.5*(sum(log(eigen(I.b.pen[index.vector.b,index.vector.b])$values)) + sum(log(eigen(make.sym.b)$values)))
log.lik.margin.b<-log.lik.b.pen(teta.b)-0.5*Re(sum(log(eigen(I.b.pen[index.vector.b,index.vector.b]%*%g.inverse.b)$values)))



if (control$print.log.lik)
  {
cat("","\n\n")    
cat("log.lik.margin.t= ",log.lik.margin.t,"\n")
cat("log.lik.margin.b= ",log.lik.margin.b,"\n\n\n")
  }


log.lik.margin<-c(log.lik.margin.t=log.lik.margin.t,log.lik.margin.b=log.lik.margin.b)
log.lik.margin.t.list[[iter.epoch]]<-log.lik.margin.t
log.lik.margin.b.list[[iter.epoch]]<-log.lik.margin.b
  

} #end for (iter.epoch in 1:control$niter.epoch)




cat("","\n\n\n")

penalty.t<-1/variance.penalty.t
penalty.b<-1/variance.penalty.b
penalty<-c(penalty.t,penalty.b)
names(penalty)<-sub("variance","penalty",names(variance.penalty))


if (control$print.estimates)
  {
cat("resulting estimates: random effects and penalties","\n\n\n")  #last updates
print(u.t)
print(u.b)
cat("","\n\n")
print(penalty)
  }


#likelihood from the 1.iteration
log.lik.margin.start<-c(log.lik.margin.t.list[[1]],log.lik.margin.b.list[[1]])
names(log.lik.margin.start)<-c("log.lik.margin.t","log.lik.margin.b")


factor.names<-colnames(Design.variables)


rm(offset.list,offset.tij.list,unlist.offset.tij.list,offset.bi,y.poisson.list,y.poisson,I.t,I.b,I.t.pen,I.b.pen,Design.variables,Z.birth,Z.time,inverse.b,inverse.t,indizes)
for (i in 1:30) gc()



##################################################################
##################################################################
#calculate variances of the resulting estimates for theta#########
##################################################################
##################################################################


#1. combine together Design.matrix=(Design.matrix.t,Design.matrix.b)
help.lapply<-lapply(1:N,FUN=function(i) kronecker(rep(1,m.list[[i]]),t(Design.matrix.b[i,])))
help.matrix<-eval(parse(text=paste("rbind(",paste(paste("help.lapply[[",1:(length(help.lapply)-1),sep=""),"]],",sep="",collapse=""),"help.lapply[[length(help.lapply)]])",sep="")))

rm(help.lapply,Design.matrix.b,help.cut,m.list)
for (i in 1:30) gc()

Design.matrix<-cbind(Design.matrix.t,help.matrix)

rm(help.matrix,Design.matrix.t)
for (i in 1:30) gc()


#2. penalized and unpenalized Fiher.type matrix for the full design
I.tb.0<-crossprod(Design.matrix,Design.matrix*exp(as.vector(Design.matrix%*%c(teta.t,teta.b))+offset.model))
Penalty.matrix.tb[1:length(teta.t),1:length(teta.t)]<-Penalty.matrix.teta.t
Penalty.matrix.tb[(length(teta.t)+1):(length(teta.t)+length(teta.b)),(length(teta.t)+1):(length(teta.t)+length(teta.b))]<-Penalty.matrix.teta.b
I.tb<-I.tb.0+Penalty.matrix.tb



#3.estimate covariances through sandwich estimator
inverse.I.tb<-qr.solve(I.tb,tol=1e-100)
co.variance.teta<-inverse.I.tb%*%I.tb.0%*%inverse.I.tb
variance.teta<-diag(co.variance.teta)
co.variance.teta.t<-co.variance.teta[1:(K.t+p*K.t),1:(K.t+p*K.t)]
co.variance.teta.b<-co.variance.teta[(K.t+p*K.t+1):length(variance.teta),(K.t+p*K.t+1):length(variance.teta)]
names(variance.teta)<-c(names(teta.t),names(teta.b))
variance.teta.t<-variance.teta[1:(K.t+p*K.t)]
variance.teta.b<-variance.teta[(K.t+p*K.t+1):length(variance.teta)]


#4.degrees.of.freedom for components
#total
#df.total<-sum(diag(inverse.I.tb%*%I.tb.0))
#components specific
df.t.Baseline<-sum(diag((inverse.I.tb%*%I.tb.0)[1:K.t,1:K.t]))
df.b.Baseline<-sum(diag((inverse.I.tb%*%I.tb.0)[(K.t+p*K.t+1):(K.t+p*K.t+K.b),(K.t+p*K.t+1):(K.t+p*K.t+K.b)]))
if (p > 0) for (k in 1:p)
  {
assign(paste("df.t.",factor.names[k+1],sep=""),sum(diag((inverse.I.tb%*%I.tb.0)[(K.t+(k-1)*K.t+1):(K.t+k*K.t),(K.t+(k-1)*K.t+1):(K.t+k*K.t)])))
assign(paste("df.b.",factor.names[k+1],sep=""),sum(diag((inverse.I.tb%*%I.tb.0)[((p+1)*K.t+K.b+(k-1)*K.b+1):((p+1)*K.t+K.b+k*K.b),((p+1)*K.t+K.b+(k-1)*K.b+1):((p+1)*K.t+K.b+k*K.b)])))
  }
#write out d.f.
df<-matrix(0,nrow=p+1,ncol=2)
df[1,]<-c(df.t.Baseline,df.b.Baseline)
if (p > 0) for (k in 1:p) df[k+1,]<-c(get(paste("df.t.",factor.names[k+1],sep="")),get(paste("df.b.",factor.names[k+1],sep="")))
dimnames(df)<-list(c("Baseline",factor.names[-1]),c("t","b"))




#############################################
#Grid for survival.time and saison.time######
#############################################
grid.t<-seq(min(data.set$time),max(data.set$time),le=1000)
grid.b<-seq(min(data.set$birth),max(data.set$birth)+1,le=1000)
names(grid.t)<-1:length(grid.t)
names(grid.b)<-1:length(grid.b)

#specify knots
if (length(grep("num.knots.t",names(control)))==0) knots.t<-default.knots(grid.t) else knots.t<-default.knots(grid.t,num.knots)
if (length(grep("num.knots.b",names(control)))==0) knots.b<-default.knots(grid.b) else knots.b<-default.knots(grid.b,num.knots)

#specify spline.bases
B.grid.t<-splineDesign(knots=knots.t,ord=bs.order,x=grid.t)
B.grid.b<-splineDesign(knots=knots.b,ord=bs.order,x=grid.b)
B.grid.b[,1:bs.degree]<-B.grid.b[,1:bs.degree]+B.grid.b[,(ncol(B.grid.b)-(bs.degree-1)):ncol(B.grid.b)]
B.grid.b<-B.grid.b[,1:(ncol(B.grid.b)-bs.degree)]
#transform the Basis in b-direction (restriction B(1)=0)
c.vector<- -B.grid.b[as.numeric(names(grid.b)[grid.b==min(grid.b)][1]),-1]/unique(B.grid.b[grid.b==min(grid.b),1])
B.grid.b<-outer(B.grid.b[,1],c.vector)+B.grid.b[,-1]



#################################
#Confidence Bands################
#################################

#for t-direction
################
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

#for b-direction
################
co.variance.teta.b.Baseline<-co.variance.teta.b[1:K.b,1:K.b]
C.b.grid<-B.grid.b
variance.b.Baseline<-apply(C.b.grid,1,FUN=function(help.row) t(help.row)%*%co.variance.teta.b.Baseline%*%help.row)
deviation.b.Baseline<-qnorm(0.975)*sqrt(variance.b.Baseline)
#for covariables
if (p > 0)
for (k in 1:p)
  {
assign(paste("co.variance.teta.b.",factor.names[k+1],sep=""),co.variance.teta.b[(K.b+(k-1)*K.b+1):(K.b+k*K.b),(K.b+(k-1)*K.b+1):(K.b+k*K.b)])
assign(paste("variance.b.",factor.names[k+1],sep=""),apply(C.b.grid,1,FUN=function(help.row) t(help.row)%*%get(paste("co.variance.teta.b.",factor.names[k+1],sep=""))%*%help.row))
assign(paste("deviation.b.",factor.names[k+1],sep=""),qnorm(0.975)*sqrt(get(paste("variance.b.",factor.names[k+1],sep=""))))
}




penalties<-1/variance.penalty.list[[length(variance.penalty.list)]]
names(penalties)<-sub("variance","penalty",names(penalties))
random.coef<-c(u.t,u.b)
var.random<-c(variance.teta.t[grep("u",names(variance.teta.t))],variance.teta.b[grep("u",names(variance.teta.b))])
names(var.random)<-sub("u","variance.u",names(var.random))




######################################################################
#varying coefficients for Baseline and covariables####################
######################################################################
alpha.t.Baseline<-B.grid.t%*%u.t[grep("Baseline",names(u.t),fixed=TRUE)]
alpha.b.Baseline<-B.grid.b%*%u.b[grep("Baseline",names(u.b),fixed=TRUE)]

if (p > 0)
{
for (k in 1:p)
  {
assign(paste("alpha.t.",factor.names[k+1],sep=""),B.grid.t%*%u.t[grep(paste("u.t.",factor.names[k+1],sep=""),names(u.t),fixed=TRUE)])  
assign(paste("alpha.b.",factor.names[k+1],sep=""),B.grid.b%*%u.b[grep(paste("u.b.",factor.names[k+1],sep=""),names(u.b),fixed=TRUE)])
}
}



#########################################################################################
#########################################################################################
#write out varying coefficient, deviation and grid vectors in frames#####################
#########################################################################################
#########################################################################################
#varying coefficients
list.t.frame<-list()
list.b.frame<-list()
t.frame<-NULL;b.frame<-NULL
if (p > 0)
{  
for (k in 1:p)
  {list.t.frame[[k]]<-get(paste("alpha.t.",factor.names[k+1],sep=""))
   list.b.frame[[k]]<-get(paste("alpha.b.",factor.names[k+1],sep=""))
   t.frame<-cbind(t.frame,list.t.frame[[k]])
   b.frame<-cbind(b.frame,list.b.frame[[k]])
  }
varying.frame<-data.frame(alpha.t.Baseline,t.frame,alpha.b.Baseline,b.frame)
names(varying.frame)<-c("alpha.t.Baseline",paste("alpha.t.",factor.names[2:(p+1)],sep=""),"alpha.b.Baseline",paste("alpha.b.",factor.names[2:(p+1)],sep=""))
} else
{varying.frame<-data.frame(alpha.t.Baseline,alpha.b.Baseline)
 names(varying.frame)<-c("alpha.t.Baseline","alpha.b.Baseline")
}



#deviations
list.t.frame<-list()
list.b.frame<-list()
t.frame<-NULL;b.frame<-NULL
if (p > 0)
{  
for (k in 1:p)
  {list.t.frame[[k]]<-get(paste("deviation.t.",factor.names[k+1],sep=""))
   list.b.frame[[k]]<-get(paste("deviation.b.",factor.names[k+1],sep=""))
   t.frame<-cbind(t.frame,list.t.frame[[k]])
   b.frame<-cbind(b.frame,list.b.frame[[k]])
  } 
deviation.frame<-data.frame(deviation.t.Baseline,t.frame,deviation.b.Baseline,b.frame)
names(deviation.frame)<-c("deviation.t.Baseline",paste("deviation.t.",factor.names[2:(p+1)],sep=""),"deviation.b.Baseline",paste("deviation.b.",factor.names[2:(p+1)],sep=""))
} else
{deviation.frame<-data.frame(deviation.t.Baseline,deviation.b.Baseline)
names(deviation.frame)<-c("deviation.t.Baseline","deviation.b.Baseline")
}

#grids
grid.frame<-data.frame(grid.t=grid.t,grid.b=grid.b)


#output elements
list(random.coef=random.coef,penalty=penalty,var.random=var.random,log.lik.margin.start=log.lik.margin.start,log.lik.margin=log.lik.margin,df=df,df.total=sum(df),niter.epoch=iter.epoch,varying.frame=varying.frame,deviation.frame=deviation.frame,grid.frame=grid.frame,p=p,factor.names=factor.names[-1])



}#end of function

 
