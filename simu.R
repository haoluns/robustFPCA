library(fdapace)

simusparse<-function(outlierexist=FALSE,kk=1,epsilon1 = 0.2,model=1){
f1<-function(x){sqrt(2)*sin(2*pi*x)}
f2<-function(x){sqrt(2)*cos(2*pi*x)}
muf <- function(x){
0.5 + sin(6*pi*x)*exp(-2*x)
}

timegrid = seq(0,1,by=0.01)

ssize=200
a1=rnorm(ssize,0,4)
a2=rnorm(ssize,0,1)




timepoints = lapply(1:ssize,function(i){
curtimegrid = sample(timegrid,10)
sort(curtimegrid)
})


observed = lapply(1:ssize,function(i){
tm=timepoints[[i]]
mu=sapply(tm,muf)
pc1add=sapply(tm,f1)*a1[i]
pc2add= sapply(tm,f2)*a2[i]
err=rnorm(1,mean=0,sd=0.01)
as.numeric(mu+pc1add+pc2add+err)
})



if(model==1){
if(outlierexist == TRUE){
for(i in 1:length(observed)){
xx=rbinom(1,1,epsilon1)
if(xx==1){
epsilon2 = 0.3
nnn=length(observed[[i]])
z1=rnorm(nnn,5,0.01)
z2=rbinom(nnn,1,epsilon2)
observed[[i]] = observed[[i]]+z1*z2
}
}
}
}



lippo=t(observed%>%do.call(rbind,.))


norm <- function(a) sqrt(sum(a^2))

# prepare data
x <- t(lippo)
te <- seq(0,1.0,by=0.01)



# prepare the spline basis (with 20 functions)
mesh <- te
n <- nrow(x)
dimension.Bspline <- 20 

library(fdapace)

res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",methodSelectK = 2,nRegGrid=101))



library(fda)


createfd<-function(te,y){
spline_basis=create.bspline.basis(rangeval=c(min(te),max(te)),nbasis=20,norder=4)
beta1 = coef(Data2fd(argvals = te,y=y,spline_basis))
pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
return(pc_fit)
}


fpc1_pace = createfd(seq(0,1,by=0.01),res_pace$phi[,1])


fpc2_pace = createfd(seq(0,1,by=0.01),res_pace$phi[,2])

fpc1_true=createfd(seq(0,1,by=0.01),sapply(seq(0,1,by=0.01),f1))
fpc2_true=createfd(seq(0,1,by=0.01),sapply(seq(0,1,by=0.01),f2))


plot(fpc1_pace)
plot(fpc1_true)
#if (eval.fd(0.3,fpc1_pace)<0)fpc1_pace=-fpc1_pace





library(fda)
timepts=seq(0,1,by=0.01);
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2;
spline_basis=create.bspline.basis(rangeval=c(0,1),nbasis=20,norder=4)

coef_mat0 = coef(Data2fd(argvals = res_pace$workGrid,y=res_pace$phi,spline_basis))
beta1=coef_mat0[,1]





### find mean function 

ptm=proc.time()


meanfit = findmean(observed=observed, timepoints=timepoints, minit=6,gamma=0,threshold=1e-5)

mu = eval.fd(seq(0,1,by=0.01), meanfit$pc_fit)



observedcenter = lapply(1:length(observed),function(i){
	 (observed[[i]]-eval.fd(timepoints[[i]], meanfit$pc_fit))[,1]
})



### first fpc
pc1s = first_FPC(coef_mat0[,1],observed=observedcenter, timepoints=timepoints,gamma=0,threshold=0.000001)


##second fpc
previous_beta = list()
previous_beta[[1]] = pc1s$beta
beta3=pc1s$beta
 
pc2s=third_FPC_conditional(coef_mat0[,2], 2, observedcenter, timepoints, betalist =previous_beta , threshold=1e-3,gamma=0)

rfpcatime=proc.time() - ptm


fpc1_r = fd(pc1s$beta, spline_basis)
fpc2_r = fd(pc2s$beta, spline_basis)





output=c(min(inprod(fpc1_pace-fpc1_true,fpc1_pace-fpc1_true),inprod(fpc1_pace+fpc1_true,fpc1_pace+fpc1_true)),
min(inprod(fpc2_pace-fpc2_true,fpc2_pace-fpc2_true),inprod(fpc2_pace+fpc2_true,fpc2_pace+fpc2_true)),

min(inprod(fpc1_r-fpc1_true,fpc1_r-fpc1_true),inprod(fpc1_r+fpc1_true,fpc1_r+fpc1_true)),
min(inprod(fpc2_r-fpc2_true,fpc2_r-fpc2_true),inprod(fpc2_r+fpc2_true,fpc2_r+fpc2_true))
)

return(output)

}





spline_basis=create.bspline.basis(rangeval=c(0,1),nbasis=20,norder=4)

library(data.table)


	



output=matrix(ncol=6,nrow=0)
for(i in 1:100){
kk=0;epsilon1= 0.3; model=1
o=try(simusparse(outlierexist=TRUE,kk,epsilon1, model))
while('try-error' %in% class(o)){
o=try(simusparse(outlierexist=TRUE,kk,epsilon1, model))
}
output=rbind(output,o)
}
colnames(output)=c('pace1','pace2','s1','s2','r1','r2')
write.table(output,paste0('model',model,'epsilon1',epsilon1*10,'kk',kk,'.csv'))


