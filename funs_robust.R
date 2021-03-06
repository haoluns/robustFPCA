

findmean = function(observed,timepoints,threshold=1e-3,minit=3,gamma=0){

w = lapply(1:length(timepoints),function(i){
nn=length(timepoints[[i]])
rep(1/nn,nn)
})

yem = do.call(c,observed)
deltaem= do.call(c,w)

xalpha = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati

})%>%do.call(rbind,.)



beta1 = lmrobust(yem,as.matrix(xalpha),deltaem)

	
pc_fit = fd(beta1, spline_basis)

beta1 = coef(pc_fit)%>%as.numeric
return(list(beta=beta1,pc_fit = pc_fit))

}


first_FPC  = function(beta1,observed,timepoints,threshold=1e-3,maxit=50,gamma=0){


delta = lapply(1:length(timepoints),function(i){
nn=length(timepoints[[i]])
rep(1/nn,nn)
})

minit=3
thresh = 1
it = 1
beta1_before = rep(0, length(beta1))
value = -1e14
pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta1 = coef(pc_fit)%>%as.numeric
R = inprod(spline_basis, spline_basis,2,2)

while (thresh >  threshold|it<minit){
	
	beta1_before = beta1
	value_before = value
	
alpha_fit = function(subj_index){
	(timei =timepoints[[subj_index]])
	deltai = delta[[subj_index]]
	xmati = eval.fd(timei, pc_fit)
	lmrobust(observed[[subj_index]],as.matrix(xmati,ncol=1),deltai)	
}

sfit= sapply(1:length(observed),alpha_fit)


yem = do.call(c,observed)
deltaem= do.call(c,delta)

xalpha = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit[subj_index]

})%>%do.call(rbind,.)

beta1 = lmrobust(yem,as.matrix(xalpha),deltaem)


pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)
beta1 = coef(pc_fit)%>%as.numeric
thresh =  max(abs(beta1_before-beta1))
it = it+1
if(it%%2==0) {print(it);print(as.numeric(thresh))}
if(it>maxit){
break
}

}
pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta1 = coef(pc_fit)%>%as.numeric
#rfit= lapply(1:length(observed),residual_fit)%>%do.call(c,.)
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta1,pc_fit = pc_fit,sfit = sfit,thresh = thresh,it=it))
}



third_FPC_conditional  = function(beta3, pc_index, observed, timepoints, betalist =previous_beta , threshold=1e-4,maxit=50,gamma=0){


delta = lapply(1:length(timepoints),function(i){
nn=length(timepoints[[i]])
rep(1/nn,nn)
})


minit = 1
#if(missing(pc_index)) pc_index= length(betalist)+3
thresh = 1
it = 1
R = inprod(spline_basis, spline_basis,2,2)
E = inprod(spline_basis, spline_basis,0,0)

pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})
value = -1e14

tempdelta <- function(x){
sum(1-x)
}

observed2 =observed
observed2[which(sapply(delta,tempdelta)<=pc_index)]=NULL
timepoints2 = timepoints
timepoints2[which(sapply(delta,tempdelta)<=pc_index)]=NULL
delta[which(sapply(delta,tempdelta)<=pc_index)]=NULL

beta3_before = rep(0, length(betalist[[1]]))
while (thresh >  threshold | it<minit){
	
	
	
	beta3_before = beta3
	value_before = value

alpha_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit)
	deltai=delta[[subj_index]]
	
	lmrobust(observed2[[subj_index]],cbind(xmat_previous,as.matrix(xmati,ncol=1)),deltai)
	#tryCatch(lmrobust(observed2[[subj_index]],cbind(xmat_previous,as.matrix(xmati,ncol=1)),deltai), error = function(e) rep(999999,ncol(xmat_previous)+1))
	
}

sfit= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)




		
N = sapply(observed2,length)%>%sum

yem = lapply(1:length(observed2), function(i){
(timei =timepoints2[[i]])
yfits_previous = lapply(1:length(pc_fits_previous), function(x){
	sfit[i,x]*eval.fd(timei,pc_fits_previous[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observed2[[i]] - yfits_previous)#/sqrt(length(timei))/sqrt(N)
})%>%do.call(c,.)


deltaem= do.call(c,delta)

xalpha = lapply(1:length(observed2),function(subj_index){

	(timei =timepoints2[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit[subj_index,ncol(sfit)])#/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


#A = xalpha
#qmat = 2*(t(A)%*%A)
#pmat = as.numeric(-2*t(yem)%*%A)
betamat=  do.call(rbind,betalist)
cmat = rbind(betamat%*%E,-betamat%*%E)

#lmrobust(yem,cbind(xmat_previous,as.matrix(xmati,ncol=1)),deltai)	


beta3 = lmrobust(yem,as.matrix(xalpha),deltaem,cmat)



pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

thresh =  max(abs(beta3_before-beta3))


it = it+1
if(it%%2==0) {print(it);print(thresh)}
if(it>maxit){
break
}

}
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, thresh = thresh,it=it))
}


