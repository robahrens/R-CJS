###################################
#	Cormack-Jolly-Seber Model 
#	Robert Ahrens
#	University of Florida
#	Version 1.0 December 2014
################################### 

#functions used in simulating data and parameter rescaling
join=function(x)
{
	tmp=""
	for (i in 1:nsample)
	{
		tmp=paste(tmp,status[x,i],sep="")
	}
	return(tmp)
}
extractdiags = function(x) tryCatch(diag(solve(x)),error = function(e) {print("Lapack routine dgesv: system is exactly singular");rep(NA,nsurv+npcap)})
logit=function(x)log(x/(1-x))
invlogit=function(x)exp(x)/(exp(x)+1)

#Global values used to generate data and estimate
nsample=5 #number of recapture events
survblock=c(1,1,1,2,2) # this is a vector for how survival should be blocked in time(the max blocks is nsample)
pcapblock=c(1,1,2,2,2) # this is a vector of how pcap shoul be blocked in time(the max blocks is nsample)
nsurv=2 #this must be specified by the user and line up with the number of blocks used in the vector above
npcap=2 #this must be specified by the user and line up with the number of blocks used in the vector above
parlab=c(paste(rep("s",nsurv),1:nsurv,sep=""),paste(rep("p",nsurv),1:npcap,sep=""))
filename="recap.inp"
createdata=TRUE
if (createdata)
{
	#specify the number of tagged individuals for data generation
	nind=100
	#set up survival and pcap for data genration these will be distributed according to the blocking
	tbsurv=c(0.7,0.8)
	tbpcap=c(0.3,0.5)
	tsurv=vector(length=nsample)
	tpcap=vector(length=nsample)
	for (i in 1:nsample)
	{
		tsurv[i]=tbsurv[survblock[i]]
		tpcap[i]=tbpcap[pcapblock[i]]
	}
	#create matricies to keep track of status of individuals
	status=matrix(rep(1,nind*nsample),ncol=nsample)#initialized matrix everyone begins all times alive
	alive=matrix(runif(nind*nsample),ncol=nsample)#coin toss to see if and individual dies
	sampled=matrix(rep(1,nind*nsample),ncol=nsample)#initialized as everyone is sampled everytime
	seen=matrix(runif(nind*nsample),ncol=nsample)#coin toss to see if an individual is sampled
	firstkill=vector(length=nind)
	firstkill=firstkill+(nsample+1)
	for (i in 1:nind)
	{
		for (j in 1:nsample)
		{
			if(seen[i,j]>tpcap[j])sampled[i,j]=0
			if(alive[i,j]>tsurv[j])
			{
				status[i,j]=0
				if (j<firstkill[i])firstkill[i]=j
			}
		}
	}
	for (i in 1:nind)
	{
		if(firstkill[i]<(nsample+1))
		{
			xx=firstkill[i]:nsample 
			status[i,xx]=0
		}
	}
	status=status*sampled
	
	# this compacts the data and creates a MARK type data set
	status2=table(sapply(1:nind,join))
	uniques=length(status2)	
	dfile=matrix(paste(names(status2),status2[1:uniques],sep=";"),ncol=1)
	write.table(dfile, file = filename,quote = FALSE,row.names = FALSE,
            col.names = FALSE)
	#This is where a data file would need to be read in as status3 and counts
}
#This is where a data file would need to be read in as status3 and counts
dfile=read.table(file = filename,colClasses = c("character","numeric"),sep=";")
status2=dfile[,1]
counts=dfile[,2]
uniques=length(counts)
nind=sum(counts)
status3=matrix(ncol=nsample,nrow=uniques)
for (i in 1:uniques)
{
	status3[i,]=as.numeric(unlist(strsplit(status2[i],"")))
}
lastrecap=apply(status3,1,FUN=function(x)max(which(x==1)))
lastrecap[which(lastrecap==-Inf)]=0

theta=c(logit(rep(0.5,nsurv)),logit(rep(0.5,npcap)))
#theta=c(logit(tbsurv),logit(tbpcap))
nll=function(theta)
{
	surv=vector(length=nsample)
	pcap=vector(length=nsample)
	for (i in 1:nsample)
	{
		surv[i]=theta[survblock[i]]
		pcap[i]=theta[pcapblock[i]+nsurv]
	}
	surv=invlogit(surv)
	pcap=invlogit(pcap)
	LL=vector(length=uniques)
	probs=vector(length=uniques)
	one=surv*pcap
	zero1=surv*(1-pcap)
	zero2=(1-surv)
	#alive looping
	for (i in which(lastrecap>0))
	{
		for (j in 1:lastrecap[i])
		{
			if (status3[i,j]==1)
			{
				if(j==1)probs[i]=one[j] else probs[i]=probs[i]*one[j]
			}
			if (status3[i,j]==0)
			{
				if(j==1)probs[i]=zero1[j] else probs[i]=probs[i]*zero1[j]
			}
		}
	}
	#dead or unseen looping 
	for (i in which(lastrecap<nsample))
	{		
			nzeros=length((lastrecap[i]+1):nsample)
			tmpprobs2=zero2[lastrecap[i]+1]
			cnt=0
			for (j in 1:nzeros)
			{
				xx=(lastrecap[i]+1):(nsample-cnt)
				tmpprobs1=cumprod(zero1[xx])[length(xx)]
				if (cnt>0)tmpprobs1=tmpprobs1*zero2[nsample-(cnt-1)]
				cnt=cnt+1
				tmpprobs2=tmpprobs2+tmpprobs1
			}
			if (probs[i]==0)probs[i]=tmpprobs2 else probs[i]=probs[i]*tmpprobs2
	}
	LL=counts*log(probs)
	return(-sum(LL))
}

fit=optim(theta,nll,hessian=TRUE)
mle=invlogit(fit$par)
diags=extractdiags(fit$hessian) 
lower95=invlogit(fit$par-1.96*sqrt(diags))
upper95=invlogit(fit$par+1.96*sqrt(diags))
k=nsurv+npcap
aic=2*fit$value+2*k 
aicc=aic+2*k*(k+1)/(nind-k-1)
svddiag=svd(fit$hessian)$d
epars=k-length(which(svddiag<0.001*max(svddiag)))
if(sum(svddiag)==0)epars=0
output=data.frame(Labels=parlab,mle=mle,lower95=lower95,upper95=upper95)
if(epars<k)print(paste("Number of estimable paramters for this configuration is",epars))
output
aic
aicc





