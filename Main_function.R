OS<-function(S,nmaj,t)
{
z=ncol(S)
ui=vector("numeric",length=z-1)
Snew=matrix(ncol=z)
m=nrow(S)
soft=vector(mode="list",length=t)
scmft=vector(mode="list",length=t)
for(i in 1:t) 
{
	sort= S[sort.list(S[,z]), ]
	Sm <- S[which(S[,ncol(S)]==(i+1)),]
	print("choosesamples")
	Scm=choosesamples(sort,Sm)		#choosesamples
	w=Scm[,(z+1)]
	Scm=Scm[,-(z+1)]	
	scmft[[i]]=Scm	
	ni=nrow(Scm)
	ui=colMeans(Scm[,-z])
	Zi=Scm[,-z]	
	for(h in 1:ni)
	{Zi[h,]=Zi[h,]-ui}
	print("pca")
	F=Pca(Zi)
	T=F[1:(nrow(F)-(z-1)),]			#after PCA		
	Eigenvectors=F[(nrow(F)-(z-2)):nrow(F),]
	V=diag(cov(T))		
	Orate=nmaj-ni
	print("MDO")
	Stemp=MDO(T,V,Orate,w)			#new examples
	soft[[i]]=Stemp
	invev=t(Eigenvectors)
	Stempo=Stemp%*%invev
	Stom=NULL
	for(l in 1:ncol(Stempo))
	{
	  mean=Stempo[,l]+ui[l]
	  Stom=cbind(Stom,mean)
	
	}
	
	cl=rep(i+1,length.out=nrow(Stom))
	Stoma=cbind(Stom,cl)			#final examples
	filename<-paste((i+1)*10,".csv",sep="")
	write.table(Stoma,filename,sep=",")
	Snew=rbind(Snew,Stoma)
	print(i)
}
Snew=Snew[-1,]
return(Snew)
}
#/////////////////////////////////////////////////

choosesamples<-function(S,Sm)
{
z=ncol(S)
Scm=matrix(,ncol=z)
weights=vector()
temp=vector(length=10)
m=nrow(S)
dist=matrix(nrow=m,ncol=3)
count=vector()
check=matrix(nrow=nrow(Sm),ncol=10)
mi=nrow(Sm)
for(j in 1:mi)
{
	print(j)
	xjc=Sm[j,]
	xj=Sm[j,-z]
	for(k in 1:m)
	{
		xk=S[k,-z]
		dist[k,1]=k	
		dist[k,2]=dist(rbind(xj,xk))
		dist[k,3]=S[k,z]
	}
	
	s= dist[sort.list(dist[,2]), ]
	
	temp=s[1:10,1]
	
	for(g in 1:10)
	{	q=temp[g]
		check[j,g]=S[q,z]	
	}
	
	num=0
	for(l in 1:10)
	{
		c=temp[l]
		if(S[c,z]==Sm[1,z])
		{num=num+1}	
	}
	
	count=c(count,num)
	
	if(num>=3)
	{
	Scm=rbind(Scm,xjc)
	weights=c(weights,num/10)
	}
	
	
	
}
Scm=Scm[-1,]
w=weights/sum(weights)
rS=cbind(Scm,w)
return(rS)
}

#\\\\\\\\\\\\\\\\\\\\\\\\\\MDO

MDO<-function(T,V,Orate,weights)
{
numa=ncol(T)
Stemp=matrix(ncol=numa)
while(Orate!=0)
{
	x=runif(1,0,1)
	l=1	
	p=weights[l]	
	while(x>p)
	{
		l=l+1
		p=p+weights[l]		
	}
	xnew=T[l,]
	X=xnew^2
	Q=X/V
	Q[!is.finite(Q)]=0
	alpha=sum(Q)
	alphaV=alpha*V	
	s=0
	xnew=vector(length=numa)
	for(j in 1:(numa-1))
	{
		d=sqrt(alphaV[j])
		r=runif(1,-1*d,d)
		qw=(r^2)/alphaV[j]
		qw[!is.finite(qw)]=0
		s=s+(qw)
		xnew[j]=r
	}
	if(s>1) 
	{next}
	lastfeaval=sqrt(alphaV[numa]*(1-s))
  xnew[numa]=sample(c(1,-1),1)*lastfeaval
  Stemp=rbind(Stemp,xnew)
Orate=Orate-1
print(Orate)
}
Stemp=Stemp[-1,]
return(Stemp)
}

#/////////////////////////PCA

Pca<-function(data)
{
  cov=cov(data)
  eval=eigen(cov)$values
  evec=eigen(cov)$vectors
  fd=data%*%evec
  fdev=rbind(fd,evec)
  return(fdev)
}
