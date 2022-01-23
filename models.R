##################  LOMAX Distribution  #####################################################

rm(list = ls())

library(optimx)


data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)


LOMAX.lik<-function(vector1,freq){
  alpha=vector1[1]
  lamda=vector1[2]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=alpha*(lamda^(alpha))
    #print(term1)
    term2=(i+lamda)^(-(alpha+1))
    #print(term2)
    
    total_value = total_value + term1*term2
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  #print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(alpha)
  
  likterm3=n*alpha*log(lamda)
  
  likterm4=(alpha+1)*(t(freq)%*%(log(y+lamda)))
  
  final_likterm<-sum(c(likterm1,likterm2,likterm3,-likterm4))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,1),LOMAX.lik, freq = freq)

#print(output11)


probability.fun_lomax<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estlamda=vector3[2]
  
  x=datafile
  
  sum2=0
  large_N=4000
  
  for(i in 1:large_N)
  {
    
    
    pdf_term1=estalpha*(estlamda^(estalpha))
    #print(term1)
    pdf_term2=(i+estlamda)^(-(estalpha+1))
    
    sum2 = sum2 + pdf_term1*pdf_term2
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  #print(C2)
  
  pdf_term3=estalpha*(estlamda^(estalpha))
  #print(term1)
  pdf_term4=(x+estlamda)^(-(estalpha+1))
  
  final_pdf_term=C2*pdf_term3*pdf_term4
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_lomax(c(1000,1),y)
outfun<-probability.fun_lomax(c(output11$par[1],output11$par[2]),y)
#print('probability_sum')
#print(sum(outfun))
#plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_twitter_lomax.csv')


################################# GLM Type-I Distribution #############################################



rm(list = ls())

library(optimx)


data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)



########### calculation of CV ######################

total_data_pts=sum(freq)
mean1=(t(freq)%*%y)/total_data_pts
sd1=sqrt((t(freq)%*%(y-mean1)^2)/total_data_pts)
print('mean1')
print(mean1)
print('sd1')
print(sd1)

data_set_total=rep(y,freq)
mean2=mean(data_set_total)
sd2=sd(data_set_total)
print('mean2')
print(mean2)
print('sd2')
print(sd2)

CV1=sd1/mean1
print('CV1')
print(CV1)
CV2=sd2/mean2
print('CV2')
print(CV2)


########################################################




version1.lik<-function(vector1,freq){
  alpha=vector1[1]
  beta=vector1[2]
  gama=vector1[3]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=alpha*((1+i/gama)^(-alpha-1))
    #print(term1)
    term2=(1+(beta/((1+log(1+i/gama))^2)))
    #print(term2)
    term3=exp(-alpha*beta*(log(1+i/gama)/(1+log(1+i/gama))))
    
    total_value = total_value + term1*term2*term3
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  #print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(alpha)
  
  likterm3=(alpha+1)*(t(freq)%*%(log(1+y/gama)))
  
  likterm4=t(freq)%*%(log(1+(beta/(1+(log(1+y/gama)))^2)))
  
  likterm5=alpha*beta*(t(freq)%*%(log(1+y/gama)/(1+log(1+y/gama))))
  
  
  final_likterm<-sum(c(likterm1,likterm2,-likterm3,likterm4,-likterm5))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,0,1),version1.lik, freq = freq)

#print(output11)


###########################################


probability.fun_version1<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estbeta=vector3[2]
  estgama=vector3[3]
  
  
  x=datafile
  
  
  sum2=0
  large_N=4000
  
  for(i in 1:large_N)
  {
    
    
    pdf_term1=estalpha*((1+i/estgama)^(-estalpha-1))
    #print(term1)
    pdf_term2=(1+(estbeta/((1+log(1+i/estgama))^2)))
    #print(term2)
    pdf_term3=exp(-estalpha*estbeta*(log(1+i/estgama)/(1+log(1+i/estgama))))
    
    
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  #print(C2)
  
  pdf_term4=estalpha*((1+x/estgama)^(-estalpha-1))
  #print(term4)
  pdf_term5=(1+(estbeta/((1+log(1+x/estgama))^2)))
  #print(term5)
  pdf_term6=exp(-estalpha*estbeta*(log(1+x/estgama)/(1+log(1+x/estgama))))
  
  
  fina1_pdf_term=C2*pdf_term4*pdf_term5*pdf_term6
  return(fina1_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_version4(c(1,0.5,1),y)
outfun<-probability.fun_version1(c(output11$par[1],output11$par[2],output11$par[3]),y)
#print('probability_sum')
#print(sum(outfun))
#plot(y,outfun)



finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_twitter_version1.csv')



################################### GLM Type-II Distribution ########################################




rm(list = ls())

library(optimx)


data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)



version2.lik<-function(vector1,freq){
  alpha=vector1[1]
  beta=vector1[2]
  gama=vector1[3]
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=(1+(beta/(1+(i/gama))))
    #print(term1)
    term2=(alpha/(i+gama))
    #print(term2)
    term3=(1+(i/gama))^(-alpha)
    #print(term3)
    term4=exp(-alpha*beta*((i/gama)/(1+(i/gama))))
    
    total_value = total_value + term1*term2*term3*term4
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  #print(C)
  
  likterm1=n*log(C)
  
  likterm2=t(freq)%*%(log((gama+y+beta*gama)/(gama+y)))
  
  likterm3=t(freq)%*%(log((alpha)/(gama+y)))
  
  likterm4=alpha*(t(freq)%*%(log(1+(y/gama))))
  
  likterm5=alpha*beta*(t(freq)%*%(y/(gama+y)))
  
  
  final_likterm<-sum(c(likterm1,likterm2,likterm3,-likterm4,-likterm5))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,0,1),version2.lik, freq = freq)

#print(output11)


###########################################


probability.fun_version2<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estbeta=vector3[2]
  estgama=vector3[3]
  
  
  x=datafile
  
  
  sum2=0
  large_N=4000
  
  for(i in 1:large_N)
  {
    
    
    pdf_term1=(1+(estbeta/(1+(i/estgama))))
    #print(term1)
    pdf_term2=(estalpha/((i+estgama)))
    #print(term2)
    pdf_term3=(1+(i/estgama))^(-estalpha)
    #print(term3)
    pdf_term4=exp(-estalpha*estbeta*((i/estgama)/(1+(i/estgama))))
    
    
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3*pdf_term4
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  #print(C2)
  
  
  pdf_term5=(1+(estbeta/(1+(x/estgama))))
  #print(term1)
  pdf_term6=(estalpha/((x+estgama)))
  #print(term2)
  pdf_term7=(1+(x/estgama))^(-estalpha)
  #print(term3)
  pdf_term8=exp(-estalpha*estbeta*((x/estgama)/(1+(x/estgama))))
  
  
  fina1_pdf_term=C2*pdf_term5*pdf_term6*pdf_term7*pdf_term8
  return(fina1_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_version2(c(1,0,1),y)
outfun<-probability.fun_version2(c(output11$par[1],output11$par[2],output11$par[3]),y)
#print('probability_sum')
#print(sum(outfun))
#plot(y,outfun)



finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_twitter_version2.csv')




########################################## GLM Type-III distribution #########################




rm(list = ls())

library(optimx)


data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)




version3.lik<-function(vector1,freq){
  alpha=vector1[1]
  beta=vector1[2]
  gama=vector1[3]
  
  
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=1+((beta*(log(1+i/gama)))/(i/gama))
    #print(term1)
    term2=((i/gama)/(1+i/gama))^(beta)
    #print(term2)
    term3=(alpha/(i+gama))
    #print(term3)
    term4=exp(-alpha*log((1+i/gama)*(((i/gama)/(1+i/gama))^beta)))
    
    total_value = total_value + term1*term2*term3*term4
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  #print(C)
  
  likterm1=n*log(C)
  
  likterm2=t(freq)%*%(log((y/gama + beta*log(1+y/gama))/(y/gama)))
  
  likterm3=beta*(t(freq)%*%(log((y/gama)/(1+y/gama))))
  
  likterm4=t(freq)%*%(log(alpha/(y+gama)))
  
  likterm5=alpha*(t(freq)%*%(log((1+y/gama)*(((y/gama)/(1+y/gama))^(beta)))))
  
  
  final_likterm<-sum(c(likterm1,likterm2,likterm3,likterm4,-likterm5))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,0,1),version3.lik, freq = freq)

#print(output11)


###########################################


probability.fun_version3<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estbeta=vector3[2]
  estgama=vector3[3]
  
  
  x=datafile
  
  
  sum2=0
  large_N=4000
  
  for(i in 1:large_N)
  {
    
    pdf_term1=1+((estbeta*(log(1+i/estgama)))/(i/estgama))
    #print(term1)
    pdf_term2=((i/estgama)/(1+i/estgama))^(estbeta)
    #print(term2)
    pdf_term3=(estalpha/(i+estgama))
    #print(term3)
    pdf_term4=exp(-estalpha*log((1+i/estgama)*(((i/estgama)/(1+i/estgama))^estbeta)))
    
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3*pdf_term4
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  print(C2)
  
  
  pdf_term5=1+((estbeta*(log(1+x/estgama)))/(x/estgama))
  #print(term1)
  pdf_term6=((x/estgama)/(1+x/estgama))^(estbeta)
  #print(term2)
  pdf_term7=(estalpha/(x+estgama))
  #print(term3)
  pdf_term8=exp(-estalpha*log((1+x/estgama)*(((x/estgama)/(1+x/estgama))^estbeta)))
  
  
  fina1_pdf_term=C2*pdf_term5*pdf_term6*pdf_term7*pdf_term8
  return(fina1_pdf_term)
  
  #return(probability_array)
  
}


#outfun<-probability.fun_version4(c(2,5,1),y)
outfun<-probability.fun_version3(c(output11$par[1],output11$par[2],output11$par[3]),y)
#print('probability_sum')
#print(sum(outfun))
#plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_twitter_version3.csv')




##############################  GLM Type-IV Distribution #############################



rm(list = ls())

library(optimx)


data<-read.csv('final_indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)


version4.lik<-function(vector1,freq){
  alpha=vector1[1]
  beta=vector1[2]
  sigma=vector1[3]
  
  n<-sum(freq)
  
  total_value=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    
    
    term1=alpha/sigma
    #print(term1)
    term2=(log((i/sigma) +1)+1+beta)/((i/sigma)+1)
    #print(term2)
    term3=((log((i/sigma)+1))^(beta))/((log((i/sigma)+1)+1)^(beta+1))
    #print(term3)
    term4=exp(-alpha*(((log((i/sigma)+1))^(beta+1))/((log((i/sigma)+1)+1)^(beta))))
    #print(term4)
    
    total_value = total_value + term1*term2*term3*term4
    #print(total_value)
    
    #print("---------")
  }
  
  C=1/total_value
  
  #print(C)
  
  likterm1=n*log(C)
  
  likterm2=n*log(alpha)
  
  likterm3=t(freq)%*%log(y+sigma)
  
  likterm4=t(freq)%*%(log(log((y/sigma)+1)+1+beta))
  
  likterm5=beta*(t(freq)%*%(log(log((y/sigma)+1))))
  
  likterm6=(beta+1)*(t(freq)%*%(log(log((y/sigma)+1)+1)))
  
  likterm7=alpha*(t(freq)%*%(((log((y/sigma)+1))^(beta+1))/((log((y/sigma)+1)+1)^(beta))))
  
  final_likterm<-sum(c(likterm1,likterm2,-likterm3,likterm4,likterm5,-likterm6,-likterm7))
  
  return(-final_likterm)
  
}

output11<-optim(c(1,0,1),version4.lik, freq = freq)

#print(output11)


###########################################


probability.fun_version4<- function(vector3,datafile){
  
  estalpha=vector3[1]
  estbeta=vector3[2]
  estsigma=vector3[3]
  
  x=datafile
  
  sum2=0
  large_N=4000
  
  for(i in 1:large_N)
  {
    
    
    pdf_term1=estalpha/estsigma
    #print(pdf_term1)
    pdf_term2=(log((i/estsigma) +1)+1+estbeta)/((i/estsigma)+1)
    #print(pdf_term2)
    pdf_term3=((log((i/estsigma)+1))^(estbeta))/((log((i/estsigma)+1)+1)^(estbeta+1))
    #print(pdf_term3)
    pdf_term4=exp(-estalpha*(((log((i/estsigma)+1))^(estbeta+1))/((log((i/estsigma)+1)+1)^(estbeta))))
    #print(pdf_term4)
    
    sum2 = sum2 + pdf_term1*pdf_term2*pdf_term3*pdf_term4
    #print(sum2)
    
    #print("---------")
  }
  
  #print('%%%%%%%%%')
  
  C2=1/sum2
  #print(sum2)
  #print(C2)
  
  pdf_term5=estalpha/estsigma
  
  pdf_term6=(log((x/estsigma)+1)+1+estbeta)/((x/estsigma)+1)
  
  pdf_term7=((log((x/estsigma)+1))^(estbeta))/((log((x/estsigma)+1)+1)^(estbeta+1))
  
  pdf_term8=exp(-estalpha*(((log((x/estsigma)+1))^(estbeta+1))/((log((x/estsigma)+1)+1)^(estbeta))))
  
  final_pdf_term=C2*pdf_term5*pdf_term6*pdf_term7*pdf_term8
  
  return(final_pdf_term)
  
  #return(probability_array)
  
}


outfun<-probability.fun_version4(c(output11$par[1],output11$par[2],output11$par[3]),y)
#print('probability_sum')
#print(sum(outfun))
#plot(y,outfun)


finaloutput<-outfun*var1
#print(finaloutput)


write.csv(finaloutput,'output_twitter_version4.csv')




##################################### Power-Law Distribution ########################

rm(list = ls())


data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata



x_min=2

var1<-sum(freq[x_min:length(freq)])

## continuous case
estalpha<-1+var1*((freq%*%log(y/x_min))^(-1))
#altestaplpha<-1+var1*((sum(log(rep(y1,freq1))))^(-1))


## discrete case
#estalpha<-1+var1*((freq%*%log(y/(x_min-0.5)))^(-1))

#estalpha<-2
#estalpha<-1.5



probability.fun<- function(gamma,datafile){
  
  estalpha1<-gamma
  
  x<-datafile
  
  large_N=4000
  
  c1<-1/sum(1/((x_min:large_N)^(estalpha1)))
  
  h1<-c1/x^(estalpha1)
  
  #h1<-((estalpha1-1)/x_min)*((x/x_min)^(-estalpha1))
  
  return(h1)
  
}


outfunction<-probability.fun(estalpha,x_min:max(y1))

finaloutput<-outfunction*var1


write.csv(finaloutput,'output_twitter_poweralw.csv')


#########################################  Pareto Type-I #############################

rm(list = ls())



data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)
#x_min=1

estxmin1<-1

n<-sum(freq)


estalpha1<-n/(t(freq)%*%(log(y)-log(estxmin1)))


probability.fun<- function(gamma,datafile){
  
  estxmin<-gamma[1]
  
  estalpha<-gamma[2]
  
  x<-datafile
  
  
  large_N=4000
  c1=1/(sum(1/((1:large_N)^(estalpha+1)))*(estalpha*(estxmin^estalpha)))
  
  term1<-estalpha
  
  term2<-estxmin^estalpha
  
  term3<-x^(estalpha+1)
  
  
  term4<-term2/term3
  
  #term5<-term1*term4
  term5<-term1*term4*c1
  
  
  return(term5)
  
}


outfunction<-probability.fun(c(estxmin1,estalpha1),1:max(y1))

finaloutput<-outfunction*var1


write.csv(finaloutput,'output_twitter_pareto.csv')



########################### Log-Normal Distribution ######################

rm(list = ls())



data<-read.csv('indegree_frequency_twitter.csv')



y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)
n=var1


########  estimated parameters  ###

estmu=(t(freq)%*%log(y))/n

estsigma=sqrt((t(freq)%*%((log(y)-estmu)^2))/n)



probability.fun<- function(gamma,datafile){
  
  estmu1<-gamma[1]
  
  estsigma1<-gamma[2]
  
  xmin1<-1
  
  x<-datafile
  

  
  h3<-1/(sqrt(2*pi)*(estsigma))
  
  #print(h3)
  
  large_N=4000
  c1<-1/(h3*sum(exp(-((log(1:large_N)-estmu1)/(sqrt(2)*estsigma1))^2)/(1:large_N)))
  
  
  
  h4<-x
  
  #print(h4)
  
  h5<-exp(-((log(x)-estmu1)/(sqrt(2)*estsigma1))^2)
  
  h6<-h5/h4
  
  #print(h6)
  
  #h7<-h3*h6
  h7<-h3*h6*c1
  
  #print(h7)
  
  return(h7)
  
}

outfunction<-probability.fun(c(estmu,estsigma),1:max(y1))


finaloutput<-outfunction*var1


write.csv(finaloutput,'output_twitter_lognormal.csv')


############################# PowerLaw-Cutoff Distribution ##########################

rm(list = ls())



data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)

x_min=1



powerlawcutof.lik<-function(vector1,freq){
  alpha<-vector1[1]
  lamda<-vector1[2]
  xmin<-1
  
  n<-sum(freq)
  sum=0
  
  large_N=4000
  
  for(i in 1:large_N)
  {
    sum = sum + ((1/(i^alpha))*(1/exp(lamda*i)))
  }
  
  C=1/sum
  
  ##############################
  
  k2<-n*log(C)
  
  
  k4<-alpha*(t(freq)%*%log(y))
  
  k5<-lamda*(t(y)%*%freq)
  
  k6<-sum(c(k2,-k4,-k5))
  
  #print(k6)
  
  return(-k6)
  
}


output11<-optim(c(1,1),powerlawcutof.lik,freq=freq)

#print(output11)


probability.fun<- function(vector2,datafile){
  
  estalpha<-vector2[1]
  
  estlamda<-vector2[2]
  
  xmin1<-1
  
  x<-datafile
  
  sum1=0
  
  large_N=4000
  
  for(i in 1:large_N)
  {
    sum1 = sum1 + ((1/(i^estalpha))*(1/exp(estlamda*i)))
  }
  
  C1=1/sum1
  
  #####################################
  
  j3<-x^(-estalpha)
  
  j4<-exp(-estlamda*x)
  
  j5<-(C1*j3*j4)
  
  return(j5)
  
}


outfun<-probability.fun(c(output11$par[1],output11$par[2]),y)

finaloutput<-outfun*var1



write.csv(finaloutput,'output_twitter_powerlaw_cutoff.csv')




#################################### Exponentiaal distribution ######################

rm(list = ls())



data<-read.csv('indegree_frequency_twitter.csv')


y1=data$C1
freq1<-data$C2

fulldata<-numeric(max(y1))

for(k in 1:max(y1))
{
  if(is.element(k,y1)==TRUE)
  {
    index=match(k,y1)
    fulldata[k]=freq1[index]
  }
  else
    fulldata[k]=0
}

y<-1:max(y1)
freq<-fulldata

var1<-sum(freq)

x_min=1


estlambda=(sum(freq))/(t(freq)%*%y)


probability.fun<- function(gamma,datafile){
  
  estlambda<-gamma
  
  xmin1<-1
  
  x<-datafile
  
  sum1=0
  
  large_N=4000
  
  
  for(i in 1:large_N)
  {
    sum1 = sum1 + (exp(-estlambda*i)*(estlambda))
  }
  
  C1=1/sum1
  
  h2<-((estlambda)*(exp(-x*estlambda)))
  
  h3<-C1*h2
  
  return(h3)
  
}

#outfunction<-probability.fun(output$par,y)
outfunction<-probability.fun(estlambda,y)


finaloutput<-outfunction*var1

write.csv(finaloutput,'output_twitter_exponential.csv')

##############################################################################################
