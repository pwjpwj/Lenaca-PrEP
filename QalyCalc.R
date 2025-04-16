QalyCalc<-function(x,r=0.03,w=w){
  Qaly<-0
  j<-0
  a<-0
  for (i in 1:length(x)){
    a<-a+1
    if (a==52) j=j+1
    if(a==52) a=0
    Qaly<-Qaly+(x[i]*w)/(1+r)^j
    
    }
  Qaly
}