CostCalc<-function(x,r=0.03,C){
  Cost<-0
  j<-0
  a<-0
  for (i in 1:length(x)){
    a<-a+1
    if (a==52) j=j+1
    if(a==52) a=0
    Cost<-Cost+x[i]*C/(1+r)^j}
  Cost

}