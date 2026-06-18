invgamma<-function(xm, se){
  a=0
  b=0
  c=0
  
  a=(xm^2)/(se^2)
  b=(se^2)/xm
  c<-rgamma(1,a,b)
}
c