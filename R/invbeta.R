invbeta<-function(xm=0,se=0){
a=0
b=0
c=0

a<-xm*((xm*(1-xm)/(se^2))-1)
b<-(xm*(1-xm)/(se^2))-1-a
c<-rbeta(1,a,b)}
c