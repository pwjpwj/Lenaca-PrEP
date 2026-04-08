#
# HIV transmission model in Spain
#
# Loading packages
library(ggplot2)

# Compartments not on PrEP
S0<-vector(mode="numeric", length=105) 
S0[1]<-0.5*100000
I1<-vector(mode="numeric", length=105) 
I1[1]<-0.5*100000
I2<-vector(mode="numeric", length=105)
I3<-vector(mode="numeric", length=105)
I4<-vector(mode="numeric", length=105)
I5<-vector(mode="numeric", length=105)
I6<-vector(mode="numeric", length=105)
I7<-vector(mode="numeric", length=105)
I8<-vector(mode="numeric", length=105)
AIDS<-vector(mode="numeric", length=105)

#Compartments on PrEP
S0p<-vector(mode="numeric", length=105)
I1p<-vector(mode="numeric", length=105)
I2p<-vector(mode="numeric", length=105)
I3p<-vector(mode="numeric", length=105)
I4p<-vector(mode="numeric", length=105)
I5p<-vector(mode="numeric", length=105)
I6p<-vector(mode="numeric", length=105)
I7p<-vector(mode="numeric", length=105)
I8p<-vector(mode="numeric", length=105)
AIDSp<-vector(mode="numeric", length=105)

#Treated
Tr<-vector(mode="numeric", length=100)

#Values 
G<-0.001/52  #Population growth rate 1/1000 anually
Woff<-0.3/52 #Proportion of people who go off PrEP 30% anually
W<-0.3/52    #Proportion of people who go on PrEP 30% anually
L<-0.05/52   #Force of infection
m<-0.001/52  #Mortality general population
I12<-1-exp(-0.25) #1 month duration in stage 1 of HIV infection
I23<-1-exp(-0.25) #1 month duration in stage 2 of HIV infection
I34<-1-exp(-0.25) #1 month duration in stage 3 of HIV infection
I45<-1-exp(-0.25*3) #3 month duration in stage 4 of HIV infection
I55<-1-exp(-0.25*3) #3 month duration in stage 5 of HIV infection
I56<-1-exp(-0.25*3) #3 month duration in stage 6 of HIV infection

#Equations each time step is one week

for (i in 2:105){
  S0[i]<-S0[i-1]+S0[i]*G+S0p[i]*Woff-(I1[i-1]+I2[i-1]+I3[i-1]+I4[i-1]+I5[i-1]
  +I6[i-1]+I7[i-1]+I8[i-1]+AIDS[i-1])*L-S0[i-1]*(W+m)
  
  I1[i]<-I1[i-1]+(I1[i-1]+I2[i-1]+I3[i-1]+I4[i-1]+I5[i-1]
                  +I6[i-1]+I7[i-1]+I8[i-1]+AIDS[i-1])*L
  -I1[i-1]*(m)-I1[i-1]*I12
  
  I2[i]<-I1[i-1]*I12-I2[i-1]*I23-I2[i-1]*(m)
  I3[i]<-I2[i-1]*I23-I3[i-1]*I34-I3[i-1]*(m)
  I4[i]<-I3[i-1]*I34-I4[i-1]*I45-I4[i-1]*(m)   
  I5[i]<-I4[i-1]*I45-I5[i-1]*I55-I5[i-1]*(m) 
}


Q<-as.data.frame(cbind(1:105,S0,I1,I2))

plot1<-ggplot(Q, aes(V1))+
    geom_line(aes(y=S0, color="black"))+
    geom_line(aes(y=I1, color="blue"))+
    geom_line(aes(y=I2, color="red"))
plot1


