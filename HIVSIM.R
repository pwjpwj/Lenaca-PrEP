HIVSIM<-function(NoSem=52, NumHSH=808862-50478, NumHSHVIH=50478, 
                 ACD4=0.065, I8CD4=0.107, I7CD4=0.179, Rest=0.65, 
                 G=0.001/52, Woff=0/52, W=0/52, m=0.001/52, 
                 Cvh=52.5/52, Ch=9/52, Cl=2.5/52, Cvl=0.2/52, 
                 nvh=0.1, nh=0.40, nl=0.30, nvl=0.20, B=0.024/52,
                 I12=1-exp(-0.25), I23=1-exp(-0.25), I34=1-exp(-0.25), 
                 I45=1-exp(-0.25*3), I55=1-exp(-0.25*3), 
                 I56=1-exp(-0.25*3), I67=1-exp(-0.25*3), 
                 I78=1-exp(-0.25*3), I89=1-exp(-0.25*3)){
  
  # Compartments not on PrEP
  # Very high risk
  S0vh<-matrix(nrow=NoSem, ncol=85) 
  S0vh[1,1]<-nvh*NumHSH
  I1vh<-matrix(nrow=NoSem, ncol=85)
  I1vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I2vh<-matrix(nrow=NoSem, ncol=85)
  I2vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I3vh<-matrix(nrow=NoSem, ncol=85)
  I3vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I4vh<-matrix(nrow=NoSem, ncol=85)
  I4vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I5vh<-matrix(nrow=NoSem, ncol=85)
  I5vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I6vh<-matrix(nrow=NoSem, ncol=85)
  I6vh[1,1]<-NumHSHVIH*nvh*Rest/6
  I7vh<-matrix(nrow=NoSem, ncol=85)
  I7vh[1,1]<-NumHSHVIH*I7CD4*nvh
  I8vh<-matrix(nrow=NoSem, ncol=85)
  I8vh[1,1]<-NumHSHVIH*I8CD4*nvh
  AIDSvh<-matrix(nrow=NoSem, ncol=85)
  AIDSvh[1,1]<-NumHSHVIH*ACD4*nvh
  
  #high risk
  S0h<-matrix(nrow=NoSem, ncol=85) 
  S0h[1,1]<-nh*NumHSH
  I1h<-matrix(nrow=NoSem, ncol=85)
  I1h[1,1]<-NumHSHVIH*nh*Rest/6
  I2h<-matrix(nrow=NoSem, ncol=85)
  I2h[1,1]<-NumHSHVIH*nh*Rest/6
  I3h<-matrix(nrow=NoSem, ncol=85)
  I3h[1,1]<-NumHSHVIH*nh*Rest/6
  I4h<-matrix(nrow=NoSem, ncol=85)
  I4h[1,1]<-NumHSHVIH*nh*Rest/6
  I5h<-matrix(nrow=NoSem, ncol=85)
  I5h[1,1]<-NumHSHVIH*nh*Rest/6
  I6h<-matrix(nrow=NoSem, ncol=85)
  I6h[1,1]<-NumHSHVIH*nh*Rest/6
  I7h<-matrix(nrow=NoSem, ncol=85)
  I7h[1,1]<-NumHSHVIH*nh*I7CD4
  I8h<-matrix(nrow=NoSem, ncol=85)
  I6h[1,1]<-NumHSHVIH*nh*I8CD4
  AIDSh<-matrix(nrow=NoSem, ncol=85)
  AIDSh[1,1]<-NumHSHVIH*nh*ACD4
  #low risk
  S0l<-matrix(nrow=NoSem, ncol=85)
  S0l[1,1]<-nl*NumHSH
  I1l<-matrix(nrow=NoSem, ncol=85)
  I1l[1,1]<-NumHSHVIH*nl*Rest/6
  I2l<-matrix(nrow=NoSem, ncol=85)
  I2l[1,1]<-NumHSHVIH*nl*Rest/6
  I3l<-matrix(nrow=NoSem, ncol=85)
  I3l[1,1]<-NumHSHVIH*nl*Rest/6
  I4l<-matrix(nrow=NoSem, ncol=85)
  I4l[1,1]<-NumHSHVIH*nl*Rest/6
  I5l<-matrix(nrow=NoSem, ncol=85)
  I5l[1,1]<-NumHSHVIH*nl*Rest/6
  I6l<-matrix(nrow=NoSem, ncol=85)
  I6l[1,1]<-NumHSHVIH*nl*Rest/6
  I7l<-matrix(nrow=NoSem, ncol=85)
  I7l[1,1]<-NumHSHVIH*nl*I7CD4
  I8l<-matrix(nrow=NoSem, ncol=85)
  I8l[1,1]<-NumHSHVIH*nl*I8CD4
  AIDSl<-matrix(nrow=NoSem, ncol=85)
  AIDSl[1,1]<-NumHSHVIH*nl*ACD4
  #very low risk
  S0vl<-matrix(nrow=NoSem, ncol=85) 
  S0vl[1,1]<-nvl*NumHSH
  I1vl<-matrix(nrow=NoSem, ncol=85)
  I1vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I2vl<-matrix(nrow=NoSem, ncol=85)
  I2vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I3vl<-matrix(nrow=NoSem, ncol=85)
  I3vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I4vl<-matrix(nrow=NoSem, ncol=85)
  I4vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I5vl<-matrix(nrow=NoSem, ncol=85)
  I5vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I6vl<-matrix(nrow=NoSem, ncol=85)
  I6vl[1,1]<-NumHSHVIH*nvl*Rest/6
  I7vl<-matrix(nrow=NoSem, ncol=85)
  I7vl[1,1]<-NumHSHVIH*nvl*I7CD4
  I8vl<-matrix(nrow=NoSem, ncol=85)
  I8vl[1,1]<-NumHSHVIH*nvl*I8CD4
  AIDSvl<-matrix(nrow=NoSem, ncol=85)
  AIDSvl[1,1]<-NumHSHVIH*nvl*ACD4
  
  #Compartments on PrEP
  #Very high risk
  S0pvh<-matrix(nrow=NoSem, ncol=85)
  I1pvh<-matrix(nrow=NoSem, ncol=85)
  I2pvh<-matrix(nrow=NoSem, ncol=85)
  I3pvh<-matrix(nrow=NoSem, ncol=85)
  I4pvh<-matrix(nrow=NoSem, ncol=85)
  I5pvh<-matrix(nrow=NoSem, ncol=85)
  I6pvh<-matrix(nrow=NoSem, ncol=85)
  I7pvh<-matrix(nrow=NoSem, ncol=85)
  I8pvh<-matrix(nrow=NoSem, ncol=85)
  AIDSpvh<-matrix(nrow=NoSem, ncol=85)
  #high risk
  S0ph<-matrix(nrow=NoSem, ncol=85)
  I1ph<-matrix(nrow=NoSem, ncol=85)
  I2ph<-matrix(nrow=NoSem, ncol=85)
  I3ph<-matrix(nrow=NoSem, ncol=85)
  I4ph<-matrix(nrow=NoSem, ncol=85)
  I5ph<-matrix(nrow=NoSem, ncol=85)
  I6ph<-matrix(nrow=NoSem, ncol=85)
  I7ph<-matrix(nrow=NoSem, ncol=85)
  I8ph<-matrix(nrow=NoSem, ncol=85)
  AIDSph<-matrix(nrow=NoSem, ncol=85)
  
  #Treated
  Tr<-matrix(nrow=NoSem, ncol=85)
  
  
  ## The equations
  
  for (j in 1:85)
  for (i in 2:NoSem){
    #Calculating the time changing variables
    
    #Prevalence of infection in very high risk group
    PrevInfvh<-sum(I1vh[i-1], I2vh[i-1], I3vh[i-1], I4vh[i-1], I5vh[i-1], I6vh[i-1], I7vh[i-1], I8vh[i-1], AIDSvh[i-1], I1pvh[i-1], I2pvh[i-1], I3pvh[i-1], I4pvh[i-1], I5pvh[i-1], I6pvh[i-1], I7pvh[i-1], I8pvh[i-1], AIDSpvh[i-1])/sum(S0vh[i-1], I1vh[i-1], I2vh[i-1], I3vh[i-1], I4vh[i-1], I5vh[i-1], I6vh[i-1], I7vh[i-1], I8vh[i-1], AIDSvh[i-1], S0pvh[i-1], I1pvh[i-1], I2pvh[i-1], I3pvh[i-1], I4pvh[i-1], I5pvh[i-1], I6pvh[i-1], I7pvh[i-1], I8pvh[i-1], AIDSpvh[i-1])
    #Prevalence of infection in high risk group
    PrevInfh<-sum(I1h[i-1], I2h[i-1], I3h[i-1], I4h[i-1], I5h[i-1], I6h[i-1], I7h[i-1], I8h[i-1], AIDSh[i-1], I1ph[i-1], I2ph[i-1], I3ph[i-1], I4ph[i-1], I5ph[i-1], I6ph[i-1], I7ph[i-1], I8ph[i-1], AIDSph[i-1])/sum(S0h[i-1], I1h[i-1], I2h[i-1], I3h[i-1], I4h[i-1], I5h[i-1], I6h[i-1], I7h[i-1], I8h[i-1], AIDSh[i-1], S0ph[i-1], I1ph[i-1], I2ph[i-1], I3ph[i-1], I4ph[i-1], I5ph[i-1], I6ph[i-1], I7ph[i-1], I8ph[i-1], AIDSph[i-1])
    #Prevalence of infection in low risk group
    PrevInfl<-sum(I1l[i-1], I2l[i-1], I3l[i-1], I4l[i-1], I5l[i-1], I6l[i-1], I7l[i-1], I8l[i-1], AIDSl[i-1])/sum(S0l[i-1], I1l[i-1], I2l[i-1], I3l[i-1], I4l[i-1], I5l[i-1], I6l[i-1], I7l[i-1], I8l[i-1], AIDSl[i-1])
    
    #Prevalence of infection in very low risk group
    PrevInfvl<-sum(I1vl[i-1], I2vl[i-1], I3vl[i-1], I4vl[i-1], I5vl[i-1], I6vl[i-1], I7vl[i-1], I8vl[i-1], AIDSvl[i-1])/sum(S0vl[i-1], I1vl[i-1], I2vl[i-1], I3vl[i-1], I4vl[i-1], I5vl[i-1], I6vl[i-1], I7vl[i-1], I8vl[i-1], AIDSvl[i-1])
    
    #Probability that a partner will be in a specific group
    Gvh<-(Cvh*nvh)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gh<-(Ch*nh)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gl<-(Cl*nl)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gvl<-(Cvl*nvl)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    
    #Probability that the partner is infectious
    P<-Gvh*PrevInfvh+Gh*PrevInfh+Gl*PrevInfl+Gvl*PrevInfvl
    
    #Force of infection
    Lvh<-Cvh*B*P
    Lh<-Ch*B*P
    Ll<-Cl*B*P
    Lvl<-Cvl*B*P
    
    #Differential Equations each time step is one week
    #Not on PrEP, very high risk group
    S0vh[i]<-S0vh[i-1]+S0vh[i]*Gvh+S0pvh[i]*Woff-(I1vh[i-1]+I2vh[i-1]+I3vh[i-1]+I4vh[i-1]+I5vh[i-1]+I6vh[i-1]+I7vh[i-1]+I8vh[i-1]+AIDSvh[i-1])*Lvh-S0vh[i-1]*(W+m)
    I1vh[i]<-I1vh[i-1]+(I1vh[i-1]+I2vh[i-1]+I3vh[i-1]+I4vh[i-1]+I5vh[i-1]+I6vh[i-1]+I7vh[i-1]+I8vh[i-1]+AIDSvh[i-1])*Lvh-I1vh[i-1]*(m)-I1vh[i-1]*I12
    I2vh[i]<-I1vh[i-1]*I12-I2vh[i-1]*I23-I2vh[i-1]*(m)
    I3vh[i]<-I2vh[i-1]*I23-I3vh[i-1]*I34-I3vh[i-1]*(m)
    I4vh[i]<-I3vh[i-1]*I34-I4vh[i-1]*I45-I4vh[i-1]*(m)   
    I5vh[i]<-I4vh[i-1]*I45-I5vh[i-1]*I56-I5vh[i-1]*(m)
    I6vh[i]<-I5vh[i-1]*I56-I6vh[i-1]*I67-I6vh[i-1]*(m)
    I7vh[i]<-I6vh[i-1]*I67-I7vh[i-1]*I78-I7vh[i-1]*(m)
    I8vh[i]<-I7vh[i-1]*I78-I8vh[i-1]*I89-I8vh[i-1]*(m)
    AIDSvh[i]<-I8vh[i-1]*I89-AIDSvh[i-1]*(m)
    
    #Not on PrEP, high risk group
    S0h[i]<-S0h[i-1]+S0h[i]*Gh+S0ph[i]*Woff-(I1h[i-1]+I2h[i-1]+I3h[i-1]+I4h[i-1]+I5h[i-1]+I6h[i-1]+I7h[i-1]+I8h[i-1]+AIDSh[i-1])*Lh-S0h[i-1]*(W+m)
    I1h[i]<-I1h[i-1]+(I1h[i-1]+I2h[i-1]+I3h[i-1]+I4h[i-1]+I5h[i-1]+I6h[i-1]+I7h[i-1]+I8h[i-1]+AIDSh[i-1])*Lh-I1h[i-1]*(m)-I1h[i-1]*I12
    I2h[i]<-I1h[i-1]*I12-I2h[i-1]*I23-I2h[i-1]*(m)
    I3h[i]<-I2h[i-1]*I23-I3h[i-1]*I34-I3h[i-1]*(m)
    I4h[i]<-I3h[i-1]*I34-I4h[i-1]*I45-I4h[i-1]*(m)   
    I5h[i]<-I4h[i-1]*I45-I5h[i-1]*I56-I5h[i-1]*(m)
    I6h[i]<-I5h[i-1]*I56-I6h[i-1]*I67-I6h[i-1]*(m)
    I7h[i]<-I6h[i-1]*I67-I7h[i-1]*I78-I7h[i-1]*(m)
    I8h[i]<-I7h[i-1]*I78-I8h[i-1]*I89-I8h[i-1]*(m)
    AIDSh[i]<-I8h[i-1]*I89-AIDSh[i-1]*(m)
    
    #Not on PrEP, Low risk group
    S0l[i]<-S0l[i-1]+S0l[i]*Gl+(I1l[i-1]+I2l[i-1]+I3l[i-1]+I4l[i-1]+I5l[i-1]+I6l[i-1]+I7l[i-1]+I8l[i-1]+AIDSl[i-1])*Ll-S0l[i-1]*(W+m)
    I1l[i]<-I1l[i-1]+(I1l[i-1]+I2l[i-1]+I3l[i-1]+I4l[i-1]+I5l[i-1]+I6l[i-1]+I7l[i-1]+I8l[i-1]+AIDSl[i-1])*Ll-I1l[i-1]*(m)-I1l[i-1]*I12
    I2l[i]<-I1l[i-1]*I12-I2l[i-1]*I23-I2l[i-1]*(m)
    I3l[i]<-I2l[i-1]*I23-I3l[i-1]*I34-I3l[i-1]*(m)
    I4l[i]<-I3l[i-1]*I34-I4l[i-1]*I45-I4l[i-1]*(m)   
    I5l[i]<-I4l[i-1]*I45-I5l[i-1]*I56-I5l[i-1]*(m)
    I6l[i]<-I5l[i-1]*I56-I6l[i-1]*I67-I6l[i-1]*(m)
    I7l[i]<-I6l[i-1]*I67-I7l[i-1]*I78-I7l[i-1]*(m)
    I8l[i]<-I7l[i-1]*I78-I8l[i-1]*I89-I8l[i-1]*(m)
    AIDSl[i]<-I8l[i-1]*I89-AIDSl[i-1]*(m)
    
    #Not on PrEP, very low risk group
    S0vl[i]<-S0vl[i-1]+S0vl[i]*G+(I1vl[i-1]+I2vl[i-1]+I3vl[i-1]+I4vl[i-1]+I5vl[i-1]+I6vl[i-1]+I7vl[i-1]+I8vl[i-1]+AIDSvl[i-1])*Lvl-S0vl[i-1]*(W+m)
    I1vl[i]<-I1vl[i-1]+(I1vl[i-1]+I2vl[i-1]+I3vl[i-1]+I4vl[i-1]+I5vl[i-1]+I6vl[i-1]+I7vl[i-1]+I8vl[i-1]+AIDSvl[i-1])*Lvl-I1vl[i-1]*(m)-I1vl[i-1]*I12
    I2vl[i]<-I1vl[i-1]*I12-I2vl[i-1]*I23-I2vl[i-1]*(m)
    I3vl[i]<-I2vl[i-1]*I23-I3vl[i-1]*I34-I3vl[i-1]*(m)
    I4vl[i]<-I3vl[i-1]*I34-I4vl[i-1]*I45-I4vl[i-1]*(m)   
    I5vl[i]<-I4vl[i-1]*I45-I5vl[i-1]*I56-I5vl[i-1]*(m)
    I6vl[i]<-I5vl[i-1]*I56-I6vl[i-1]*I67-I6vl[i-1]*(m)
    I7vl[i]<-I6vl[i-1]*I67-I7vl[i-1]*I78-I7vl[i-1]*(m)
    I8vl[i]<-I7vl[i-1]*I78-I8vl[i-1]*I89-I8vl[i-1]*(m)
    AIDSvl[i]<-I8vl[i-1]*I89-AIDSvl[i-1]*(m)
    
    #On PrEP, very high risk group
    S0pvh[i]<-S0pvh[i-1]+S0pvh[i]*Gvh+S0vh[i]*Woff-(I1pvh[i-1]+I2pvh[i-1]+I3pvh[i-1]+I4pvh[i-1]+I5pvh[i-1]+I6pvh[i-1]+I7pvh[i-1]+I8pvh[i-1]+AIDSpvh[i-1])*Lvh-S0pvh[i-1]*(W+m)
    I1pvh[i]<-I1pvh[i-1]+(I1pvh[i-1]+I2pvh[i-1]+I3pvh[i-1]+I4pvh[i-1]+I5pvh[i-1]+I6pvh[i-1]+I7pvh[i-1]+I8pvh[i-1]+AIDSpvh[i-1])*Lvh-I1pvh[i-1]*(m)-I1pvh[i-1]*I12
    I2pvh[i]<-I1pvh[i-1]*I12-I2pvh[i-1]*I23-I2pvh[i-1]*(m)
    I3pvh[i]<-I2pvh[i-1]*I23-I3pvh[i-1]*I34-I3pvh[i-1]*(m)
    I4pvh[i]<-I3pvh[i-1]*I34-I4pvh[i-1]*I45-I4pvh[i-1]*(m)   
    I5pvh[i]<-I4pvh[i-1]*I45-I5pvh[i-1]*I56-I5pvh[i-1]*(m)
    I6pvh[i]<-I5pvh[i-1]*I56-I6pvh[i-1]*I67-I6pvh[i-1]*(m)
    I7pvh[i]<-I6pvh[i-1]*I67-I7pvh[i-1]*I78-I7pvh[i-1]*(m)
    I8pvh[i]<-I7pvh[i-1]*I78-I8pvh[i-1]*I89-I8pvh[i-1]*(m)
    AIDSpvh[i]<-I8pvh[i-1]*I89-AIDSpvh[i-1]*(m)  
    #On PrEP, high risk group
    S0ph[i]<-S0ph[i-1]+S0ph[i]*Gh+S0h[i]*Woff-(I1ph[i-1]+I2ph[i-1]+I3ph[i-1]+I4ph[i-1]+I5ph[i-1]+I6ph[i-1]+I7ph[i-1]+I8ph[i-1]+AIDSph[i-1])*Lh-S0ph[i-1]*(W+m)
    I1ph[i]<-I1ph[i-1]+(I1ph[i-1]+I2ph[i-1]+I3ph[i-1]+I4ph[i-1]+I5ph[i-1]+I6ph[i-1]+I7ph[i-1]+I8ph[i-1]+AIDSph[i-1])*Lh-I1ph[i-1]*(m)-I1ph[i-1]*I12
    I2ph[i]<-I1ph[i-1]*I12-I2ph[i-1]*I23-I2ph[i-1]*(m)
    I3ph[i]<-I2ph[i-1]*I23-I3ph[i-1]*I34-I3ph[i-1]*(m)
    I4ph[i]<-I3ph[i-1]*I34-I4ph[i-1]*I45-I4ph[i-1]*(m)   
    I5ph[i]<-I4ph[i-1]*I45-I5ph[i-1]*I56-I5ph[i-1]*(m)
    I6ph[i]<-I5ph[i-1]*I56-I6ph[i-1]*I67-I6ph[i-1]*(m)
    I7ph[i]<-I6ph[i-1]*I67-I7ph[i-1]*I78-I7ph[i-1]*(m)
    I8ph[i]<-I7ph[i-1]*I78-I8ph[i-1]*I89-I8ph[i-1]*(m)
    AIDSph[i]<-I8ph[i-1]*I89-AIDSph[i-1]*(m)  
  }
Q<-as.data.frame(cbind(1:NoSem,S0vh,S0h,S0l,S0vl,S0pvh,S0ph,I1vh, I1h, I1l, I1vl, I1pvh,I1ph,I2vh, I2h, 
                       I2l, I2vl, I2pvh,I2ph,I3vh, I3h, I3l, I3vl, I3pvh,I3ph, I4vh, I4h, I4l, 
                       I4vl, I4pvh, I4ph,I5vh, I5h, I5l, I5vl, I5pvh, I5ph,I6vh, I6h, I6l, I6vl, 
                       I6pvh,I6ph,I7vh, I7h, I7l, I7vl, I7pvh, I7ph, I8vh, I8h, I8l, I8vl, I8pvh, 
                       I8ph,AIDSvh, AIDSh, AIDSl, AIDSvl, AIDSpvh,AIDSph))
}