HIVSIMPREP<-function(Age=get("Age", globalenv()) , NoSem=get("NoSem", globalenv()),
                 NumHSH=get("NumHSH", globalenv()), NumHSHVIH=get("NumHSHVIH", globalenv()), 
                 ACD4=get("ACD4", globalenv()), I8CD4=get("I8CD4", globalenv()),
                 I7CD4=get("I7CD4", globalenv()), Rest=get("Rest", globalenv()),
                 G=get("G", globalenv()), Woff=get("Woff", globalenv()), 
                 W=get("W", globalenv()), m=get("m", globalenv()),
                 mVIH235=get("mVIH235", globalenv()), mVIH35=get("mVIH35", globalenv()),
                 mAIDS=get("mAIDS", globalenv()),
                 mTr=get("mTr", globalenv()), Cvh=get("Cvh", globalenv()), 
                 Ch=get("Ch", globalenv()), Cl=get("Cl", globalenv()), 
                 Cvl=get("Cvl", globalenv()), 
                 nvh=get("nvh", globalenv()), nh=get("nvh", globalenv()), 
                 nl=get("nl", globalenv()), nvl=get("nvl", globalenv()), 
                 B=get("B", globalenv()), I12=get("I12", globalenv()), 
                 I23=get("I23", globalenv()), I34=get("I34", globalenv()), 
                 I45=get("I45", globalenv()), I55=get("I55", globalenv()), 
                 I56=get("I56", globalenv()), I67=get("I67", globalenv()), 
                 I78=get("I78", globalenv()), I89=get("I89", globalenv()),
                 PTr=get("PTr", globalenv()), years=get("years", globalenv()),
                 PropTR=get("PropTR", globalenv()),
                 RiskRed=get("RiskRed", globalenv())){
  
  #The list to store the yearly matrices
  Q<-vector("list", years)
    # Compartments not on PrEP s
    # Very high risk
    S0vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8vh<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #high risk
    S0h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7h<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8h<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSh<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #low risk
    S0l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7l<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8l<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSl<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #very low risk
    S0vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8vl<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSvl<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #Treated
    Tr<-matrix(data=0, nrow=NoSem, ncol=Age)
    #Dead
    Dead<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #Force of infection
    Lvh_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Lh_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Ll_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Lvl_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    for (j in 1:Age){
    #Defining starting state for the non prep compartments
    S0vh[1,j]<-nvh*NumHSH*DistEdad[j]
    I1vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I2vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I3vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I4vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I5vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I6vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    I7vh[1,j]<-NumHSHVIH*I7CD4*nvh*DistEdad[j]
    I8vh[1,j]<-NumHSHVIH*I8CD4*nvh*DistEdad[j]
    AIDSvh[1,j]<-NumHSHVIH*nh*ACD4*DistEdad[j]
    
    S0h[1,j]<-nh*NumHSH*DistEdad[j]
    I1h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I2h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I3h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I4h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I5h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I6h[1,j]<-NumHSHVIH*nh*Rest/6*DistEdad[j]
    I7h[1,j]<-NumHSHVIH*nh*I7CD4*DistEdad[j]
    I8h[1,j]<-NumHSHVIH*nh*I8CD4*DistEdad[j]
    AIDSh[1,j]<-NumHSHVIH*nh*ACD4*DistEdad[j]
    
    S0l[1,j]<-nl*NumHSH*DistEdad[j]
    I1l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I2l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I3l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I4l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I5l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I6l[1,j]<-NumHSHVIH*nl*Rest/6*DistEdad[j]
    I7l[1,j]<-NumHSHVIH*nl*I7CD4*DistEdad[j]
    I8l[1,j]<-NumHSHVIH*nl*I8CD4*DistEdad[j]
    AIDSl[1,j]<-NumHSHVIH*nl*ACD4*DistEdad[j]
    
    S0vl[1,j]<-nvl*NumHSH*DistEdad[j]
    I1vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I2vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I3vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I4vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I5vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I6vl[1,j]<-NumHSHVIH*nvl*Rest/6*DistEdad[j]
    I7vl[1,j]<-NumHSHVIH*nvl*I7CD4*DistEdad[j]
    I8vl[1,j]<-NumHSHVIH*nvl*I8CD4*DistEdad[j]
    AIDSvl[1,j]<-NumHSHVIH*nvl*ACD4*DistEdad[j]
    
    Tr[1,j]<-NumHSHVIH*PropTR*DistEdad[j]
    }
    #Compartments on PrEP
    #Very high risk
  
    S0pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8pvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSpvh<-matrix(data=0, nrow=NoSem, ncol=Age)
    #high risk
    S0ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I1ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I2ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I3ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I4ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I5ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I6ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I7ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    I8ph<-matrix(data=0, nrow=NoSem, ncol=Age)
    AIDSph<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    
for (t in 1:years){
  for (j in 1:Age){
    if (t>=2 & j<85){
      S0vh[1,j+1]<-S0vh[52,j]
      I1vh[1,j+1]<-I1vh[52,j]
      I2vh[1,j+1]<-I2vh[52,j]
      I3vh[1,j+1]<-I3vh[52,j]
      I4vh[1,j+1]<-I4vh[52,j]
      I5vh[1,j+1]<-I5vh[52,j]
      I6vh[1,j+1]<-I6vh[52,j]
      I7vh[1,j+1]<-I7vh[52,j]
      I8vh[1,j+1]<-I8vh[52,j]
      AIDSvh[1,j+1]<-AIDSvh[52,j]
      
      S0h[1,j+1]<-S0h[52,j]
      I1h[1,j+1]<-I1h[52,j]
      I2h[1,j+1]<-I2h[52,j]
      I3h[1,j+1]<-I3h[52,j]
      I4h[1,j+1]<-I4h[52,j]
      I5h[1,j+1]<-I5h[52,j]
      I6h[1,j+1]<-I6h[52,j]
      I7h[1,j+1]<-I7h[52,j]
      I8h[1,j+1]<-I8h[52,j]
      AIDSh[1,j+1]<-AIDSh[52,j]
      
      S0l[1,j+1]<-S0l[52,j]
      I1l[1,j+1]<-I1l[52,j]
      I2l[1,j+1]<-I2l[52,j]
      I3l[1,j+1]<-I3l[52,j]
      I4l[1,j+1]<-I4l[52,j]
      I5l[1,j+1]<-I5l[52,j]
      I6l[1,j+1]<-I6l[52,j]
      I7l[1,j+1]<-I7l[52,j]
      I8l[1,j+1]<-I8l[52,j]
      AIDSl[1,j+1]<-AIDSl[52,j]
      
      S0vl[1,j+1]<-S0vl[52,j]
      I1vl[1,j+1]<-I1vl[52,j]
      I2vl[1,j+1]<-I2vl[52,j]
      I3vl[1,j+1]<-I3vl[52,j]
      I4vl[1,j+1]<-I4vl[52,j]
      I5vl[1,j+1]<-I5vl[52,j]
      I6vl[1,j+1]<-I6vl[52,j]
      I7vl[1,j+1]<-I7vl[52,j]
      I8vl[1,j+1]<-I8vl[52,j]
      AIDSvl[1,j+1]<-AIDSvl[52,j]
      
      #Compartments on PrEP
      #Very high risk
      
      S0pvh[1,j+1]<-S0pvh[52,j]
      I1pvh[1,j+1]<-I1pvh[52,j]
      I2pvh[1,j+1]<-I2pvh[52,j]
      I3pvh[1,j+1]<-I3pvh[52,j]
      I4pvh[1,j+1]<-I4pvh[52,j]
      I5pvh[1,j+1]<-I5pvh[52,j]
      I6pvh[1,j+1]<-I6pvh[52,j]
      I7pvh[1,j+1]<-I7pvh[52,j]
      I8pvh[1,j+1]<-I8pvh[52,j]
      AIDSpvh[1,j+1]<-AIDSpvh[52,j]
      
      #high risk
      S0ph[1,j+1]<-S0ph[52,j]
      I1ph[1,j+1]<-I1ph[52,j]
      I2ph[1,j+1]<-I2ph[52,j]
      I3ph[1,j+1]<-I3ph[52,j]
      I4ph[1,j+1]<-I4ph[52,j]
      I5ph[1,j+1]<-I5ph[52,j]
      I6ph[1,j+1]<-I6ph[52,j]
      I7ph[1,j+1]<-I7ph[52,j]
      I8ph[1,j+1]<-I8ph[52,j]
      AIDSph[1,j+1]<-AIDSph[52,j]
      
      #Treated
      Tr[1,j+1]<-Tr[52,j]  
    }
    
    
  for (i in 2:NoSem){
      #Calculating the time changing variables
      
      #Prevalence of infection in very high risk group
      PrevInfvh<-sum(I1vh[i-1,], I2vh[i-1,], I3vh[i-1,], I4vh[i-1,], I5vh[i-1,], I6vh[i-1,], I7vh[i-1,], I8vh[i-1,], AIDSvh[i-1,], I1pvh[i-1,], I2pvh[i-1,], I3pvh[i-1,], I4pvh[i-1,], I5pvh[i-1,], I6pvh[i-1,], I7pvh[i-1,], I8pvh[i-1,], AIDSpvh[i-1,])/sum(S0vh[i-1,], I1vh[i-1,], I2vh[i-1,], I3vh[i-1,], I4vh[i-1,], I5vh[i-1,], I6vh[i-1,], I7vh[i-1,], I8vh[i-1,], AIDSvh[i-1,], S0pvh[i-1,], I1pvh[i-1,], I2pvh[i-1,], I3pvh[i-1,], I4pvh[i-1,], I5pvh[i-1,], I6pvh[i-1,], I7pvh[i-1,], I8pvh[i-1,], AIDSpvh[i-1,])
      #Prevalence of infection in high risk group
      PrevInfh<-sum(I1h[i-1,], I2h[i-1,], I3h[i-1,], I4h[i-1,], I5h[i-1,], I6h[i-1,], I7h[i-1,], I8h[i-1,], AIDSh[i-1,], I1ph[i-1,], I2ph[i-1,], I3ph[i-1,], I4ph[i-1,], I5ph[i-1,], I6ph[i-1,], I7ph[i-1,], I8ph[i-1,], AIDSph[i-1,])/sum(S0h[i-1,], I1h[i-1,], I2h[i-1,], I3h[i-1,], I4h[i-1,], I5h[i-1,], I6h[i-1,], I7h[i-1,], I8h[i-1,], AIDSh[i-1,], S0ph[i-1,], I1ph[i-1,], I2ph[i-1,], I3ph[i-1,], I4ph[i-1,], I5ph[i-1,], I6ph[i-1,], I7ph[i-1,], I8ph[i-1,], AIDSph[i-1,])
      #Prevalence of infection in low risk group
      PrevInfl<-sum(I1l[i-1,], I2l[i-1,], I3l[i-1,], I4l[i-1,], I5l[i-1,], I6l[i-1,], I7l[i-1,], I8l[i-1,], AIDSl[i-1,])/sum(S0l[i-1,], I1l[i-1,], I2l[i-1,], I3l[i-1,], I4l[i-1,], I5l[i-1,], I6l[i-1,], I7l[i-1,], I8l[i-1,], AIDSl[i-1,])
      
      #Prevalence of infection in very low risk group
      PrevInfvl<-sum(I1vl[i-1,], I2vl[i-1,], I3vl[i-1,], I4vl[i-1,], I5vl[i-1,], I6vl[i-1,], I7vl[i-1,], I8vl[i-1,], AIDSvl[i-1,])/sum(S0vl[i-1,], I1vl[i-1,], I2vl[i-1,], I3vl[i-1,], I4vl[i-1,], I5vl[i-1,], I6vl[i-1,], I7vl[i-1,], I8vl[i-1,], AIDSvl[i-1,])
      
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
      S0vh[i,j]<-S0vh[i-1,j]+(S0vh[i-1,j]+S0pvh[i-1,j])*G+S0pvh[i-1,j]*Woff-S0vh[i-1,j]*(Lvh+W+m[j,4])
      I1vh[i,j]<-I1vh[i-1,j]+S0vh[i-1,j]*Lvh-I1vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])-I1vh[i-1,j]*I12
      I2vh[i,j]<-I2vh[i-1,j]+I1vh[i-1,j]*I12-I2vh[i-1,j]*I23-I2vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I3vh[i,j]<-I3vh[i-1,j]+I2vh[i-1,j]*I23-I3vh[i-1,j]*I34-I3vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I4vh[i,j]<-I4vh[i-1,j]+I3vh[i-1,j]*I34-I4vh[i-1,j]*I45-I4vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])   
      I5vh[i,j]<-I5vh[i-1,j]+I4vh[i-1,j]*I45-I5vh[i-1,j]*I56-I5vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I6vh[i,j]<-I6vh[i-1,j]+I5vh[i-1,j]*I56-I6vh[i-1,j]*I67-I6vh[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I7vh[i,j]<-I7vh[i-1,j]+I6vh[i-1,j]*I67-I7vh[i-1,j]*I78-I7vh[i-1,j]*(m[j,5]+mVIH35+PTr[2])
      I8vh[i,j]<-I8vh[i-1,j]+I7vh[i-1,j]*I78-I8vh[i-1,j]*I89-I8vh[i-1,j]*(m[j,5]+mVIH235+PTr[2])
      AIDSvh[i,j]<-AIDSvh[i-1,j]+I8vh[i-1,j]*I89-AIDSvh[i-1,j]*(m[j,5]+mAIDS+PTr[1])
      
      #Not on PrEP, high risk group
  
      S0h[i,j]<-S0h[i-1,j]+(S0h[i-1,j]+S0ph[i-1,j])*G+S0ph[i-1,j]*Woff-S0h[i-1,j]*(Lh+W+m[j,4])
      I1h[i,j]<-I1h[i-1,j]+S0h[i-1,j]*Lh-I1h[i-1,j]*(m[j,5]+mVIH35+PTr[3])-I1h[i-1,j]*I12
      I2h[i,j]<-I2h[i-1,j]+I1h[i-1,j]*I12-I2h[i-1,j]*I23-I2h[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I3h[i,j]<-I3h[i-1,j]+I2h[i-1,j]*I23-I3h[i-1,j]*I34-I3h[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I4h[i,j]<-I4h[i-1,j]+I3h[i-1,j]*I34-I4h[i-1,j]*I45-I4h[i-1,j]*(m[j,5]+mVIH35+PTr[3])   
      I5h[i,j]<-I5h[i-1,j]+I4h[i-1,j]*I45-I5h[i-1,j]*I56-I5h[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I6h[i,j]<-I6h[i-1,j]+I5h[i-1,j]*I56-I6h[i-1,j]*I67-I6h[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I7h[i,j]<-I7h[i-1,j]+I6h[i-1,j]*I67-I7h[i-1,j]*I78-I7h[i-1,j]*(m[j,5]+mVIH35+PTr[2])
      I8h[i,j]<-I8h[i-1,j]+I7h[i-1,j]*I78-I8h[i-1,j]*I89-I8h[i-1,j]*(m[j,5]+mVIH235+PTr[2])
      AIDSh[i,j]<-AIDSh[i-1,j]+I8h[i-1,j]*I89-AIDSh[i-1,j]*(m[j,5]+mAIDS+PTr[1])
      
      #Not on PrEP, Low risk group
  
      S0l[i,j]<-S0l[i-1,j]+S0l[i-1,j]*G-S0l[i-1,j]*Ll-S0l[i-1,j]*(m[j,4])
      I1l[i,j]<-I1l[i-1,j]+S0l[i-1,j]*Ll-I1l[i-1,j]*(m[j,5]+mVIH35+PTr[3])-I1l[i-1,j]*I12
      I2l[i,j]<-I2l[i-1,j]+I1l[i-1,j]*I12-I2l[i-1,j]*I23-I2l[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I3l[i,j]<-I3l[i-1,j]+I2l[i-1,j]*I23-I3l[i-1,j]*I34-I3l[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I4l[i,j]<-I4l[i-1,j]+I3l[i-1,j]*I34-I4l[i-1,j]*I45-I4l[i-1,j]*(m[j,5]+mVIH35+PTr[3])   
      I5l[i,j]<-I5l[i-1,j]+I4l[i-1,j]*I45-I5l[i-1,j]*I56-I5l[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I6l[i,j]<-I6l[i-1,j]+I5l[i-1,j]*I56-I6l[i-1,j]*I67-I6l[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I7l[i,j]<-I7l[i-1,j]+I6l[i-1,j]*I67-I7l[i-1,j]*I78-I7l[i-1,j]*(m[j,5]+mVIH35+PTr[2])
      I8l[i,j]<-I8l[i-1,j]+I7l[i-1,j]*I78-I8l[i-1,j]*I89-I8l[i-1,j]*(m[j,5]+mVIH235+PTr[2])
      AIDSl[i,j]<-AIDSl[i-1,j]+I8l[i-1,j]*I89-AIDSl[i-1,j]*(m[j,5]+mAIDS+PTr[1])
      
      #Not on PrEP, very low risk group
      S0vl[i,j]<-S0vl[i-1,j]+S0vl[i-1,j]*G-S0vl[i-1,j]*Lvl-S0vl[i-1,j]*(m[j,4])
      I1vl[i,j]<-I1vl[i-1,j]+S0vl[i-1,j]*Lvl-I1vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])-I1vl[i-1,j]*I12
      I2vl[i,j]<-I2vl[i-1,j]+I1vl[i-1,j]*I12-I2vl[i-1,j]*I23-I2vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I3vl[i,j]<-I3vl[i-1,j]+I2vl[i-1,j]*I23-I3vl[i-1,j]*I34-I3vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I4vl[i,j]<-I4vl[i-1,j]+I3vl[i-1,j]*I34-I4vl[i-1,j]*I45-I4vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])   
      I5vl[i,j]<-I5vl[i-1,j]+I4vl[i-1,j]*I45-I5vl[i-1,j]*I56-I5vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I6vl[i,j]<-I6vl[i-1,j]+I5vl[i-1,j]*I56-I6vl[i-1,j]*I67-I6vl[i-1,j]*(m[j,5]+mVIH35+PTr[3])
      I7vl[i,j]<-I7vl[i-1,j]+I6vl[i-1,j]*I67-I7vl[i-1,j]*I78-I7vl[i-1,j]*(m[j,5]+mVIH35+PTr[2])
      I8vl[i,j]<-I8vl[i-1,j]+I7vl[i-1,j]*I78-I8vl[i-1,j]*I89-I8vl[i-1,j]*(m[j,5]+mVIH235+PTr[2])
      AIDSvl[i,j]<-AIDSvl[i-1,j]+I8vl[i-1,j]*I89-AIDSvl[i-1,j]*(m[j,5]+mAIDS+PTr[1])
      
      #On PrEP, very high risk group
      PTrp=0.95
  
      S0pvh[i,j]<-S0pvh[i-1,j]+S0vh[i-1,j]*W-S0pvh[i-1,j]*RiskRed*Lvh-S0pvh[i-1,j]*(Woff+m[j,4])
      I1pvh[i,j]<-I1pvh[i-1,j]+S0pvh[i-1,j]*RiskRed*Lvh-I1pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)-I1pvh[i-1,j]*I12
      I2pvh[i,j]<-I2pvh[i-1,j]+I1pvh[i-1,j]*I12-I2pvh[i-1,j]*I23-I2pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I3pvh[i,j]<-I3pvh[i-1,j]+I2pvh[i-1,j]*I23-I3pvh[i-1,j]*I34-I3pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I4pvh[i,j]<-I4pvh[i-1,j]+I3pvh[i-1,j]*I34-I4pvh[i-1,j]*I45-I4pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)   
      I5pvh[i,j]<-I5pvh[i-1,j]+I4pvh[i-1,j]*I45-I5pvh[i-1,j]*I56-I5pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I6pvh[i,j]<-I6pvh[i-1,j]+I5pvh[i-1,j]*I56-I6pvh[i-1,j]*I67-I6pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I7pvh[i,j]<-I7pvh[i-1,j]+I6pvh[i-1,j]*I67-I7pvh[i-1,j]*I78-I7pvh[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I8pvh[i,j]<-I8pvh[i-1,j]+I7pvh[i-1,j]*I78-I8pvh[i-1,j]*I89-I8pvh[i-1,j]*(m[j,5]+mVIH235+PTrp)
      AIDSpvh[i,j]<-AIDSpvh[i-1,j]+I8pvh[i-1,j]*I89-AIDSpvh[i-1,j]*(m[j,5]+mAIDS+PTrp)  
      
      #On PrEP, high risk group
  
      S0ph[i,j]<-S0ph[i-1,j]+S0h[i-1,j]*W-S0ph[i-1,j]*RiskRed*Lh-S0ph[i-1,j]*(Woff+m[j,4])
      I1ph[i,j]<-I1ph[i-1,j]+S0ph[i-1,j]*RiskRed*Lh-I1ph[i-1,j]*(m[j,5]+mVIH35+PTrp)-I1ph[i-1,j]*I12
      I2ph[i,j]<-I2ph[i-1,j]+I1ph[i-1,j]*I12-I2ph[i-1,j]*I23-I2ph[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I3ph[i,j]<-I3ph[i-1,j]+I2ph[i-1,j]*I23-I3ph[i-1,j]*I34-I3ph[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I4ph[i,j]<-I4ph[i-1,j]+I3ph[i-1,j]*I34-I4ph[i-1,j]*I45-I4ph[i-1,j]*(m[j,5]+mVIH35+PTrp)   
      I5ph[i,j]<-I5ph[i-1,j]+I4ph[i-1,j]*I45-I5ph[i-1,j]*I56-I5ph[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I6ph[i,j]<-I6ph[i-1,j]+I5ph[i-1,j]*I56-I6ph[i-1,j]*I67-I6ph[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I7ph[i,j]<-I7ph[i-1,j]+I6ph[i-1,j]*I67-I7ph[i-1,j]*I78-I7ph[i-1,j]*(m[j,5]+mVIH35+PTrp)
      I8ph[i,j]<-I8ph[i-1,j]+I7ph[i-1,j]*I78-I8ph[i-1,j]*I89-I8ph[i-1,j]*(m[j,5]+mVIH235+PTrp)
      AIDSph[i,j]<-AIDSph[i-1,j]+I8ph[i-1,j]*I89-AIDSph[i-1,j]*(m[j,5]+mAIDS+PTrp)  
      
      PTrp=0.95
      #On treatment
      Tr[i,j]<-Tr[i-1,j]+sum(I1vh[i-1,j], I2vh[i-1,j], I3vh[i-1,j],
                             I4vh[i-1,j], I5vh[i-1,j], I6vh[i-1,j],
                             I1h[i-1,j], I2h[i-1,j], I3h[i-1,j],
                             I4h[i-1,j], I5h[i-1,j], I6h[i-1,j],
                             I1l[i-1,j], I2l[i-1,j], I3l[i-1,j],
                             I4l[i-1,j], I5l[i-1,j], I6l[i-1,j],
                             I1vl[i-1,j], I2vl[i-1,j], I3vl[i-1,j],
                             I4vl[i-1,j], I5vl[i-1,j], I6vl[i-1,j])*PTr[3]
      +sum(I1pvh[i-1,j], I2pvh[i-1,j], I3pvh[i-1,j],
           I4pvh[i-1,j], I5pvh[i-1,j], I6pvh[i-1,j],
           I1ph[i-1,j], I2ph[i-1,j], I3ph[i-1,j],
           I4ph[i-1,j], I5ph[i-1,j], I6ph[i-1,j], I7pvh[i-1,j], I8pvh[i-1,j],
           I7ph[i-1,j], I8ph[i-1,j], AIDSpvh[i-1,j], AIDSph[i-1,j])*PTrp
      +sum(I7vh[i-1,j], I8vh[i-1,j],I7h[i-1,j], I8h[i-1,j], I7l[i-1,j],
           I8l[i-1,j],I7vl[i-1,j], I8vl[i-1,j])*PTr[2]
      +sum(AIDSvh[i-1,j], AIDSh[i-1,j], AIDSl[i-1,j],AIDSvl[i-1,j])*PTr[1]
      -Tr[i-1,j]*m[j,4]
      
      #Deaths
      Dead[i,j]<-Dead[i-1,j]+sum(I1vh[i-1,j], I2vh[i-1,j], I3vh[i-1,j],
                             I4vh[i-1,j], I5vh[i-1,j], I6vh[i-1,j],
                             I1h[i-1,j], I2h[i-1,j], I3h[i-1,j],
                             I4h[i-1,j], I5h[i-1,j], I6h[i-1,j],
                             I1l[i-1,j], I2l[i-1,j], I3l[i-1,j],
                             I4l[i-1,j], I5l[i-1,j], I6l[i-1,j],
                             I1vl[i-1,j], I2vl[i-1,j], I3vl[i-1,j],
                             I4vl[i-1,j], I5vl[i-1,j], I6vl[i-1,j],
                             I1pvh[i-1,j], I2pvh[i-1,j], I3pvh[i-1,j],
                             I4pvh[i-1,j], I5pvh[i-1,j], I6pvh[i-1,j],
                             I1ph[i-1,j], I2ph[i-1,j], I3ph[i-1,j],
                             I4ph[i-1,j], I5ph[i-1,j], I6ph[i-1,j])*(m[j,5]+mVIH35)
      +sum(I7vh[i-1,j], I8vh[i-1,j],I7h[i-1,j], I8h[i-1,j], I7l[i-1,j],
           I8l[i-1,j],I7vl[i-1,j], I8vl[i-1,j], I7pvh[i-1,j], I8pvh[i-1,j],
           I7ph[i-1,j], I8ph[i-1,j])*(m[j,5]+mVIH235)
      +sum(AIDSvh[i-1,j], AIDSh[i-1,j], AIDSl[i-1,j],AIDSvl[i-1,j],
           AIDSpvh[i-1,j], AIDSph[i-1,j])*(m[j,5]+mAIDS)
      +Tr[i-1,j]*m[j,4]
      
      #New AIDS cases per year
      CalAIDS<-matrix(data=0, nrow=NoSem, ncol=Age)
      CalAIDS[i,j]<-sum(I8vh[i-1,j]*I89, I8h[i-1,j]*I89, I8l[i-1,j]*I89, I8vl[i-1,j]*I89, I8pvh[i-1,j]*I89, 
                        I8ph[i-1,j]*I89)
      
      #Force of infection
      Lvh_mat[i,j]<-Lvh
      Lh_mat[i,j]<-Lh
      Ll_mat[i,j]<-Ll
      Lvl_mat[i,j]<-Lvl
      
      
      #New HIV infections
      NewHIV<-matrix(data=0, nrow=NoSem, ncol=Age)
      NewHIV[i,j]<-sum(S0vh[i-1,j]*Lvh, S0h[i-1,j]*Lh, S0l[i-1,j]*Ll, S0vl[i-1,j]*Lvl, S0pvh[i-1,j]*Lvh*RiskRed, 
                       S0pvh[i-1,j]*Lh*RiskRed)
      
      #New HIV treated
    }
  }

P<-as.data.frame(cbind(1:NoSem,S0vh, S0h, S0l,S0vl,S0pvh,S0ph,I1vh, I1h, I1l, I1vl, I1pvh,I1ph,I2vh, I2h, 
                         I2l, I2vl, I2pvh,I2ph,I3vh, I3h, I3l, I3vl, I3pvh,I3ph, I4vh, I4h, I4l, 
                         I4vl, I4pvh, I4ph,I5vh, I5h, I5l, I5vl, I5pvh, I5ph,I6vh, I6h, I6l, I6vl, 
                         I6pvh,I6ph,I7vh, I7h, I7l, I7vl, I7pvh, I7ph, I8vh, I8h, I8l, I8vl, I8pvh, 
                         I8ph,AIDSvh, AIDSh, AIDSl, AIDSvl, AIDSpvh,AIDSph,Tr,Dead, CalAIDS, Lvh_mat, Lh_mat,
                       Ll_mat, Lvl_mat, NewHIV))

Q[[t]]<-P
}

Q

}