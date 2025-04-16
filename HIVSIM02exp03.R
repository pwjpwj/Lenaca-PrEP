HIVSIM<-function(Age=Age , NoSem=NoSem,
                 NumHSH=NumHSH, NumHSHVIH=NumHSHVIH, ACD4=ACD4, I8CD4=I8CD4,
                 I7CD4=I7CD4, Rest=Rest, G=G, Woff=Woff, W=W, m=m,
                 mVIH235=mVIH235, mVIH35=mVIH35, mAIDS=mAIDS, mTr=mTr,
                 Cvh=Cvh, Ch=Ch, Cl=Cl, Cvl=Cvl, nvh=nvh, nh=nh, 
                 nl=nl, nvl=nvl, B=B, I12=I12, I23=I23, I34=I34, 
                 I45=I45, I55=I55, I56=I56, I67=I67, I78=I78, I89=I89,
                 PTr=PTr, years=years, PropTR=PropTR, DistEdad=DistEdad,
                 RiskRed=RiskRed){
  
  #The list to store the yearly arrays
  Q<-vector("list", years)
  
  Compartments<-array(data=0, c(NoSem, Age,10,4,2))
  #first dimension weeks, second dimension age, third dimension S0---AIDS, forth dimension risk group and fifth prep

    #Treated
    Tr<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #Dead
    Dead<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #Force of infection
    Lvh_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Lh_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Ll_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    Lvl_mat<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #New AIDS cases per year
    CalAIDS<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #New HIV infections
    NewHIV<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #New HIV Treatments
    NewTreatments<-matrix(data=0, nrow=NoSem, ncol=Age)
    
    #Probability that a partner will be in a specific group
    Gvh<-(Cvh*nvh)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gh<-(Ch*nh)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gl<-(Cl*nl)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    Gvl<-(Cvl*nvl)/(Cvh*nvh+Ch*nh+Cl+nl*Cvl+nvl)
    
    for (j in 1:Age){
      for (k in 1:4){
        RG<-c(nvh, nh, nl, nvl)
    #Defining starting state for the non prep compartments
    
    Compartments[1,j,1,k,1]<-RG[k]*NumHSH*DistEdad[j]#  S0vh[1,j]<-nvh*NumHSH*DistEdad[j]
    Compartments[1,j,2:7,k,1]<-NumHSHVIH*RG[k]*Rest/6*DistEdad[j]#  I1vh[1,j]<-NumHSHVIH*nvh*Rest/6*DistEdad[j]
    Compartments[1,j,8:9,k,1]<-NumHSHVIH*I7CD4*RG[k]*DistEdad[j] #I7vh[1,j]<-NumHSHVIH*I7CD4*nvh*DistEdad[j]
    Compartments[1,j,10,k,1]<-NumHSHVIH*RG[k]*ACD4*DistEdad[j] #AIDSvh[1,j]<-NumHSHVIH*nvh*ACD4*DistEdad[j]
      }
   
    
    Tr[1,j]<-NumHSHVIH*PropTR*DistEdad[j]
    }
    #Compartments on PrEP Already initiated as 0 when the array was initiated. Compartments of prep
    #are only 2, high and very high risk of transmission.

    
    
for (t in 1:years){
  for (j in 1:Age){
    if (t>=2 & j<85){
      for(k in 1:4){
        for(prep in 1:2){
      Compartments[1,j+1,1,k,prep] <-Compartments[52,j+1,1,k,prep]
      Compartments[1,j+1,2:10,k,prep] <-Compartments[52,j+1,2,k,prep]
        }
      }
     
      #Treated
      Tr[1,j+1]<-Tr[52,j]  
    }
    
  for (i in 2:NoSem){
      #Calculating the time changing variables
      
      #Prevalence of infection in very high risk group
      PrevInfvh<-sum(Compartments[i-1,,2:10,1,])/sum(Compartments[i-1,,,1,])
      #Prevalence of infection in high risk group
      PrevInfh<-sum(Compartments[i-1,,2:10,2,])/sum(Compartments[i-1,,,2,])
      #Prevalence of infection in low risk group
      PrevInfl<-sum(Compartments[i-1,,2:10,3,])/sum(Compartments[i-1,,,3,])
      #Prevalence of infection in very low risk group
      PrevInfvl<-sum(Compartments[i-1,,2:10,4,])/sum(Compartments[i-1,,,4,])
      
      #Probability that the partner is infectious
      P<-Gvh*PrevInfvh+Gh*PrevInfh+Gl*PrevInfl+Gvl*PrevInfvl
      
      #Force of infection
      Lvh<-Cvh*B*P
      Lh<-Ch*B*P
      Ll<-Cl*B*P
      Lvl<-Cvl*B*P
      
      
      FOI<-c(Lvh, Lh, Ll, Lvl)
      
      #Differential Equations each time step is one week
      #Not on PrEP, very high risk group
      for (k in 1:4){
        
      Compartments[i,j,1,k,1]<-Compartments[i-1,j,1,k,1]+Compartments[i-1,j,1,k,1]*G+Compartments[i-1,j,1,k,2]*Woff-Compartments[i-1,j,1,k,1]*(FOI[k]+W+m[j,4])
      Compartments[i,j,2,k,1]<-Compartments[i-1,j,2,k,1]+Compartments[i-1,j,1,k,1]*FOI[k]-Compartments[i-1,j,2,k,1]*(m[j,5]+mVIH35+PTr[3])-Compartments[i-1,j,2,k,1]*I12
      Compartments[i,j,3,k,1]<-Compartments[i-1,j,3,k,1]+Compartments[i-1,j,2,k,1]*I12-Compartments[i-1,j,3,k,1]*I23-Compartments[i-1,j,3,k,1]*(m[j,5]+mVIH35+PTr[3])
      Compartments[i,j,4,k,1]<-Compartments[i-1,j,4,k,1]+Compartments[i-1,j,3,k,1]*I23-Compartments[i-1,j,4,k,1]*I34-Compartments[i-1,j,4,k,1]*(m[j,5]+mVIH35+PTr[3])
      Compartments[i,j,5,k,1]<-Compartments[i-1,j,5,k,1]+Compartments[i-1,j,4,k,1]*I34-Compartments[i-1,j,5,k,1]*I45-Compartments[i-1,j,5,k,1]*(m[j,5]+mVIH35+PTr[3])   
      Compartments[i,j,6,k,1]<-Compartments[i-1,j,6,k,1]+Compartments[i-1,j,5,k,1]*I45-Compartments[i-1,j,6,k,1]*I56-Compartments[i-1,j,6,k,1]*(m[j,5]+mVIH35+PTr[3])
      Compartments[i,j,7,k,1]<-Compartments[i-1,j,7,k,1]+Compartments[i-1,j,6,k,1]*I56-Compartments[i-1,j,7,k,1]*I67-Compartments[i-1,j,7,k,1]*(m[j,5]+mVIH35+PTr[3])
      Compartments[i,j,8,k,1]<-Compartments[i-1,j,8,k,1]+Compartments[i-1,j,7,k,1]*I67-Compartments[i-1,j,8,k,1]*I78-Compartments[i-1,j,8,k,1]*(m[j,5]+mVIH35+PTr[2])
      Compartments[i,j,9,k,1]<-Compartments[i-1,j,9,k,1]+Compartments[i-1,j,8,k,1]*I78-Compartments[i-1,j,9,k,1]*I89-Compartments[i-1,j,9,k,1]*(m[j,5]+mVIH235+PTr[2])
      Compartments[i,j,10,k,1]<-Compartments[i-1,j,10,k,1]+Compartments[i-1,j,9,k,1]*I89-Compartments[i-1,j,10,k,1]*(m[j,5]+mAIDS+PTr[1])
      }
      
      # On PrEP
      for (k in 1:2){
          Compartments[i,j,1,k,2]<-Compartments[i-1,j,1,k,2]+Compartments[i-1,j,1,k,1]*W+Compartments[i-1,j,1,k,1]*Woff-Compartments[i-1,j,1,k,2]*(FOI[k]*RiskRed+W+m[j,4])
          Compartments[i,j,2,k,2]<-Compartments[i-1,j,2,k,2]+Compartments[i-1,j,1,k,2]*FOI[k]*RiskRed-Compartments[i-1,j,2,k,2]*(m[j,5]+mVIH35+PTr[3])-Compartments[i-1,j,2,k,2]*I12
          Compartments[i,j,3,k,2]<-Compartments[i-1,j,3,k,2]+Compartments[i-1,j,2,k,2]*I12-Compartments[i-1,j,3,k,2]*I23-Compartments[i-1,j,3,k,2]*(m[j,5]+mVIH35+PTr[3])
          Compartments[i,j,4,k,2]<-Compartments[i-1,j,4,k,2]+Compartments[i-1,j,3,k,2]*I23-Compartments[i-1,j,4,k,2]*I34-Compartments[i-1,j,4,k,2]*(m[j,5]+mVIH35+PTr[3])
          Compartments[i,j,5,k,2]<-Compartments[i-1,j,5,k,2]+Compartments[i-1,j,4,k,2]*I34-Compartments[i-1,j,5,k,2]*I45-Compartments[i-1,j,5,k,2]*(m[j,5]+mVIH35+PTr[3])   
          Compartments[i,j,6,k,2]<-Compartments[i-1,j,6,k,2]+Compartments[i-1,j,5,k,2]*I45-Compartments[i-1,j,6,k,2]*I56-Compartments[i-1,j,6,k,2]*(m[j,5]+mVIH35+PTr[3])
          Compartments[i,j,7,k,2]<-Compartments[i-1,j,7,k,2]+Compartments[i-1,j,6,k,2]*I56-Compartments[i-1,j,7,k,2]*I67-Compartments[i-1,j,7,k,2]*(m[j,5]+mVIH35+PTr[3])
          Compartments[i,j,8,k,2]<-Compartments[i-1,j,8,k,2]+Compartments[i-1,j,7,k,2]*I67-Compartments[i-1,j,8,k,2]*I78-Compartments[i-1,j,8,k,2]*(m[j,5]+mVIH35+PTr[2])
          Compartments[i,j,9,k,2]<-Compartments[i-1,j,9,k,2]+Compartments[i-1,j,8,k,2]*I78-Compartments[i-1,j,9,k,2]*I89-Compartments[i-1,j,9,k,2]*(m[j,5]+mVIH235+PTr[2])
          Compartments[i,j,10,k,2]<-Compartments[i-1,j,10,k,2]+Compartments[i-1,j,9,k,2]*I89-Compartments[i-1,j,10,k,2]*(m[j,5]+mAIDS+PTr[1])
      }
      
      
      #On treatment
      Tr[i,j]<-Tr[i-1,j]+sum(Compartments[i-1,j,2:7,1:4,1], Compartments[i-1,j,1:7,1:2,2])*PTr[3]
             +sum(Compartments[i-1,j,8:9,1:4,1], Compartments[i-1,j,8:9,1:2,1])*PTr[2]
      +sum(Compartments[i-1,j,10,1:4,1], Compartments[i-1,j,10,1:2,1])*PTr[1]
      -Tr[i-1,j]*m[j,4]
      
      #Deaths
      Dead[i,j]<-Dead[i-1,j]+sum(Compartments[i-1,j,2:7,1:4,1], Compartments[i-1,j,1:7,1:2,2])*(m[j,5]+mVIH35)
      +sum(Compartments[i-1,j,8:9,1:4,1], Compartments[i-1,j,8:9,1:2,1])*(m[j,5]+mVIH235)
      +sum(Compartments[i-1,j,10,1:4,1], Compartments[i-1,j,10,1:2,1])*(m[j,5]+mAIDS)
      -Tr[i-1,j]*m[j,4]
     
      #New AIDS cases per year
      CalAIDS[i,j]<-sum(Compartments[i-1,j,9,1:4,1]*I89, Compartments[i-1,j,9,1:2,2]*I89)
      
      #Force of infection
      Lvh_mat[i,j]<-Lvh
      Lh_mat[i,j]<-Lh
      Ll_mat[i,j]<-Ll
      Lvl_mat[i,j]<-Lvl
        
        
      #New HIV infections
      for (k in 1:4)
      NewHIV[i,j]<-sum(Compartments[i-1,j,1,k,1]*FOI[k], Compartments[i-1,j,1,k,2]*FOI[k]*RiskRed)
        
      #New HIV treated
      NewTreatments[i,j]<-sum(Compartments[i-1,j,2:7,1:4,1], Compartments[i-1,j,1:7,1:2,2])*PTr[3]
      +sum(Compartments[i-1,j,8:9,1:4,1], Compartments[i-1,j,8:9,1:2,1])*PTr[2]
      +sum(Compartments[i-1,j,10,1:4,1], Compartments[i-1,j,10,1:2,1])*PTr[1]
    }
  }

P<-as.data.frame(cbind(1:NoSem,Compartments[,,1,1,1], Compartments[,,1,2,1], Compartments[,,1,3,1],Compartments[,,1,4,1],Compartments[,,1,1,2],Compartments[,,1,2,2],
                       Compartments[,,2,1,1], Compartments[,,2,2,1], Compartments[,,2,3,1],Compartments[,,2,4,1],Compartments[,,2,1,2],Compartments[,,2,2,2],
                       Compartments[,,3,1,1], Compartments[,,3,2,1], Compartments[,,3,3,1],Compartments[,,3,4,1],Compartments[,,3,1,2],Compartments[,,3,2,2],
                       Compartments[,,4,1,1], Compartments[,,4,2,1], Compartments[,,4,3,1],Compartments[,,4,4,1],Compartments[,,4,1,2],Compartments[,,4,2,2],
                       Compartments[,,5,1,1], Compartments[,,5,2,1], Compartments[,,5,3,1],Compartments[,,5,4,1],Compartments[,,5,1,2],Compartments[,,5,2,2],
                       Compartments[,,6,1,1], Compartments[,,6,2,1], Compartments[,,6,3,1],Compartments[,,6,4,1],Compartments[,,6,1,2],Compartments[,,6,2,2],
                       Compartments[,,7,1,1], Compartments[,,7,2,1], Compartments[,,7,3,1],Compartments[,,7,4,1],Compartments[,,7,1,2],Compartments[,,7,2,2],
                       Compartments[,,8,1,1], Compartments[,,8,2,1], Compartments[,,8,3,1],Compartments[,,8,4,1],Compartments[,,8,1,2],Compartments[,,8,2,2],
                       Compartments[,,9,1,1], Compartments[,,9,2,1], Compartments[,,9,3,1],Compartments[,,9,4,1],Compartments[,,9,1,2],Compartments[,,9,2,2],
                       Compartments[,,10,1,1], Compartments[,,10,2,1], Compartments[,,10,3,1],Compartments[,,10,4,1],Compartments[,,10,1,2],Compartments[,,10,2,2],
                       Tr,Dead, CalAIDS, Lvh_mat, Lh_mat
                       , Ll_mat, Lvl_mat, NewHIV, NewTreatments))

Q[[t]]<-P
}

Q


}