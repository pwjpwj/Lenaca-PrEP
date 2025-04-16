
HIVSIMPSA<-function(repeats){
# Dinamic HIV transmission model in HSH in Spain
##Values

for (k in 1:repeats){
  NumHSHVIH<-invgamma(Bestparam[1],1000)
  NumHSH<-892955-NumHSHVIH#(INE)
  #Encuesta hospitalaria
  
  #Distribucion de los pacientes por rango de CD4 
  ACD4<-invbeta(Bestparam[2]*0.1, Bestparam[2]*0.1)   #<200
  I8CD4<-invbeta(Bestparam[3]*0.1,Bestparam[3]*0.1) #200-350
  I7CD4<-invbeta(Bestparam[4]*0.1,Bestparam[4]*0.1)#350-500
  Rest<-invbeta(Bestparam[5]*0.1, Bestparam[5]*0.1)   #>500
  PropTR<-invbeta(Bestparam[6],0.1)
  G<-invbeta(Bestparam[7]/52,Bestparam[7]/52) #Population growth rate 1/1000 anually
  Woff<-0/52 #Proportion of people who go off PrEP 30% anually
  W<-0/52    #Proportion of people who go on PrEP 30% anually
  m<-read.csv("MortPobGen.csv", header=TRUE, sep=";")  #Mortality general population
  mVIH235<-invbeta(Bestparam[8], 0.000061) #Alejos et al Medicine 2016
  mVIH35<-invbeta(Bestparam[9], 0.000048)#Alejos et al Medicine 2016
  mAIDS<-invbeta(Bestparam[10]/52, Bestparam[10]/52)    #Parastu et al PLOSone 2018
  Cvh<-invgamma(Bestparam[11]/52, Bestparam[11]/52) #Partner change rate in the very high group 18-123/year Nichols et al. 2016
  Ch<-invbeta(Bestparam[12]/52, Bestparam[12]/52)     #Partner change rate in the high group 5-18/year Nichols et al. 2016
  Cl<-invbeta(Bestparam[13]/52, Bestparam[13]/52)    #Partner change rate in the low group 0.5-5/year Nichols et al. 2016
  Cvl<-invbeta(Bestparam[14]/52, Bestparam[14]/52)  #Partner change rate in the very low group 0.04-0.5/year Nichols et al. 2016
  
  nvh<-invbeta(Bestparam[15], 0.1)     #Proportion of the total population in each group 5-15% Nichols et al. 2016
  nh<-invbeta(Bestparam[16], 0.01)     # previo 0.4---30-50% Nichols et al. 2016
  nl<-invbeta(Bestparam[17], 0.1)     # previo 0.3---10-35% Nichols et al. 2016
  nvl<-invbeta(Bestparam[18],0.01)    # previo 0.1---4-46%
  
  B<-invbeta(Bestparam[19], 0.024)  #per partnership transmissibility Nichols et al 2016.
  I12<-invbeta(1-exp(-Bestparam[20]),1-exp(-Bestparam[20])) #1 month duration in stage 1 of HIV infection
  I23<-invbeta(1-exp(-Bestparam[21]),1-exp(-Bestparam[21])) #1 month duration in stage 2 of HIV infection
  I34<-invbeta(1-exp(-Bestparam[22]),1-exp(-Bestparam[22])) #1 month duration in stage 3 of HIV infection
  I45<-invbeta(1-exp(-Bestparam[23]/3),1-exp(-Bestparam[23]/3)) #3 month duration in stage 4 of HIV infection
  I55<-invbeta(1-exp(-Bestparam[25]/3),1-exp(-Bestparam[25]/3))  #3 month duration in stage 5 of HIV infection
  I56<-invbeta(1-exp(-Bestparam[26]/3),1-exp(-Bestparam[26]/3))  #3 month duration in stage 6 of HIV infection
  I67<-invbeta(1-exp(-Bestparam[27]/3),1-exp(-Bestparam[27]/3))  #Chronic duration in stage 7 of HIV infection
  I78<-invbeta(1-exp(-Bestparam[28]/(52*3)),1-exp(-Bestparam[28]/(52*3))) #Chronic duration in stage 9 of HIV infection
  I89<-invbeta(1-exp(-Bestparam[29]/(52*3.5)),1-exp(-Bestparam[29]/(52*3.5)))  #Chronic duration in stage 10 of HIV infection
  PTr<-c(invbeta(Bestparam[30],Bestparam[30]),invbeta(Bestparam[31],Bestparam[31]),
         invbeta(Bestparam[32],Bestparam[32])) #Nuñez et al AIDS 2018     
  ##Simulating
}

Q<-HIVSIM(Age, NoSem, NumHSH, NumHSHVIH, 
          ACD4, I8CD4, I7CD4, Rest, 
          G, Woff, W, 
          m, mAIDS, mVIH235, mVIH35, mTr, 
          Cvh, Ch, Cl, Cvl, 
          nvh, nh, nl, nvl, B,
          I12, I23, I34, 
          I45, I55, 
          I56, I67, 
          I78, I89, PTr, years, PropTR, DistEdad)

Q

}
