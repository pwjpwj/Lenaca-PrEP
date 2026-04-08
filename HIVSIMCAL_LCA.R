HIVSIMCAL<-function(repeats){
# Dinamic HIV transmission model in HSH in Spain
##Values
#RiskRed<-0.1
DistEdad<-vector(length=85)
DistEdad[1:10]<-0.114106384/10
DistEdad[11:20]<-0.116645332/10
DistEdad[21:30]<-0.281029342/10
DistEdad[31:40]<-0.30552628/10
DistEdad[41:50]<-0.112127906/10
DistEdad[51:85]<-0.070564755/35
  
Age<-85
years<-5
NoSem<-52
R<-vector("list", repeats)
Params<-matrix(data=0, nrow=32, ncol=repeats)

for (k in 1:repeats){
NumHSHVIH<-rnorm(1,33558, 4000)
Params[1,k]<-NumHSHVIH

NumHSH<-892955-NumHSHVIH#(INE)
  #Encuesta hospitalaria

#Distribucion de los pacientes por rango de CD4 (Informe VIH SIDA 2023)
j<-rdirichlet(1,c(0.20,0.221,0.22, 0.359))
ACD4<-j[1,1]#*0.2 #<200
Params[2,k]<-ACD4
I8CD4<-j[1,2]#*0.1#200-350
Params[3,k]<-I8CD4
I7CD4<-j[1,3]#*0.1#350-500
Params[4,k]<-I7CD4
Rest<-j[1,4]#*0.1 #>500
Params[5,k]<-Rest
PropTR<-invbeta(0.89,0.1)#Proportion on treatment (cotinuo atención VIH 2021-2022)
Params[6,k]<-PropTR


G<-invbeta(0.001/52,0.001/52) #Population growth rate 1/1000 anually
Params[7,k]<-G 
Woff<-0/52 #Proportion of people who go off PrEP 30% anually
W<-0/52    #Proportion of people who go on PrEP 30% anually
#L<-0.05/52   #Force of infection
m<-read.csv("MortPobGen.csv", header=TRUE, sep=";")  #Mortality general population
mVIH235<-invbeta(0.00012692, 0.000061) #Alejos et al Medicine 2016
Params[8,k]<-mVIH235
mVIH35<-invbeta(0.0000538, 0.000048)   #Alejos et al Medicine 2016
Params[9,k]<-mVIH35
mAIDS<-invbeta(0.04/52, 0.04/52)    #Parastu et al PLOSone 2018
Params[10,k]<-mAIDS
Cvh<-invgamma(52.5/52, 52.5/52) #Partner change rate in the very high group 18-123/year Nichols et al. 2016
Params[11,k]<-Cvh
Ch<-invbeta(9/52, 9/52)     #Partner change rate in the high group 5-18/year Nichols et al. 2016
Params[12,k]<-Ch
Cl<-invbeta(2.5/52,2.5/52)    #Partner change rate in the low group 0.5-5/year Nichols et al. 2016
Params[13,k]<-Cl
Cvl<-invbeta(0.2/52, 0.2/52)  #Partner change rate in the very low group 0.04-0.5/year Nichols et al. 2016
Params[14,k]<-Cvl
j1<-rdirichlet(1,c(0.1,0.4,0.3,0.2))
nvh<-j1[1,1]     #Proportion of the total population in each group 5-15% Nichols et al. 2016
Params[15,k]<-nvh
nh<-j1[1,2]     #30-50% Nichols et al. 2016
Params[16,k]<-nh
nl<-j1[1,3]     #10-35% Nichols et al. 2016
Params[17,k]<-nl
nvl<-j1[1,4]    #4-46%
Params[18,k]<-nvl

B<-invbeta(0.024, 0.024)  #per partnership transmissibility Nichols et al 2016.
Params[19,k]<-B
I12<-invbeta(1-exp(-0.25),1-exp(-0.25)) #1 month duration in stage 1 of HIV infection
Params[20,k]<-I12
I23<-invbeta(1-exp(-0.25),1-exp(-0.25)) #1 month duration in stage 2 of HIV infection
Params[21,k]<-I23
I34<-invbeta(1-exp(-0.25),1-exp(-0.25)) #1 month duration in stage 3 of HIV infection
Params[22,k]<-I34
I45<-invbeta(1-exp(-0.25/3),1-exp(-0.25/3)) #3 month duration in stage 4 of HIV infection
Params[23,k]<-I45
I55<-invbeta(1-exp(-0.25/3),1-exp(-0.25/3))  #3 month duration in stage 5 of HIV infection
Params[25,k]<-I55
I56<-invbeta(1-exp(-0.25/3),1-exp(-0.25/3))  #3 month duration in stage 6 of HIV infection
Params[26,k]<-I56
I67<-invbeta(1-exp(-0.25/3),1-exp(-0.25/3))  #Chronic duration in stage 7 of HIV infection
Params[27,k]<-I67
I78<-invbeta(1-exp(-0.25/(52*3)),1-exp(-0.25/(52*3))) #Chronic duration in stage 9 of HIV infection
Params[28,k]<-I78
I89<-invbeta(1-exp(-0.25/(52*3.5)),1-exp(-0.25/(52*3.5)))  #Chronic duration in stage 10 of HIV infection
Params[29,k]<-I89

PTr<-c(invbeta(0.002644231,0.002644231),invbeta(0.001682692,0.001682692),
       invbeta(0.000528846,0.000528846)) #Nuñez et al AIDS 2018     
Params[30,k]<-PTr[1]
Params[31,k]<-PTr[2]
Params[32,k]<-PTr[3]

#RiskRed<-0.86

#200-500->500


##Simulating


Q<-HIVSIM(Age, NoSem, NumHSH, NumHSHVIH, 
          ACD4, I8CD4, I7CD4, Rest, 
          G, Woff, W, 
          m, mAIDS, mVIH235, mVIH35, 
          Cvh, Ch, Cl, Cvl, 
          nvh, nh, nl, nvl, B,
          I12, I23, I34, 
          I45, I55, 
          I56, I67, 
          I78, I89, PTr, years, PropTR, DistEdad,
          RiskRed)

head(Q[[3]][1:years])
QQ<-array(unlist(Q), dim=c(52,6291,20))


#Arranging the array into a megamatrix

Q1<-NULL
for (t in 1:years){
  Q1<-rbind(Q1, as.data.frame(cbind(QQ[1:52,1,t],rowSums(QQ[1:52,2:511,t]),rowSums(QQ[1:52,512:1021,t]),rowSums(QQ[1:52,1022:1531,t]), rowSums(QQ[,1532:2041,t]), rowSums(QQ[,2042:2551,t]), rowSums(QQ[,2552:3061,t]), rowSums(QQ[,3062:3571,t]), rowSums(QQ[,3572:4081,t]), rowSums(QQ[,4082:4591,t]), rowSums(QQ[,4592:5101,t]), rowSums(QQ[,5102:5611,t]),
                                    rowSums(QQ[,5612:5696,t]), rowSums(QQ[,5697:5780,t]), rowSums(QQ[,5781:5865,t]),
                                    rowSums(QQ[,5866:5950,t]), rowSums(QQ[,5951:6035,t]), rowSums(QQ[,6036:6120,t]), rowSums(QQ[,6121:6206,t]), rowSums(QQ[,6207:6291,t]))))
}


colnames(Q1)<-c("Weeks", "Susceptible", "Inf1", "Inf2", "Inf3", "Inf4", 
                "Inf5", "Inf6", "Inf7", "Inf8", "AIDS", "Treated", "Dead", "NewAIDS",
                "Lvh", "Lh", "Ll", "Lvl","NewHIV", "New_HIV_treatments")
Q1$Weeks<-1:(years*52)

R[[k]]<-Q1
}
Result<-vector("list",2)
Result[[1]]<-R 
Result[[2]]<-Params
Result
}

