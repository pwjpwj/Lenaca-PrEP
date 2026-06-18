SettingArray<-function(NoSem=52, Age=85){
# Compartments not on PrEP s
# Very high risk
S0vh<-matrix(data=0, nrow=NoSem, ncol=Age)
S0vh[1,1]<-nvh*NumHSH
I1vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I1vh[1,1]<-NumHSHVIH*nvh*Rest/6
I2vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I2vh[1,1]<-NumHSHVIH*nvh*Rest/6
I3vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I3vh[1,1]<-NumHSHVIH*nvh*Rest/6
I4vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I4vh[1,1]<-NumHSHVIH*nvh*Rest/6
I5vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I5vh[1,1]<-NumHSHVIH*nvh*Rest/6
I6vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I6vh[1,1]<-NumHSHVIH*nvh*Rest/6
I7vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I7vh[1,1]<-NumHSHVIH*I7CD4*nvh
I8vh<-matrix(data=0,nrow=NoSem, ncol=Age)
I8vh[1,1]<-NumHSHVIH*I8CD4*nvh
AIDSvh<-matrix(data=0,nrow=NoSem, ncol=Age)
AIDSvh[1,1]<-NumHSHVIH*ACD4*nvh

#high risk
S0h<-matrix(data=0,nrow=NoSem, ncol=Age) 
S0h[1,1]<-nh*NumHSH
I1h<-matrix(data=0,nrow=NoSem, ncol=Age)
I1h[1,1]<-NumHSHVIH*nh*Rest/6
I2h<-matrix(data=0,nrow=NoSem, ncol=Age)
I2h[1,1]<-NumHSHVIH*nh*Rest/6
I3h<-matrix(data=0,nrow=NoSem, ncol=Age)
I3h[1,1]<-NumHSHVIH*nh*Rest/6
I4h<-matrix(data=0,nrow=NoSem, ncol=Age)
I4h[1,1]<-NumHSHVIH*nh*Rest/6
I5h<-matrix(data=0,nrow=NoSem, ncol=Age)
I5h[1,1]<-NumHSHVIH*nh*Rest/6
I6h<-matrix(data=0,nrow=NoSem, ncol=Age)
I6h[1,1]<-NumHSHVIH*nh*Rest/6
I7h<-matrix(data=0,nrow=NoSem, ncol=Age)
I7h[1,1]<-NumHSHVIH*nh*I7CD4
I8h<-matrix(data=0,nrow=NoSem, ncol=Age)
I6h[1,1]<-NumHSHVIH*nh*I8CD4
AIDSh<-matrix(data=0,nrow=NoSem, ncol=Age)
AIDSh[1,1]<-NumHSHVIH*nh*ACD4
#low risk
S0l<-matrix(data=0,nrow=NoSem, ncol=Age)
S0l[1,1]<-nl*NumHSH
I1l<-matrix(data=0,nrow=NoSem, ncol=Age)
I1l[1,1]<-NumHSHVIH*nl*Rest/6
I2l<-matrix(data=0,nrow=NoSem, ncol=Age)
I2l[1,1]<-NumHSHVIH*nl*Rest/6
I3l<-matrix(data=0,nrow=NoSem, ncol=Age)
I3l[1,1]<-NumHSHVIH*nl*Rest/6
I4l<-matrix(data=0,nrow=NoSem, ncol=Age)
I4l[1,1]<-NumHSHVIH*nl*Rest/6
I5l<-matrix(data=0,nrow=NoSem, ncol=Age)
I5l[1,1]<-NumHSHVIH*nl*Rest/6
I6l<-matrix(data=0,nrow=NoSem, ncol=Age)
I6l[1,1]<-NumHSHVIH*nl*Rest/6
I7l<-matrix(data=0,nrow=NoSem, ncol=Age)
I7l[1,1]<-NumHSHVIH*nl*I7CD4
I8l<-matrix(data=0,nrow=NoSem, ncol=Age)
I8l[1,1]<-NumHSHVIH*nl*I8CD4
AIDSl<-matrix(data=0,nrow=NoSem, ncol=Age)
AIDSl[1,1]<-NumHSHVIH*nl*ACD4
#very low risk
S0vl<-matrix(data=0, nrow=NoSem, ncol=Age) 
S0vl[1,1]<-nvl*NumHSH
I1vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I1vl[1,1]<-NumHSHVIH*nvl*Rest/6
I2vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I2vl[1,1]<-NumHSHVIH*nvl*Rest/6
I3vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I3vl[1,1]<-NumHSHVIH*nvl*Rest/6
I4vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I4vl[1,1]<-NumHSHVIH*nvl*Rest/6
I5vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I5vl[1,1]<-NumHSHVIH*nvl*Rest/6
I6vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I6vl[1,1]<-NumHSHVIH*nvl*Rest/6
I7vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I7vl[1,1]<-NumHSHVIH*nvl*I7CD4
I8vl<-matrix(data=0, nrow=NoSem, ncol=Age)
I8vl[1,1]<-NumHSHVIH*nvl*I8CD4
AIDSvl<-matrix(data=0, nrow=NoSem, ncol=Age)
AIDSvl[1,1]<-NumHSHVIH*nvl*ACD4

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

#Treated
Tr<-matrix(data=0, nrow=NoSem, ncol=Age)

grouplist<-list(S0vh, I1vh, I2vh, I3vh, I4vh, I5vh, I6vh, I7vh, I8vh, AIDSvh,
                S0h, I1h, I2h, I3h, I4h, I5h, I6h, I7h, I8h, AIDSh,
                S0l, I1l, I2l, I3l, I4l, I5l, I6l, I7l, I8l, AIDSl,
                S0vl, I1vl, I2vl, I3vl, I4vl, I5vl, I6vl, I7vl, I8vl, AIDSvl,
                S0pvh, I1pvh, I2pvh, I3pvh, I4pvh, I5pvh, I6pvh, I7pvh, I8pvh, AIDSpvh,
                S0ph, I1ph, I2ph, I3ph, I4ph, I5ph, I6ph, I7ph, I8ph, AIDSph, Tr)

groupnames<-list("S0vh", "I1vh", "I2vh","I3vh","I4vh","I5vh","I6vh","I7vh",
                 "I8vh","AIDSvh","S0h","I1h","I2h","I3h","I4h","I5h","I6h",
                 "I7h","I8h","AIDSh","S0l","I1l","I2l","I3l","I4l","I5l",
                 "I6l","I7l","I8l","AIDSl","S0vl","I1vl","I2vl","I3vl","I4vl",
                 "I5vl","I6vl","I7vl","I8vl","AIDSvl","S0pvh","I1pvh","I2pvh",
                 "I3pvh","I4pvh","I5pvh","I6pvh","I7pvh","I8pvh","AIDSpvh",
                 "S0ph","I1ph","I2ph","I3ph","I4ph","I5ph","I6ph","I7ph",
                 "I8ph","AIDSph","Tr")

D<-array(grouplist, dim=c(52,85,61), dimnames=list(NULL, NULL, groupnames))

D
}
