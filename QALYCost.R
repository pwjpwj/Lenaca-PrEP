QALYCost<-function(Q){
  
  QQ<-array(unlist(Q), dim=c(52,5781,years))
  Q1<-NULL
  for (t in 1:years){
    Q1<-rbind(Q1, as.data.frame(cbind(QQ[1:52,1,t],rowSums(QQ[1:52,2:511,t]),
                                      rowSums(QQ[1:52,512:1021,t]),
                                      rowSums(QQ[1:52,1022:1531,t]),
                                      rowSums(QQ[,1532:2041,t]),
                                      rowSums(QQ[,2042:2551,t]),
                                      rowSums(QQ[,2552:3061,t]),
                                      rowSums(QQ[,3062:3571,t]),
                                      rowSums(QQ[,3572:4081,t]),
                                      rowSums(QQ[,4082:4591,t]), 
                                      rowSums(QQ[,4592:5101,t]),
                                      rowSums(QQ[,5102:5186,t]),
                                      rowSums(QQ[,5187:5271,t]),
                                      rowSums(QQ[,5272:5356,t]),
                                      rowSums(QQ[,5357:5441,t]),
                                      rowSums(QQ[,5442:5526,t]),
                                      rowSums(QQ[,5527:5611,t]),
                                      rowSums(QQ[,5612:5696,t]),
                                      rowSums(QQ[,5697:5781,t]))))
  }
  colnames(Q1)[c(1:13,18,19)]<-c("Weeks", "Susceptible", "Inf1", "Inf2", "Inf3", "Inf4", "Inf5", "Inf6", "Inf7", "Inf8", "AIDS", "Treated", "Dead", "CalAIDS","NewHIV")
  Q1$Weeks<-1:(years*52)
  QALYs<-0
  Costs<-0
  x<-Q1[,c("Inf2", "Inf3", "Inf4", "Inf5","Inf6", "Inf7")]
  for (i in 1:length(x)){QALYs<-QALYs+QalyCalc(x[,i],w=WM350)}
  QALYs<-QALYs+QalyCalc(Q1[,"Inf1"],w=WHIV)
  QALYs<-QALYs+QalyCalc(Q1[,"Susceptible"],w=WSusceptible)
  QALYs<-QALYs+QalyCalc(Q1[,"Inf8"],w=WHIV)
  QALYs<-QALYs+QalyCalc(Q1[,"AIDS"],w=WAIDS)
  QALYs<-QALYs+QalyCalc(Q1[,"Treated"],w=WTREAT)
  
  Costs<-Costs+CostCalc(Q1[,"AIDS"],C=AIDScost)
  Costs<-Costs+CostCalc(Q1[,"Treated"],C=Treatcost)
  
cbind(QALYs, Costs)

}