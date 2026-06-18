QALYCostP<-function(QP){
  
  QQP<-array(unlist(QP), dim=c(52,5271,years))

  
  Q1P<-NULL
  for (t in 1:years){
    Q1P<-rbind(Q1P, as.data.frame(cbind(QQP[1:52,1,t],rowSums(QQP[1:52,2:341,t]),rowSums(QQP[1:52,342:511,t]),rowSums(QQP[1:52,512:1021,t]), rowSums(QQP[,1022:1531,t]), rowSums(QQP[,1532:2041,t]), rowSums(QQP[,2042:2551,t]), rowSums(QQP[,2552:3061,t]), rowSums(QQP[,3062:3571,t]), rowSums(QQP[,3572:4081,t]), rowSums(QQP[,4082:4591,t]), rowSums(QQP[,4592:4676,t]), rowSums(QQP[,4677:4761,t]))))
  }
  colnames(Q1P)<-c("Weeks", "Susceptible","PrEP", "Inf1", "Inf2", "Inf3", "Inf4", "Inf5", "Inf6", "Inf7", "AIDS", "Treated", "Dead")
  Q1P$Weeks<-1:(years*52)
  
  QALYsP<-0
  
  x<-Q1P[,c("Inf2", "Inf3", "Inf4", "Inf5","Inf6", "Inf7")]
  for (i in 1:length(x)){QALYsP<-QALYsP+QalyCalc(x[,i],w=WM350)}
  QALYsP<-QALYsP+QalyCalc(Q1P[,"Inf1"],w=WHIV)
  QALYsP<-QALYsP+QalyCalc(Q1P[,"Susceptible"],w=WSusceptible)
  QALYsP<-QALYsP+QalyCalc(Q1P[,"PrEP"],w=WSusceptible)
  QALYsP<-QALYsP+QalyCalc(Q1P[,"AIDS"],w=WAIDS)
  QALYsP<-QALYsP+QalyCalc(Q1P[,"Treated"],w=WTREAT)
  
  
  QALYsP
  CostsP<-0
  CostsP<-CostsP+CostCalc(Q1P[,"AIDS"],C=AIDScost)
  CostsP<-CostsP+CostCalc(Q1P[,"Treated"],C=Treatcost)
  CostsP<-CostsP+CostCalc(Q1P[,"PrEP"],C=PrEPcost)
  
cbind(QALYsP, CostsP)

}
