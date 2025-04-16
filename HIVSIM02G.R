HIVSIM <- function(Age=get("Age", globalenv()) , NoSem=get("NoSem", globalenv()),
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
                   PropTR=getAnywhere(PropTR), DistEdad=getAnywhere(DistEdad),
                   RiskRed=getAnywhere(RiskRed)) {
  
  # Initialize state arrays: (NoSem, Age, Compartment index)
  # 12 compartments + Dead compartment
  compartments <- c("S0", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "AIDS", "Tr", "Dead")
  comp_matrix <- array(0, dim = c(NoSem, Age, length(compartments), 4))  
  # VH: 1, H: 2, L: 3, VL: 4  
  
  # Initial state setup with vectorized calculations
  dist_weight <- DistEdad * NumHSH
  comp_matrix[1, , 1, 1] <- nvh * dist_weight  
  comp_matrix[1, , 1, 2] <- nh * dist_weight  
  comp_matrix[1, , 1, 3] <- nl * dist_weight  
  comp_matrix[1, , 1, 4] <- nvl * dist_weight  
  
  # Vectorized setting of infected compartments for each risk group
  init_inf_factors <- NumHSHVIH * Rest / 6
  for (i in 2:8) {  # I1 to I8 compartments
    comp_matrix[1, , i, 1] <- init_inf_factors * nvh * DistEdad  
    comp_matrix[1, , i, 2] <- init_inf_factors * nh * DistEdad  
    comp_matrix[1, , i, 3] <- init_inf_factors * nl * DistEdad  
    comp_matrix[1, , i, 4] <- init_inf_factors * nvl * DistEdad  
  }
  
  # Loop optimization by processing all ages simultaneously
  for (t in 1:years) {  
    if (t > 1) {  
      for (j in 1:(Age - 1)) {  
        comp_matrix[1, j + 1, , ] <- comp_matrix[NoSem, j, , ]  
      }
    }
    
    # Force of infection calculations (example)
    FOI <- rowSums(comp_matrix[, , 2:9, ])  # Summing compartments contributing to infection
  }
  
  return(comp_matrix)  
}