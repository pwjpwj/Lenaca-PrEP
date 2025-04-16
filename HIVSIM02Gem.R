HIVSIM <- function(Age = get("Age", globalenv()), 
                   NoSem = get("NoSem", globalenv()),
                   NumHSH = get("NumHSH", globalenv()), 
                   NumHSHVIH = get("NumHSHVIH", globalenv()),
                   ACD4 = get("ACD4", globalenv()), 
                   Rest = get("Rest", globalenv()),  # Removed I8CD4, I7CD4 as they are not used
                   G = get("G", globalenv()), 
                   Woff = get("Woff", globalenv()),
                   W = get("W", globalenv()), 
                   m = get("m", globalenv()),
                   mVIH235 = get("mVIH235", globalenv()), 
                   mVIH35 = get("mVIH35", globalenv()),
                   mAIDS = get("mAIDS", globalenv()),
                   mTr = get("mTr", globalenv()), 
                   Cvh = get("Cvh", globalenv()),
                   Ch = get("Ch", globalenv()), 
                   Cl = get("Cl", globalenv()),
                   Cvl = get("Cvl", globalenv()),
                   nvh = get("nvh", globalenv()), 
                   nh = get("nh", globalenv()),
                   nl = get("nl", globalenv()), 
                   nvl = get("nvl", globalenv()),
                   B = get("B", globalenv()), 
                   I12 = get("I12", globalenv()),
                   I23 = get("I23", globalenv()), 
                   I34 = get("I34", globalenv()),
                   I45 = get("I45", globalenv()), 
                   I56 = get("I56", globalenv()),
                   I67 = get("I67", globalenv()),
                   I78 = get("I78", globalenv()), 
                   I89 = get("I89", globalenv()),
                   PTr = get("PTr", globalenv()), 
                   years = get("years", globalenv()),
                   PropTR = getAnywhere(PropTR), 
                   DistEdad = getAnywhere(DistEdad),
                   RiskRed = getAnywhere(RiskRed)) {
  
  # ... (Q, compartments, Tr, Dead, L_mat, CalAIDS, NewHIV, NewTreatments initialization - same as before)
  
  # ... (Initialize compartments and starting state - same as before)
  
  for (t in 1:years) {
    for (j in 1:Age) {
      # ... (Age progression - same as before)
      
      for (i in 2:NoSem) {
        # ... (Prevalence, G, L, P calculations - same as before)
        
        # Differential equations (CORRECTED I TRANSITION RATES)
        for (risk in c("vh", "h", "l", "vl")) {
          for (prep in c("", "p")) {
            for (stage in c("S0", paste0("I", 1:8), "AIDS")) {
              compartment_index <- which(dimnames(compartments)[[5]] == stage)
              prev_stage_index <- if (stage == "S0") NA else which(dimnames(compartments)[[5]] == (if (stage == "I1") "S0" else paste0("I", as.numeric(gsub("I", "", stage)) - 1)))
              next_stage_index <- if (stage == "AIDS") NA else which(dimnames(compartments)[[5]] == paste0("I", as.numeric(gsub("I", "", stage)) + 1))
              
              I_transition_rate <- function(n) { # Corrected I transition rate function
                switch(as.character(n),  # Use as.character for switch
                       "1" = I12,
                       "2" = I23,
                       "3" = I34,
                       "4" = I45,
                       "5" = I56,
                       "6" = I67,
                       "7" = I78,
                       "8" = I89,
                       0 # Default case (important!)
                )
              }
              
              compartments[i, j, risk, prep, stage] <- compartments[i - 1, j, risk, prep, stage] +
                if (stage == "S0" && prep == "") compartments[i - 1, j, risk, prep, stage] * G[risk] + if (prep == "p") compartments[i - 1, j, risk, "", stage] * W else 0 +
                if (stage == "S0" && prep == "p") compartments[i - 1, j, risk, "", stage] * W - compartments[i - 1, j, risk, prep, stage] * (if (prep == "") L[risk] + W + m[j, 4] else Woff + m[j, 4]) +
                if (!is.na(prev_stage_index)) compartments[i - 1, j, risk, prep, prev_stage_index] * (if (stage != "S0") I_transition_rate(compartment_index - 1) else 1) -
                if (!is.na(next_stage_index)) compartments[i - 1, j, risk, prep, stage] * I_transition_rate(compartment_index) -
                compartments[i - 1, j, risk, prep, stage] * (m[j, 5] + (if (stage %in% paste0("I", 1:6)) mVIH35 else if (stage %in% c("I7", "I8")) mVIH235 else mAIDS) + PTr[if (stage == "AIDS") 1 else if (stage %in% c("I7", "I8")) 2 else 3])
              
            }
          }
        }
        
        # Corrected Tr update
        Tr[i, j] <- Tr[i - 1, j] + sum(compartments[i - 1, j, , , "I1"]) * PTr[3] +  # Only I1 contributes to Tr
          sum(compartments[i - 1, j, , , "I2"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I3"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I4"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I5"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I6"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I7"]) * PTr[3] +
          sum(compartments[i - 1, j, , , "I8"]) * PTr[3]
        
        L_mat[i, j, ] <- L
        NewHIV[i, j] <- sum(compartments[i - 1, j, , "", "S0"] * L * c(1, 1, 1, 1), compartments[i - 1, j, , "p", "S0"] * L * RiskRed * c(1, 1, 1, 1))
        
        # ... (Dead, CalAIDS, NewTreatments calculations - similar array operations)
        Dead[i, j] <- Dead[i - 1, j] + sum(compartments[i - 1, j, , , "AIDS"])*(mAIDS + mTr) + sum(compartments[i - 1, j, , , c(paste0("I", 1:8))])*(mVIH35 + mVIH235) + sum(compartments[i - 1, j, , , "S0"])*m[j, 4]
        CalAIDS[i, j] <- CalAIDS[i - 1, j] + sum(compartments[i - 1, j, , , "AIDS"])
        NewTreatments[i, j] <- NewTreatments[i - 1, j] + sum(compartments[i - 1, j, , , c(paste0("I", 1:8))]) * PTr[3]
      }
    }
    # Corrected and complete Q[[t]] assignment
    Q[[t]] <- as.data.frame(aperm(compartments[,,,,,1:9], c(1,3,4,5,2)))
    Q[[t]]$Tr <- Tr
    Q[[t]]$Dead <- Dead
    Q[[t]]$CalAIDS <- CalAIDS
    Q[[t]]$Lvh_mat <- L_mat[,, "vh"]
    Q[[t]]$Lh_mat <- L_mat[,, "h"]
    Q[[t]]$Ll_mat <- L_mat[,, "l"]
    Q[[t]]$Lvl_mat <- L_mat[,, "vl"]
    Q[[t]]$NewHIV <- NewHIV
    Q[[t]]$NewTreatments <- NewTreatments
    Q[[t]]$Time <- 1:NoSem
    
  }
  Q
}