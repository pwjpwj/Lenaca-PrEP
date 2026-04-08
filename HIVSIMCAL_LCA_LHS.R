HIVSIMCAL_LHS <- function(U_block, m, DistEdad, RiskRed) {
  # ---- constants 
  Age   <- 85L
  years <- 5L
  NoSem <- 52L
  
  # U_block: matrix of uniforms in (0,1), nrow = repeats, ncol = 27
  U_block <- as.matrix(U_block)
  repeats <- nrow(U_block)
  stopifnot(ncol(U_block) == 28L)
  
    # ---- outputs ----
  R      <- vector("list", repeats)
  Params <- matrix(0, nrow = 28L, ncol = repeats)
  
  Slices <- make_slices_HIVSIM(Age)
  stopifnot(Slices$ncol_expected == 5356L)
  
  # ---- main loop over simulations ----
  for (k in seq_len(repeats)) {
    # --- parameters from LHS row ---
    u_k   <- U_block[k, ]
    par_k <- draw_params_from_u(u_k)
    Params[, k] <- par_k$param_vec
    
    # --- simulate (same call as your original code) ---
    Q <- HIVSIM(
      Age, NoSem, par_k$NumHSH, par_k$NumHSHVIH,
      par_k$ACD4, par_k$I7CD4, par_k$Rest,
      par_k$G, par_k$Woff, par_k$W,
      m, par_k$mVIH235, par_k$mVIH35, par_k$mAIDS,
      par_k$Cvh, par_k$Ch, par_k$Cl, par_k$Cvl,
      par_k$nvh, par_k$nh, par_k$nl, par_k$nvl, par_k$B,
      par_k$I12, par_k$I23, par_k$I34,
      par_k$I45,
      par_k$I56, par_k$I67,
      par_k$I78, par_k$PTr, years, par_k$PropTR,
      DistEdad, RiskRed
    )
    
    QQ <- array(
      as.numeric(unlist(lapply(Q, as.matrix))),
      dim = c(52L, Slices$ncol_expected, years)
    )
    
    stopifnot(dim(QQ)[2] == Slices$ncol_expected)
    
    # --- build output matrix  ---
    nrows <- years * 52L
    Qmat  <- matrix(0, nrow = nrows, ncol = 19L)
    colnames(Qmat) <- c(
      "Weeks", "Susceptible", "Inf1", "Inf2", "Inf3", "Inf4",
      "Inf5", "Inf6", "Inf7", "AIDS", "Treated", "Dead", "NewAIDS",
      "Lvh", "Lh", "Ll", "Lvl", "New_HIV_inf", "New_HIV_treatments"
    )
    
    r0 <- 0L
    for (t in seq_len(years)) {
      rr <- r0 + seq_len(52L)
      Qmat[rr, 1L]  <- rr
      Qmat[rr, 2L]  <- rs_at(QQ, Slices$S0,   t)
      Qmat[rr, 3L]  <- rs_at(QQ, Slices$I1,   t)
      Qmat[rr, 4L]  <- rs_at(QQ, Slices$I2,   t)
      Qmat[rr, 5L]  <- rs_at(QQ, Slices$I3,   t)
      Qmat[rr, 6L]  <- rs_at(QQ, Slices$I4,   t)
      Qmat[rr, 7L]  <- rs_at(QQ, Slices$I5,   t)
      Qmat[rr, 8L]  <- rs_at(QQ, Slices$I6,   t)
      Qmat[rr, 9L]  <- rs_at(QQ, Slices$I7,   t)
      Qmat[rr,10L]  <- rs_at(QQ, Slices$AIDS, t)
      
      Qmat[rr,11L]  <- rs_at(QQ, Slices$Treated, t)
      Qmat[rr,12L]  <- rs_at(QQ, Slices$Dead,    t)
      Qmat[rr,13L]  <- rs_at(QQ, Slices$NewAIDS, t)
      
      Qmat[rr,14L]  <- rs_at(QQ, Slices$FI_vh, t)
      Qmat[rr,15L]  <- rs_at(QQ, Slices$FI_h,  t)
      Qmat[rr,16L]  <- rs_at(QQ, Slices$FI_l,  t)
      Qmat[rr,17L]  <- rs_at(QQ, Slices$FI_vl, t)
      
      Qmat[rr,18L]  <- rs_at(QQ, Slices$NewHIV,        t)
      Qmat[rr,19L]  <- rs_at(QQ, Slices$NewTreatments, t)
      r0 <- r0 + 52L
    }
    
    R[[k]] <- as.data.frame(Qmat, check.names = FALSE)
  }
  
  # keep structure similar to original: list of trajectories + parameters
  list(R, Params)
}
