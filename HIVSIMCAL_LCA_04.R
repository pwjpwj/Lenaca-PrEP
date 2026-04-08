HIVSIMCAL <- function(repeats, m, DistEdad, RiskRed) {
  # ---- constants 
  Age   <- 85L
  years <- 5L
  NoSem <- 52L
  
  # Predefine column slices once (fixed model layout)
  make_slices_HIVSIM <- function(Age){
    stopifnot(Age > 0)
    
    w   <- 1L
    blk <- 6L * Age          # 4 non-PrEP + 2 PrEP per state
    st  <- function(s) {     # s = 1..10 (S0, I1..I8, AIDS)
      start <- 2L + (s-1L)*blk
      start:(start + blk - 1L)
    }
    
    slices <- list(
      Weeks = 1L,
      S0    = st(1),
      I1    = st(2),
      I2    = st(3),
      I3    = st(4),
      I4    = st(5),
      I5    = st(6),
      I6    = st(7),
      I7    = st(8),
      I8    = st(9),
      AIDS  = st(10)
    )
    
    after_states <- 1L + 10L*blk  # column index of last state block end
    takeA <- function(offset){ (after_states + offset*Age + 1L):(after_states + (offset+1L)*Age) }
    
    slices$Treated       <- takeA(0)
    slices$Dead          <- takeA(1)
    slices$NewAIDS       <- takeA(2)  # CalAIDS
    slices$FI_vh         <- takeA(3)
    slices$FI_h          <- takeA(4)
    slices$FI_l          <- takeA(5)
    slices$FI_vl         <- takeA(6)
    slices$NewHIV        <- takeA(7)
    slices$NewTreatments <- takeA(8)
    
    slices$ncol_expected <- 1L + 10L*blk + 9L*Age
    slices
  }
  
  # ---- outputs ----
  R      <- vector("list", repeats)
  Params <- matrix(0, nrow = 32L, ncol = repeats)
  
  # small helper: fast row-sums for a slice at time t
  rs_at <- function(QQ, idx, t) {
    if (length(idx) == 1L) QQ[, idx, t]
    else rowSums(QQ[, idx, t, drop = FALSE])
  }
  
  # ---- main loop over simulations ----
  for (k in seq_len(repeats)) {
    # --- parameter draws (unchanged logic; keep your priors/distributions) ---
    NumHSHVIH <- rnorm(1, 80575, 5000);           Params[1, k] <- NumHSHVIH
    NumHSH    <- 892955 - NumHSHVIH
    
    j <- rdirichlet(1, c(0.20, 0.221, 0.22, 0.359))
    ACD4 <- j[1, 1]; Params[2, k] <- ACD4
    I8CD4<- j[1, 2]; Params[3, k] <- I8CD4
    I7CD4<- j[1, 3]; Params[4, k] <- I7CD4
    Rest <- j[1, 4]; Params[5, k] <- Rest
    
    PropTR <- invbeta(0.807, 0.1);                 Params[6,  k] <- PropTR
    G      <- invbeta(0.001/52, 0.001/52);        Params[7,  k] <- G
    Woff <- 0.014/52;  W <- 0.089/52
    
    mVIH235 <- invbeta(0.0002246, 0.0002246);     Params[8,  k] <- mVIH235
    mVIH35  <- invbeta(0.0000594, 0.0000594);     Params[9,  k] <- mVIH35
    mAIDS   <- invbeta(0.0006632, 0.0006632);     Params[10, k] <- mAIDS
    
    Cvh <- invgamma(52.5/52, 52.5/52);            Params[11, k] <- Cvh
    Ch  <- invbeta(9/52,   9/52);                 Params[12, k] <- Ch
    Cl  <- invbeta(2.5/52, 2.5/52);               Params[13, k] <- Cl
    Cvl <- invbeta(0.2/52, 0.2/52);               Params[14, k] <- Cvl
    
    j1 <- rdirichlet(1, c(0.1, 0.4, 0.3, 0.2))
    nvh <- j1[1,1]; Params[15, k] <- nvh
    nh  <- j1[1,2]; Params[16, k] <- nh
    nl  <- j1[1,3]; Params[17, k] <- nl
    nvl <- j1[1,4]; Params[18, k] <- nvl
    
    B   <- invbeta(0.024, 0.024);                 Params[19, k] <- B
    I12 <- invbeta(1 - exp(-0.25),             1 - exp(-0.25));             Params[20, k] <- I12
    I23 <- invbeta(1 - exp(-0.25),             1 - exp(-0.25));             Params[21, k] <- I23
    I34 <- invbeta(1 - exp(-0.25),             1 - exp(-0.25));             Params[22, k] <- I34
    I45 <- invbeta(1 - exp(-0.25/3),           1 - exp(-0.25/3));           Params[23, k] <- I45
    I55 <- invbeta(1 - exp(-0.25/3),           1 - exp(-0.25/3));           Params[25, k] <- I55
    I56 <- invbeta(1 - exp(-0.25/3),           1 - exp(-0.25/3));           Params[26, k] <- I56
    I67 <- invbeta(1 - exp(-0.25/3),           1 - exp(-0.25/3));           Params[27, k] <- I67
    I78 <- invbeta(1 - exp(-0.25/(52*3)),      1 - exp(-0.25/(52*3)));      Params[28, k] <- I78
    I89 <- invbeta(1 - exp(-0.25/(52*3.5)),    1 - exp(-0.25/(52*3.5)));    Params[29, k] <- I89
    
    PTr <- c(
      invbeta(0.002644231, 0.002644231),
      invbeta(0.001682692, 0.001682692),
      invbeta(0.000528846, 0.000528846)
    )
    Params[30, k] <- PTr[1]
    Params[31, k] <- PTr[2]
    Params[32, k] <- PTr[3]
    
    # --- simulate ---
    Q <- HIVSIM(
      Age, NoSem, NumHSH, NumHSHVIH,
      ACD4, I8CD4, I7CD4, Rest,
      G, Woff, W,
      m, mVIH235, mVIH35, mAIDS,    
      Cvh, Ch, Cl, Cvl,
      nvh, nh, nl, nvl, B,
      I12, I23, I34,
      I45, I55,
      I56, I67,
      I78, I89, PTr, years, PropTR, DistEdad,
      RiskRed
    )
    
    Slices <- make_slices_HIVSIM(Age)
    stopifnot(Slices$ncol_expected == 5866L)
    
    QQ <- QQ <- array(
      as.numeric(unlist(lapply(Q, as.matrix))),
      dim = c(52L, Slices$ncol_expected, years)
    )
    
    stopifnot(dim(QQ)[2] == Slices$ncol_expected)
    
    # --- build output matrix  ---
    nrows <- years * 52L
    Qmat  <- matrix(0, nrow = nrows, ncol = 20L)
    colnames(Qmat) <- c(
      "Weeks", "Susceptible", "Inf1", "Inf2", "Inf3", "Inf4",
      "Inf5", "Inf6", "Inf7", "Inf8", "AIDS", "Treated", "Dead", "NewAIDS",
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
      Qmat[rr,10L]  <- rs_at(QQ, Slices$I8,   t)
      Qmat[rr,11L]  <- rs_at(QQ, Slices$AIDS, t)
      
      Qmat[rr,12L]  <- rs_at(QQ, Slices$Treated, t)
      Qmat[rr,13L]  <- rs_at(QQ, Slices$Dead,    t)
      Qmat[rr,14L]  <- rs_at(QQ, Slices$NewAIDS, t)
      
      Qmat[rr,15L]  <- rs_at(QQ, Slices$FI_vh, t)
      Qmat[rr,16L]  <- rs_at(QQ, Slices$FI_h,  t)
      Qmat[rr,17L]  <- rs_at(QQ, Slices$FI_l,  t)
      Qmat[rr,18L]  <- rs_at(QQ, Slices$FI_vl, t)
      
      Qmat[rr,19L]  <- rs_at(QQ, Slices$NewHIV,        t)
      Qmat[rr,20L]  <- rs_at(QQ, Slices$NewTreatments, t)
      r0 <- r0 + 52L
    }
    
    # convert once at the end
    Q1 <- as.data.frame(Qmat, check.names = FALSE)
    R[[k]] <- Q1
  }
  
  list(R, Params)
}
