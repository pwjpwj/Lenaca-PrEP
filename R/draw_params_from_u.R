draw_params_from_u <- function(u) {
  stopifnot(length(u) == 28L)
  i <- 1L
  
  # 1) NumHSHVIH ~ same as rnorm(1, 80575, 5000)
  NumHSHVIH <- qnorm(u[i], mean = 80575, sd = 5000); i <- i + 1L
  NumHSH    <- 892955 - NumHSHVIH
  
  # 2) Dirichlet for ACD4, I7CD4, Rest (I8CD4 merged into ACD4: 0.20+0.221=0.421)
  alpha1 <- c(0.421, 0.22, 0.359)
  j      <- rdirichlet_u(u[i:(i+2L)], alpha1*10); i <- i + 3L
  ACD4   <- j[1]; I7CD4 <- j[2]; Rest <- j[3]
  
  # 3) Beta-type parameters using invbeta_u (xm,se)
  PropTR   <- invbeta_u(u[i],   xm = 0.89,          se = 0.1);          i <- i + 1L
  G        <- invbeta_u(u[i],   xm = 0.001/52,       se = 0.001/52);     i <- i + 1L
  Woff     <- 0.0038
  W        <- 0.00192
  
  mVIH235  <- invbeta_u(u[i],   xm = 0.0002246,      se = 0.0002246*0.2);    i <- i + 1L
  mVIH35   <- invbeta_u(u[i],   xm = 0.0000594,      se = 0.0000594*0.2);    i <- i + 1L
  mAIDS    <- invbeta_u(u[i],   xm = 0.0006632,      se = 0.0006632*0.2);    i <- i + 1L
  
  # 4) Cvh, Ch, Cl, Cvl
  Cvh      <- invgamma_u(u[i],  xm = 52.5/52,        se = 52.5/52);      i <- i + 1L
  Ch       <- invbeta_u(u[i],   xm = 9/52,           se = 9/52);         i <- i + 1L
  Cl       <- invbeta_u(u[i],   xm = 2.5/52,         se = 2.5/52);       i <- i + 1L
  Cvl      <- invbeta_u(u[i],   xm = 0.2/52,         se = 0.2/52);       i <- i + 1L
  
  # 5) Dirichlet for nvh, nh, nl, nvl
  alpha2   <- c(0.1, 0.4, 0.3, 0.2)
  j1       <- rdirichlet_u(u[i:(i+3L)], alpha2); i <- i + 4L
  nvh      <- j1[1]; nh <- j1[2]; nl <- j1[3]; nvl <- j1[4]
  
  # 6) B and progression probs I12..I78
  B        <- invbeta_u(u[i],   xm = 0.024,          se = 0.024*0.2);                     i <- i + 1L
  
  I12      <- invbeta_u(u[i],   xm = 1 - exp(-0.25),            se = (1 - exp(-0.25))*0.2);            i <- i + 1L
  I23      <- invbeta_u(u[i],   xm = 1 - exp(-0.25),            se = (1 - exp(-0.25))*0.2);            i <- i + 1L
  I34      <- invbeta_u(u[i],   xm = 1 - exp(-0.25),            se = (1 - exp(-0.25))*0.2);            i <- i + 1L
  I45      <- invbeta_u(u[i],   xm = 1 - exp(-0.25/3),          se = (1 - exp(-0.25/3))*0.2);          i <- i + 1L
  I56      <- invbeta_u(u[i],   xm = 1 - exp(-0.25/3),          se = (1 - exp(-0.25/3))*0.2);          i <- i + 1L
  I67      <- invbeta_u(u[i],   xm = 1 - exp(-0.25/3),          se = (1 - exp(-0.25/3))*0.2);          i <- i + 1L
  I78      <- invbeta_u(u[i],   xm = 1 - exp(-0.25/(52*3)),     se = (1 - exp(-0.25/(52*3)))*0.2);     i <- i + 1L
  
  # 7) PTr (3 treatment probabilities)
  PTr1     <- invbeta_u(u[i],   xm = 0.002644231,    se = 0.002644231); i <- i + 1L
  PTr2     <- invbeta_u(u[i],   xm = 0.001682692,    se = 0.001682692); i <- i + 1L
  PTr3     <- invbeta_u(u[i],   xm = 0.000528846,    se = 0.000528846); i <- i + 1L
  PTr      <- c(PTr1, PTr2, PTr3)
  
  # Pack into Params vector in same order as your original code
  param_vec <- numeric(28L)
  param_vec[1]  <- NumHSHVIH
  param_vec[2]  <- ACD4
  param_vec[3]  <- I7CD4
  param_vec[4]  <- Rest
  param_vec[5]  <- PropTR
  param_vec[6]  <- G
  param_vec[7]  <- mVIH235
  param_vec[8]  <- mVIH35
  param_vec[9]  <- mAIDS
  param_vec[10] <- Cvh
  param_vec[11] <- Ch
  param_vec[12] <- Cl
  param_vec[13] <- Cvl
  param_vec[14] <- nvh
  param_vec[15] <- nh
  param_vec[16] <- nl
  param_vec[17] <- nvl
  param_vec[18] <- B
  param_vec[19] <- I12
  param_vec[20] <- I23
  param_vec[21] <- I34
  param_vec[22] <- I45
  param_vec[23] <- I56
  param_vec[24] <- I67
  param_vec[25] <- I78
  param_vec[26] <- PTr1
  param_vec[27] <- PTr2
  param_vec[28] <- PTr3
  
  list(
    NumHSHVIH = NumHSHVIH,
    NumHSH    = NumHSH,
    ACD4      = ACD4,
    I7CD4     = I7CD4,
    Rest      = Rest,
    PropTR    = PropTR,
    G         = G,
    Woff      = Woff,
    W         = W,
    mVIH235   = mVIH235,
    mVIH35    = mVIH35,
    mAIDS     = mAIDS,
    Cvh       = Cvh,
    Ch        = Ch,
    Cl        = Cl,
    Cvl       = Cvl,
    nvh       = nvh,
    nh        = nh,
    nl        = nl,
    nvl       = nvl,
    B         = B,
    I12       = I12,
    I23       = I23,
    I34       = I34,
    I45       = I45,
    I56       = I56,
    I67       = I67,
    I78       = I78,
    PTr       = PTr,
    param_vec = param_vec
  )
}
