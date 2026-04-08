make_slices_HIVSIM <- function(Age){
  stopifnot(Age > 0)
  
  blk <- 6L * Age          # 4 non-PrEP + 2 PrEP per state
  st  <- function(s) {     # s = 1..9 (S0, I1..I7, AIDS)
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
    AIDS  = st(9)
  )
  
  after_states <- 1L + 9L*blk  # column index of last state block end
  takeA <- function(offset){
    (after_states + offset*Age + 1L):(after_states + (offset+1L)*Age)
  }
  
  slices$Treated       <- takeA(0)
  slices$Dead          <- takeA(1)
  slices$NewAIDS       <- takeA(2)  # CalAIDS
  slices$FI_vh         <- takeA(3)
  slices$FI_h          <- takeA(4)
  slices$FI_l          <- takeA(5)
  slices$FI_vl         <- takeA(6)
  slices$NewHIV        <- takeA(7)
  slices$NewTreatments <- takeA(8)
  
  slices$ncol_expected <- 1L + 9L*blk + 9L*Age
  slices
}
