sim_weekly_to_yearly <- function(sim_df, years = 5L, weeks_per_year = 52L) {
  # Ensure ordered
  sim_df <- sim_df[order(sim_df$Weeks), ]
  
  # Year index: 1..5
  sim_df$Year <- ceiling(sim_df$Weeks / weeks_per_year)
  
  # 1) New AIDS cases: sum weekly counts per year
  newAIDS_year <- tapply(sim_df$NewAIDS, sim_df$Year, sum)
  
  # 2) New HIV infections: sum weekly New_HIV_inf per year
  newInf_year  <- tapply(sim_df$New_HIV_inf, sim_df$Year, sum)
  
  # 3) AIDS deaths:
  #    'Dead' is a flow (weekly death count), not a stock. Sum directly per year.
  AIDSdeaths_year <- tapply(sim_df$Dead, sim_df$Year, sum)
  
  # Build tidy yearly data.frame
  data.frame(
    Year           = as.integer(names(newAIDS_year)),
    NewAIDS_sim    = as.numeric(newAIDS_year),
    AIDSdeaths_sim = as.numeric(AIDSdeaths_year),
    NewInf_sim     = as.numeric(newInf_year)
  )
}
