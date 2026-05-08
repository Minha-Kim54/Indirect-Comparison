define_covariates <- function() {
  schema_AB1 <- list(                 
    C1 = list(type="continuous", mean=55, sd=15),
    C2 = list(type="continuous", mean=125, sd=5),
    C3 = list(type="continuous", mean=35, sd=10)
  )
  Sigma_AB1 <- matrix(c(              
    1.0, 0.1, 0.4,
    0.1, 1.0, 0.2,
    0.4, 0.2, 1.0
  ), 3, 3, byrow=TRUE,
  dimnames=list(names(schema_AB1), names(schema_AB1)))
  
  schema_CB2 <- list(                 
    C1 = list(type="continuous", mean=65, sd=10),
    C2 = list(type="continuous", mean=135, sd=15),
    C3 = list(type="continuous", mean=25, sd=5)
  )
  Sigma_CB2 <- matrix(c(
    1.0, 0.1, 0.4,
    0.1, 1.0, 0.2,
    0.4, 0.2, 1.0
  ), 3, 3, byrow=TRUE,
  dimnames=list(names(schema_CB2), names(schema_CB2)))
  
  schema_DB3 <- list(                
    C1 = list(type="continuous", mean=65, sd=10),
    C2 = list(type="continuous", mean=135, sd=15),
    C3 = list(type="continuous", mean=25, sd=5)
  )
  Sigma_DB3 <- matrix(c(
    1.0, 0.1, 0.4,
    0.1, 1.0, 0.2,
    0.4, 0.2, 1.0
  ), 3, 3, byrow=TRUE,
  dimnames=list(names(schema_DB3), names(schema_DB3)))
  
  schema_CD <- list(                
    C1 = list(type="continuous", mean=65, sd=10),
    C2 = list(type="continuous", mean=135, sd=15),
    C3 = list(type="continuous", mean=25, sd=5)
  )
  Sigma_CD <- matrix(c(
    1.0, 0.1, 0.4,
    0.1, 1.0, 0.2,
    0.4, 0.2, 1.0
  ), 3, 3, byrow=TRUE,
  dimnames=list(names(schema_CD), names(schema_CD)))
  

  coef_common <- list(
    g_C1 =  0.01,
    g_C2 = 0.01,
    g_C3 = 0.01,
    delta = -0.02
  )
  
  list(schema_AB1=schema_AB1, Sigma_AB1=Sigma_AB1,
       schema_CB2=schema_CB2, Sigma_CB2=Sigma_CB2,
       schema_DB3=schema_DB3, Sigma_DB3=Sigma_DB3,
       schema_CD=schema_CD, Sigma_CD=Sigma_CD,
       coef_common=coef_common)
}

# True log OR mapping for AB (based on Weber Table 3)
true_logOR_AB <- function(effect_class = c("high","moderate","low","no")){
  effect_class <- match.arg(effect_class)
  switch(effect_class,
         high     = -0.7,
         moderate = -0.45,
         low      = -0.28,
         no       = -0.22)
}

# True log OR mapping for DB
true_logOR_DB <- function(effect_class = c("high", "moderate", "low", "no")){
  effect_class <- match.arg(effect_class)
  switch(effect_class,
         high     = -0.44,
         moderate = -0.3,
         low      = -0.25,
         no       = -0.22)
}

# True log OR mapping for DC
# Under consistency: derived as logOR_DB - logOR_CB
# Under inconsistency: fixed at -0.42 across all effect classes
true_logOR_DC <- function(effect_class = c("high", "moderate", "low", "no"),
                          consistency=TRUE){
  effect_class <- match.arg(effect_class)
  
  if (consistency){
    true_logOR_DB(effect_class) - (-0.22)
  } else {
    switch(effect_class,
           high     = -0.42,
           moderate = -0.42,
           low      = -0.42,
           no       = -0.42)
  }
}

# Construct a single scenario object (called internally by scenario_grid())
make_scenario <- function(
    n_BC = 3,
    n_BD = 3,
    n_CD = 3,
    N_per_arm_AB = c(200,500),
    N_per_arm_BC = c(200,500),
    N_per_arm_BD = c(200,500),
    N_per_arm_CD = c(200,500),
    effect_class = c("high","moderate","low","no"),
    true_logOR_CB = -0.22,
    # Target event probabilities (see Table 4 in paper)
    p_A = 0.22,       # Probability in arm A (AB study)
    p_B_AB = 0.36,    # Probability in arm B (AB study)
    p_B_CB = 0.35,    # Probability in arm B (BC study)
    p_C = 0.30,       # Probability in arm C (BC study)
    p_B_DB = 0.40,    # Probability in arm B (BD study)
    p_D = 0.30,       # Probability in arm D (BD study)
    p_C_DC_inconsistent = 0.34,  # Probability in arm C (DC study, inconsistency scenario)
    # p_D_DC_inconsistent is derived from p_C_DC_inconsistent and logOR_DC = -0.42
    # p_D = plogis(qlogis(p_C_DC_inconsistent) + (-0.42))
    consistency = TRUE,
    kappa_BC = 0,
    kappa_BD = 0,
    kappa_CD = 0,
    em_scale = c(0,1),
    seed = 1000
){
  # In CD studies, beta_0 is calibrated on the C arm (control).
  # D's event rate is then determined automatically via beta_0 + logOR_DC.
  # - Consistent:   p_C_DC = p_C (0.30, same as BC study, see Table 4)
  # - Inconsistent: p_C_DC = 0.34 (see Table 4)
  p_C_DC <- if (consistency) p_C else p_C_DC_inconsistent
  
  effect_class <- match.arg(effect_class)
  list(
    n_BC = n_BC,
    n_BD = n_BD,
    n_CD = n_CD,
    N_per_arm_AB = N_per_arm_AB,
    N_per_arm_BC = N_per_arm_BC,
    N_per_arm_BD = N_per_arm_BD,
    N_per_arm_CD = N_per_arm_CD,
    effect_class = effect_class,
    consistency = consistency,
    true_logOR_CB = true_logOR_CB,
    p_A = p_A,
    p_B_AB = p_B_AB,
    p_B_CB = p_B_CB,
    p_C = p_C,
    p_B_DB = p_B_DB,
    p_D = p_D,
    p_C_DC = p_C_DC,
    kappa_BC = kappa_BC,
    kappa_BD = kappa_BD,
    kappa_CD = kappa_CD,
    em_scale = em_scale,
    seed = seed
  )
}


# -----------------------------------------------------------------------
# Scenario grid
# -----------------------------------------------------------------------

# 1) Scenario grid constructor
scenario_grid <- function(
    network_type = c("basic_3", "star_4", "one_loop_closed"),
    n_studies       = c(1, 2, 3, 5, 7, 10, 15, 20),
    consistency     = c(TRUE, FALSE),
    N_per_arm_AB    = c(200, 500),
    N_per_arm       = c(200, 500),
    effect_class    = c("high"),
    true_logOR_CB   = -0.22,
    p_A             = 0.22,
    p_B_AB          = 0.36,
    p_B_CB          = 0.35,
    p_C             = 0.30,
    p_B_DB          = 0.40,
    p_D             = 0.30,
    p_C_DC_inconsistent = 0.34,  # p_D_DC is derived internally in make_scenario()
    kappa           = c(0, 1),
    em_scale        = c(0, 1),
    base_seed       = 1000
){
  
  network_type <- match.arg(network_type)
  
  # Determine number of studies per comparison arm based on network type
  if (network_type == "basic_3") {
    studies_grid <- data.frame(
      n_BC = n_studies,
      n_BD = 0,
      n_CD = 0
    )
  } else if (network_type == "star_4") {
    studies_grid <- data.frame(
      n_BC = n_studies,
      n_BD = n_studies,
      n_CD = 0
    )
  } else if (network_type == "one_loop_closed") {
    studies_grid <- data.frame(
      n_BC = n_studies,
      n_BD = n_studies,
      n_CD = n_studies
    )
  }
  
  other_grid <- expand.grid(
    consistency  = consistency,
    N_per_arm_AB = N_per_arm_AB,
    N_per_arm    = N_per_arm,
    effect_class = effect_class,
    true_logOR_CB = true_logOR_CB,
    kappa        = kappa,
    em_scale     = em_scale,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  
  # Cross-join studies_grid and other_grid
  grid <- do.call(rbind, lapply(1:nrow(studies_grid), function(i) {
    cbind(studies_grid[i, , drop = FALSE], other_grid)
  }))
  
  scenarios <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    gi <- grid[i, ]
    scenarios[[i]] <- make_scenario(
      n_BC = gi$n_BC,
      n_BD = gi$n_BD,
      n_CD = gi$n_CD,
      consistency          = gi$consistency,
      N_per_arm_AB         = gi$N_per_arm_AB,
      N_per_arm_BC         = gi$N_per_arm,
      N_per_arm_BD         = gi$N_per_arm,
      N_per_arm_CD         = gi$N_per_arm,
      effect_class         = gi$effect_class,
      true_logOR_CB        = gi$true_logOR_CB,
      p_A                  = p_A,
      p_B_AB               = p_B_AB,
      p_B_CB               = p_B_CB,
      p_C                  = p_C,
      p_B_DB               = p_B_DB,
      p_D                  = p_D,
      p_C_DC_inconsistent  = p_C_DC_inconsistent,
      kappa_BC             = gi$kappa,
      kappa_BD             = gi$kappa,
      kappa_CD             = gi$kappa,
      em_scale             = gi$em_scale,
      seed                 = base_seed + i
    )
  }
  
  # Auto-generate descriptive scenario names
  names(scenarios) <- sprintf(
    "네트워크:%s_연구:%d_표본_AB:%d_표본_BC/BD/CD:%d_Effect_size:%s_Kappa:%.1f_EM:%d_Con%s",
    network_type,
    grid$n_BC,
    grid$N_per_arm_AB,
    grid$N_per_arm,
    substr(grid$effect_class, 1, 1),
    grid$kappa,
    grid$em_scale,
    ifelse(grid$consistency, "T", "F")
  )
  attr(scenarios, "grid") <- grid
  scenarios
}

# 2) Batch dataset generation across all scenarios (data generation only)
run_dataset_grid <- function(scenarios, covdefs = define_covariates){
  out <- lapply(scenarios, function(scn) generate_dataset_AB_BC(scn, covdefs = covdefs))
  # Summary config table
  cfgs <- do.call(rbind, lapply(out, function(x) as.data.frame(t(unlist(x$config)))))
  rownames(cfgs) <- names(out)
  list(datasets = out, config_table = cfgs)
}