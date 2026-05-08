
# data_generation
# -----------------------------------------------------------------------
# Functions for covariate blending, data generation, and AgD conversion.
# -----------------------------------------------------------------------

# --- 1. Blend covariate schemas and covariance matrices via kappa -------

blend_schema_sigma <- function(
    schema_AB1, Sigma_AB1,
    schema_CB2, Sigma_CB2,
    schema_DB3, Sigma_DB3,
    schema_CD,  Sigma_CD,
    kappa_BC = 0, kappa_BD = 0, kappa_CD = 0
) {
  # BC: kappa_BC=1 -> same as AB, kappa_BC=0 -> same as CB2
  sc_BC <- list(
    C1 = list(type="continuous",
              mean = kappa_BC*schema_AB1$C1$mean + (1-kappa_BC)*schema_CB2$C1$mean,
              sd   = kappa_BC*schema_AB1$C1$sd   + (1-kappa_BC)*schema_CB2$C1$sd),
    C2 = list(type="continuous",
              mean = kappa_BC*schema_AB1$C2$mean + (1-kappa_BC)*schema_CB2$C2$mean,
              sd   = kappa_BC*schema_AB1$C2$sd   + (1-kappa_BC)*schema_CB2$C2$sd),
    C3 = list(type="continuous",
              mean = kappa_BC*schema_AB1$C3$mean + (1-kappa_BC)*schema_CB2$C3$mean,
              sd   = kappa_BC*schema_AB1$C3$sd   + (1-kappa_BC)*schema_CB2$C3$sd)
  )
  S_BC <- kappa_BC*Sigma_AB1 + (1-kappa_BC)*Sigma_CB2
  S_BC <- as.matrix(Matrix::nearPD(S_BC, corr=FALSE)$mat)
  R_BC <- cov2cor(S_BC)
  
  # BD: kappa_BD=1 -> same as AB, kappa_BD=0 -> same as DB3
  sc_BD <- list(
    C1 = list(type="continuous",
              mean = kappa_BD*schema_AB1$C1$mean + (1-kappa_BD)*schema_DB3$C1$mean,
              sd   = kappa_BD*schema_AB1$C1$sd   + (1-kappa_BD)*schema_DB3$C1$sd),
    C2 = list(type="continuous",
              mean = kappa_BD*schema_AB1$C2$mean + (1-kappa_BD)*schema_DB3$C2$mean,
              sd   = kappa_BD*schema_AB1$C2$sd   + (1-kappa_BD)*schema_DB3$C2$sd),
    C3 = list(type="continuous",
              mean = kappa_BD*schema_AB1$C3$mean + (1-kappa_BD)*schema_DB3$C3$mean,
              sd   = kappa_BD*schema_AB1$C3$sd   + (1-kappa_BD)*schema_DB3$C3$sd)
  )
  S_BD <- kappa_BD*Sigma_AB1 + (1-kappa_BD)*Sigma_DB3
  S_BD <- as.matrix(Matrix::nearPD(S_BD, corr=FALSE)$mat)
  R_BD <- cov2cor(S_BD)
  
  # CD: kappa_CD=1 -> same as AB, kappa_CD=0 -> same as CD baseline
  sc_CD <- list(
    C1 = list(type="continuous",
              mean = kappa_CD*schema_AB1$C1$mean + (1-kappa_CD)*schema_CD$C1$mean,
              sd   = kappa_CD*schema_AB1$C1$sd   + (1-kappa_CD)*schema_CD$C1$sd),
    C2 = list(type="continuous",
              mean = kappa_CD*schema_AB1$C2$mean + (1-kappa_CD)*schema_CD$C2$mean,
              sd   = kappa_CD*schema_AB1$C2$sd   + (1-kappa_CD)*schema_CD$C2$sd),
    C3 = list(type="continuous",
              mean = kappa_CD*schema_AB1$C3$mean + (1-kappa_CD)*schema_CD$C3$mean,
              sd   = kappa_CD*schema_AB1$C3$sd   + (1-kappa_CD)*schema_CD$C3$sd)
  )
  S_CD <- kappa_CD*Sigma_AB1 + (1-kappa_CD)*Sigma_CD
  S_CD <- as.matrix(Matrix::nearPD(S_CD, corr=FALSE)$mat)
  R_CD <- cov2cor(S_CD)
  
  list(schema_BC=sc_BC, R_BC=R_BC,
       schema_BD=sc_BD, R_BD=R_BD,
       schema_CD=sc_CD, R_CD=R_CD)
}


# --- 2. Generate covariates via Gaussian copula -------------------------

generate_covariates <- function(N, schema, R) {
  Z <- mvtnorm::rmvnorm(N, sigma = R)
  colnames(Z) <- names(schema)
  tibble::as_tibble(setNames(lapply(names(schema), function(v) {
    spec <- schema[[v]]; z <- Z[, v]
    if (spec$type == "binary") as.integer(z > qnorm(1 - spec$p))
    else spec$mean + spec$sd * z
  }), names(schema)))
}


# --- 3. Calibrate beta_0 to match target event rate in control arm -----

# is_active: whether the control arm is an active treatment (not the reference B).
# AB/BC/BD, control arm is B (reference) -> no interaction term in beta_0. (is_active<-F)
# CD, control arm is C (active) -> interaction term included in beta_0 (is_active<-T)

calibrate_intercept <- function(
    trial, X,
    p_A, p_B_AB, p_B_CB, p_C, p_B_DB, p_D, p_C_DC,
    em_scale, coef_common
) {
  N_per_arm <- nrow(X) / 2
  idx <- 1:N_per_arm
  
  if (trial == "AB") {
    target_p <- p_B_AB; is_active <- FALSE
  } else if (trial == "BC") {
    target_p <- p_B_CB; is_active <- FALSE
  } else if (trial == "BD") {
    target_p <- p_B_DB; is_active <- FALSE
  } else if (trial == "CD") {
    target_p <- p_C_DC; is_active <- TRUE   # calibrate on C arm
  } else {
    stop("Unknown trial type: ", trial)
  }
  
  X_means <- colMeans(X[idx, , drop = FALSE])
  covariate_contrib <- coef_common$g_C1 * X_means["C1"] +
    coef_common$g_C2 * X_means["C2"] +
    coef_common$g_C3 * X_means["C3"]
  
  if (is_active) {
    covariate_contrib <- covariate_contrib +
      em_scale * coef_common$delta * X_means["C1"]
  }
  
  beta_0 <- qlogis(target_p) - covariate_contrib
  return(beta_0)
}


# --- 4. Generate IPD for a single study (trial: AB / BC / BD / CD) ---

simulate_one_study <- function(
    N_per_arm,
    trial      = c("AB","BC","BD","CD"),
    effect_class  = "high",
    true_logOR_CB = -0.22,
    consistency   = TRUE,
    p_A, p_B_AB, p_B_CB, p_C, p_B_DB, p_D, p_C_DC,
    em_scale   = 0,
    R, schema,
    coef_common
) {
  trial <- match.arg(trial)
  
  # Set arms and treatment log OR per trial type
  if (trial == "AB") {
    arms     <- c("B", "A")
    beta_trt <- true_logOR_AB(effect_class)
  } else if (trial == "BC") {
    arms     <- c("B", "C")
    beta_trt <- true_logOR_CB
  } else if (trial == "BD") {
    arms     <- c("B", "D")
    beta_trt <- true_logOR_DB(effect_class)
  } else {   # CD
    arms     <- c("C", "D")
    beta_trt <- true_logOR_DC(effect_class, consistency)
  }
  
  # Generate covariates (Gaussian copula) and randomise arms
  X          <- generate_covariates(2*N_per_arm, schema, R)
  X          <- X[sample(1:(2*N_per_arm)), ]   # mimic RCT randomisation
  T_vec      <- rep(arms, each = N_per_arm)
  
  # Calibrate intercept on control arm
  beta_0 <- calibrate_intercept(
    trial, X,
    p_A, p_B_AB, p_B_CB, p_C, p_B_DB, p_D, p_C_DC,
    em_scale, coef_common
  )
  
  # Compute event probabilities
  delta_eff        <- coef_common$delta * em_scale
  covariate_linear <- with(X,
                           coef_common$g_C1 * C1 + coef_common$g_C2 * C2 + coef_common$g_C3 * C3)
  
  logit_p <- ifelse(
    T_vec == arms[1], # control arm(B)
    beta_0 + covariate_linear,
    beta_0 + beta_trt + covariate_linear + delta_eff * X$C1   # Treatment arm
  )
  
  p <- plogis(logit_p)
  y <- rbinom(2*N_per_arm, 1, p)
  
  result <- tibble::tibble(studyn = NA_integer_, trtn = T_vec, y = y) %>%
    dplyr::bind_cols(X)
  attr(result, "beta_0") <- beta_0
  return(result)
}


# --- 5. Aggregate IPD to AgD -------------------------------------------

ipd_to_agd <- function(ipd_df, studyn) {
  ipd_df %>%
    dplyr::mutate(studyn = studyn) %>%
    dplyr::group_by(studyn, trtn) %>%
    dplyr::summarise(
      r       = sum(y),
      N       = dplyr::n(),
      mean_C1 = mean(C1), sd_C1 = sd(C1),
      mean_C2 = mean(C2), sd_C2 = sd(C2),
      mean_C3 = mean(C3), sd_C3 = sd(C3),
      .groups = "drop"
    )
}


# --- 6. Generate full dataset for one scenario -------------------------
#
# Returns IPD and AgD for all comparison arms (AB, BC, BD, CD).
# Note: p_D_DC is not passed — D's event rate in the CD study is
#       implicitly determined by beta_0 (calibrated on C) and logOR_DC.

generate_dataset_AB_BC <- function(scn, covdefs = define_covariates()) {
  set.seed(scn$seed)
  
  # --- AB (single study, uses AB1 distribution) ---
  S_AB <- as.matrix(Matrix::nearPD(covdefs$Sigma_AB1, corr=FALSE)$mat)
  R_AB <- cov2cor(S_AB)
  
  ipd_AB <- simulate_one_study(
    scn$N_per_arm_AB, "AB",
    scn$effect_class, scn$true_logOR_CB, scn$consistency,
    scn$p_A, scn$p_B_AB, scn$p_B_CB, scn$p_C,
    scn$p_B_DB, scn$p_D, scn$p_C_DC,
    scn$em_scale, R_AB, covdefs$schema_AB1, covdefs$coef_common
  ) %>% dplyr::mutate(studyn = 1L)
  beta_0_AB <- attr(ipd_AB, "beta_0")
  agd_AB    <- ipd_to_agd(ipd_AB, studyn = 1L)
  
  # --- BC (K studies, blended distribution via kappa_BC) ---
  mix_BC <- blend_schema_sigma(
    covdefs$schema_AB1, covdefs$Sigma_AB1,
    covdefs$schema_CB2, covdefs$Sigma_CB2,
    covdefs$schema_DB3, covdefs$Sigma_DB3,
    covdefs$schema_CD,  covdefs$Sigma_CD,
    kappa_BC = scn$kappa_BC, kappa_BD = 0, kappa_CD = 0
  )
  
  beta_0_BC_list <- numeric(scn$n_BC)
  ipd_BC_list    <- vector("list", scn$n_BC)
  agd_BC_list    <- vector("list", scn$n_BC)
  for (i in seq_len(scn$n_BC)) {
    sid   <- 100L + i
    ipd_i <- simulate_one_study(
      scn$N_per_arm_BC, "BC",
      scn$effect_class, scn$true_logOR_CB, scn$consistency,
      scn$p_A, scn$p_B_AB, scn$p_B_CB, scn$p_C,
      scn$p_B_DB, scn$p_D, scn$p_C_DC,
      scn$em_scale, mix_BC$R_BC, mix_BC$schema_BC, covdefs$coef_common
    ) %>% dplyr::mutate(studyn = sid)
    beta_0_BC_list[i] <- attr(ipd_i, "beta_0")
    agd_BC_list[[i]]  <- ipd_to_agd(ipd_i, studyn = sid)
    ipd_BC_list[[i]]  <- ipd_i
  }
  ipd_BC <- dplyr::bind_rows(ipd_BC_list)
  agd_BC <- dplyr::bind_rows(agd_BC_list)
  
  # --- BD (K studies, blended distribution via kappa_BD) ---
  mix_BD <- blend_schema_sigma(
    covdefs$schema_AB1, covdefs$Sigma_AB1,
    covdefs$schema_CB2, covdefs$Sigma_CB2,
    covdefs$schema_DB3, covdefs$Sigma_DB3,
    covdefs$schema_CD,  covdefs$Sigma_CD,
    kappa_BC = 0, kappa_BD = scn$kappa_BD, kappa_CD = 0
  )
  
  beta_0_BD_list <- numeric(scn$n_BD)
  ipd_BD_list    <- vector("list", scn$n_BD)
  agd_BD_list    <- vector("list", scn$n_BD)
  for (i in seq_len(scn$n_BD)) {
    sid   <- 200L + i
    ipd_i <- simulate_one_study(
      scn$N_per_arm_BD, "BD",
      scn$effect_class, true_logOR_DB(scn$effect_class), scn$consistency,
      scn$p_A, scn$p_B_AB, scn$p_B_CB, scn$p_C,
      scn$p_B_DB, scn$p_D, scn$p_C_DC,
      scn$em_scale, mix_BD$R_BD, mix_BD$schema_BD, covdefs$coef_common
    ) %>% dplyr::mutate(studyn = sid)
    beta_0_BD_list[i] <- attr(ipd_i, "beta_0")
    agd_BD_list[[i]]  <- ipd_to_agd(ipd_i, studyn = sid)
    ipd_BD_list[[i]]  <- ipd_i
  }
  ipd_BD <- dplyr::bind_rows(ipd_BD_list)
  agd_BD <- dplyr::bind_rows(agd_BD_list)
  
  # --- CD (K studies, blended distribution via kappa_CD) ---
  mix_CD <- blend_schema_sigma(
    covdefs$schema_AB1, covdefs$Sigma_AB1,
    covdefs$schema_CB2, covdefs$Sigma_CB2,
    covdefs$schema_DB3, covdefs$Sigma_DB3,
    covdefs$schema_CD,  covdefs$Sigma_CD,
    kappa_BC = 0, kappa_BD = 0, kappa_CD = scn$kappa_CD
  )
  
  beta_0_CD_list <- numeric(scn$n_CD)
  ipd_CD_list    <- vector("list", scn$n_CD)
  agd_CD_list    <- vector("list", scn$n_CD)
  for (i in seq_len(scn$n_CD)) {
    sid   <- 300L + i
    ipd_i <- simulate_one_study(
      scn$N_per_arm_CD, "CD",
      scn$effect_class, true_logOR_DC(scn$effect_class, scn$consistency),
      scn$consistency,
      scn$p_A, scn$p_B_AB, scn$p_B_CB, scn$p_C,
      scn$p_B_DB, scn$p_D, scn$p_C_DC,
      scn$em_scale, mix_CD$R_CD, mix_CD$schema_CD, covdefs$coef_common
    ) %>% dplyr::mutate(studyn = sid)
    beta_0_CD_list[i] <- attr(ipd_i, "beta_0")
    agd_CD_list[[i]]  <- ipd_to_agd(ipd_i, studyn = sid)
    ipd_CD_list[[i]]  <- ipd_i
  }
  ipd_CD <- dplyr::bind_rows(ipd_CD_list)
  agd_CD <- dplyr::bind_rows(agd_CD_list)
  
  list(
    config   = scn,
    ipd_AB   = ipd_AB,   ipd_BC = ipd_BC,
    ipd_BD   = ipd_BD,   ipd_CD = ipd_CD,
    agd_AB   = agd_AB,   agd_BC = agd_BC,
    agd_BD   = agd_BD,   agd_CD = agd_CD,
    beta_0_AB = beta_0_AB,
    beta_0_BC = mean(beta_0_BC_list),
    beta_0_BD = if (scn$n_BD > 0) mean(beta_0_BD_list) else NA_real_,
    beta_0_CD = if (scn$n_CD > 0) mean(beta_0_CD_list) else NA_real_
  )
}