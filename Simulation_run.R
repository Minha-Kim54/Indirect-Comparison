### 1. Single replication ###
run_rep_core <- function(r, scenarios, scen_names, base_seed, jags_ctrl, covdefs, include_bucher, include_bayes) {
  
  reject_null <- function(LCL, UCL) as.integer((LCL > 0) | (UCL < 0))
  
  # ===== Effect modifier가 있을 때 marginal logOR_AC 계산 =====
  calculate_marginal_logOR_AC <- function(ipd_AB, ipd_BC, covdefs, effect_class,true_logOR_CB, em_scale) {
    coef <- covdefs$coef_common
    
    mean_C1_AB <- mean(ipd_AB$C1)
    mean_C1_BC <- mean(ipd_BC$C1)
    
    # Marginal effect
    marginal_effect_AB <- true_logOR_AB(effect_class) + coef$delta * em_scale * mean_C1_AB
    marginal_effect_CB <- true_logOR_CB + coef$delta * em_scale * mean_C1_BC
    
    marginal_logOR_AC <- marginal_effect_AB - marginal_effect_CB
    
    return(marginal_logOR_AC)
  }
  
  rows <- vector("list", length(scenarios))
  for (j in seq_along(scenarios)) {
    scn <- scenarios[[j]]
    # AB, BC, effect_class, kappa_BC, em_scale만으로 seed 결정
    effect_num <- switch(scn$effect_class,
                         "high" = 1,
                         "moderate" = 2, 
                         "low" = 3,
                         "no" = 4)
    
    seed_offset <- scn$N_per_arm_AB + 
      scn$N_per_arm_BC * 10 + 
      scn$n_BC * 100 + 
      scn$n_BD * 150 + 
      scn$n_CD * 200 + 
      scn$kappa_BC * 100 + 
      scn$em_scale * 100 +
      as.integer(scn$consistency) * 250 +
      effect_num * 1000
    
    scn$seed <- base_seed + r*1000 + seed_offset
    jags_ctrl$seed <- base_seed + r*100 + j
    #Generate dataset
    ds <- generate_dataset_AB_BC(scn, covdefs = covdefs)
    
    beta_0_info<-tibble(
      beta_0_AB = ds$beta_0_AB,
      beta_0_BC = ds$beta_0_BC,
      beta_0_BD = ds$beta_0_BD,
      beta_0_CD = ds$beta_0_CD
    )
    
    # true indirect effect AC
    true_AC <- true_logOR_AB(scn$effect_class)-scn$true_logOR_CB
    
    # Covariate coefficient estimates (logistic regression on IPD)
    cov_summary <- tibble(
      seed_used = scn$seed,
      #AB study
      est_C1_coef_AB = {
        ds_temp <- ds$ipd_AB %>%
          mutate(
            C1 = C1,
            C2 = C2,
            C3 = C3
          )
        fit_AB <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_AB)["C1"]
      },
      est_C2_coef_AB = {
        ds_temp <- ds$ipd_AB %>%
          mutate(
            C1_std = C1,
            C2_std = C2,
            C3_std = C3
          )
        fit_AB <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_AB)["C2"]
      },
      est_C3_coef_AB = {
        ds_temp <- ds$ipd_AB %>%
          mutate(
            C1 = C1,
            C2 = C2,
            C3 = C3
          )
        fit_AB <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_AB)["C3"]
      },
      
      # BC study
      est_C1_coef_BC = {
        ds_temp <- ds$ipd_BC %>%
          mutate(
            C1 = C1,
            C2 = C2,
            C3 = C3
          )
        fit_BC <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_BC)["C1"]
      },
      est_C2_coef_BC = {
        ds_temp <- ds$ipd_BC %>%
          mutate(
            C1 = C1,
            C2 = C2,
            C3 = C3
          )
        fit_BC <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_BC)["C2"]
      },
      est_C3_coef_BC = {
        ds_temp <- ds$ipd_BC %>%
          mutate(
            C1 = C1,
            C2 = C2,
            C3 = C3
          )
        fit_BC <- glm(y ~ trtn + C1 + C2 + C3, 
                      family = binomial(), data = ds_temp)
        coef(fit_BC)["C3"]
      },
      
      # BD study 
      est_C1_coef_BD = {
        if (nrow(ds$ipd_BD) > 0) {
          ds_temp <- ds$ipd_BD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_BD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_BD)["C1"]
        } else {
          NA_real_
        }
      },
      est_C2_coef_BD = {
        if (nrow(ds$ipd_BD) > 0) {
          ds_temp <- ds$ipd_BD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_BD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_BD)["C2"]
        } else {
          NA_real_
        }
      },
      est_C3_coef_BD = {
        if (nrow(ds$ipd_BD) > 0) {
          ds_temp <- ds$ipd_BD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_BD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_BD)["C3"]
        } else {
          NA_real_
        }
      },
      
      # CD study
      est_C1_coef_CD = {
        if (nrow(ds$ipd_CD) > 0) {
          ds_temp <- ds$ipd_CD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_CD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_CD)["C1"]
        } else {
          NA_real_
        }
      },
      est_C2_coef_CD = {
        if (nrow(ds$ipd_CD) > 0) {
          ds_temp <- ds$ipd_CD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_CD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_CD)["C2"]
        } else {
          NA_real_
        }
      },
      est_C3_coef_CD = {
        if (nrow(ds$ipd_CD) > 0) {
          ds_temp <- ds$ipd_CD %>%
            mutate(C1 = C1, C2 = C2, C3 = C3)
          fit_CD <- glm(y ~ trtn + C1 + C2 + C3, 
                        family = binomial(), data = ds_temp)
          coef(fit_CD)["C3"]
        } else {
          NA_real_
        }
      },
      
      #Covariate means and SD
      AB_mean_C1 = mean(ds$ipd_AB$C1),
      AB_sd_C1 = sd(ds$ipd_AB$C1),
      AB_mean_C2 = mean(ds$ipd_AB$C2),
      AB_sd_C2 = sd(ds$ipd_AB$C2),
      AB_mean_C3 = mean(ds$ipd_AB$C3),
      AB_sd_C3 = sd(ds$ipd_AB$C3),
      
      BC_mean_C1 = mean(ds$ipd_BC$C1),
      BC_sd_C1 = sd(ds$ipd_BC$C1),
      BC_mean_C2 = mean(ds$ipd_BC$C2),
      BC_sd_C2 = sd(ds$ipd_BC$C2),
      BC_mean_C3 = mean(ds$ipd_BC$C3),
      BC_sd_C3 = sd(ds$ipd_BC$C3),
      
      BD_mean_C1 = if (nrow(ds$ipd_BD) > 0) mean(ds$ipd_BD$C1) else NA_real_,
      BD_sd_C1 = if (nrow(ds$ipd_BD) > 0) sd(ds$ipd_BD$C1) else NA_real_,
      BD_mean_C2 = if (nrow(ds$ipd_BD) > 0) mean(ds$ipd_BD$C2) else NA_real_,
      BD_sd_C2 = if (nrow(ds$ipd_BD) > 0) sd(ds$ipd_BD$C2) else NA_real_,
      BD_mean_C3 = if (nrow(ds$ipd_BD) > 0) mean(ds$ipd_BD$C3) else NA_real_,
      BD_sd_C3 = if (nrow(ds$ipd_BD) > 0) sd(ds$ipd_BD$C3) else NA_real_,
      
      CD_mean_C1 = if (nrow(ds$ipd_CD) > 0) mean(ds$ipd_CD$C1) else NA_real_,
      CD_sd_C1 = if (nrow(ds$ipd_CD) > 0) sd(ds$ipd_CD$C1) else NA_real_,
      CD_mean_C2 = if (nrow(ds$ipd_CD) > 0) mean(ds$ipd_CD$C2) else NA_real_,
      CD_sd_C2 = if (nrow(ds$ipd_CD) > 0) sd(ds$ipd_CD$C2) else NA_real_,
      CD_mean_C3 = if (nrow(ds$ipd_CD) > 0) mean(ds$ipd_CD$C3) else NA_real_,
      CD_sd_C3 = if (nrow(ds$ipd_CD) > 0) sd(ds$ipd_CD$C3) else NA_real_
      
    )
    
    # observed event rates per arm
    ab_summary <- ds$agd_AB %>%
      group_by(trtn) %>%
      summarise(total_r = sum(r), total_N = sum(N), .groups = "drop")
    
    pA <- ab_summary$total_r[ab_summary$trtn == "A"] / ab_summary$total_N[ab_summary$trtn == "A"]
    pB_AB <- ab_summary$total_r[ab_summary$trtn == "B"] / ab_summary$total_N[ab_summary$trtn == "B"]
    
    bc_summary <- ds$agd_BC %>%
      group_by(trtn) %>%
      summarise(total_r = sum(r), total_N = sum(N), .groups = "drop")
    
    pB_BC <- bc_summary$total_r[bc_summary$trtn == "B"] / bc_summary$total_N[bc_summary$trtn == "B"]
    pC <- bc_summary$total_r[bc_summary$trtn == "C"] / bc_summary$total_N[bc_summary$trtn == "C"]
    
    if (nrow(ds$agd_BD) > 0) {
      bd_summary <- ds$agd_BD %>%
        group_by(trtn) %>%
        summarise(total_r = sum(r), total_N = sum(N), .groups = "drop")
      
      pB_BD <- bd_summary$total_r[bd_summary$trtn == "B"] / bd_summary$total_N[bd_summary$trtn == "B"]
      pD <- bd_summary$total_r[bd_summary$trtn == "D"] / bd_summary$total_N[bd_summary$trtn == "D"]
    } else {
      pB_BD <- NA_real_
      pD <- NA_real_
    }
    
    if (nrow(ds$agd_CD) > 0) {
      cd_summary <- ds$agd_CD %>%
        group_by(trtn) %>%
        summarise(total_r = sum(r), total_N = sum(N), .groups = "drop")
      
      pC_DC <- cd_summary$total_r[cd_summary$trtn == "C"] / cd_summary$total_N[cd_summary$trtn == "C"]
      pD_DC <- cd_summary$total_r[cd_summary$trtn == "D"] / cd_summary$total_N[cd_summary$trtn == "D"]
    } else {
      pC_DC <- NA_real_
      pD_DC <- NA_real_
    }
    
    all_events <- sum(ds$agd_AB$r) + sum(ds$agd_BC$r)
    all_samples <- sum(ds$agd_AB$N) + sum(ds$agd_BC$N)
    
    if (nrow(ds$agd_BD) > 0) {
      all_events <- all_events + sum(ds$agd_BD$r)
      all_samples <- all_samples + sum(ds$agd_BD$N)
    }
    
    if (nrow(ds$agd_CD) > 0) {
      all_events <- all_events + sum(ds$agd_CD$r)
      all_samples <- all_samples + sum(ds$agd_CD$N)
    }
    
    event_rates_rep <- tibble(
      pA = pA,
      pB_AB = pB_AB,
      pB_BC = pB_BC,
      pC = pC,
      pB_BD = pB_BD,
      pD = pD,
      pC_DC = pC_DC,
      pD_DC = pD_DC,
      overall_rate = all_events / all_samples  # 실제 전체 발생률
    )
    
    
    out_list <- list()
    
    if (include_bucher) {
      b <- do_bucher_agd_agd(ds$agd_AB, ds$agd_BC) %>%
        dplyr::mutate(SEorSD = SE,
                      rep = r, scen_name = scen_names[j],
                      consistency = scn$consistency,
                      n_BC = scn$n_BC, n_BD = scn$n_BD, n_CD = scn$n_CD,
                      N_AB = scn$N_per_arm_AB, N_BC = scn$N_per_arm_BC, N_BD = scn$N_per_arm_BD,
                      N_CD = scn$N_per_arm_CD,
                      effect_class = scn$effect_class, 
                      kappa_BC = scn$kappa_BC, kappa_BD = scn$kappa_BD, kappa_CD = scn$kappa_CD, 
                      em_scale = scn$em_scale,
                      true_AC = true_AC,
                      reject = reject_null(LCL, UCL)) %>%
        bind_cols(event_rates_rep) %>%
        bind_cols(cov_summary)%>%
        bind_cols(beta_0_info)
      out_list$bucher <- b}
    
    if (include_bayes) {
      bay <- run_bayes_nma_gemtc(ds$agd_AB, ds$agd_BC, ds$agd_BD, ds$agd_CD,
                                 n.chains=jags_ctrl$n.chains,
                                 n.adapt=jags_ctrl$n.adapt,
                                 n.iter=jags_ctrl$n.iter,
                                 thin=jags_ctrl$thin,
                                 em_scale=scn$em_scale,
                                 seed=jags_ctrl$seed,
                                 beta_0_AB = ds$beta_0_AB,
                                 beta_0_BC = ds$beta_0_BC,
                                 beta_0_BD = ds$beta_0_BD,
                                 beta_0_CD = ds$beta_0_CD) %>%
        dplyr::rename(SEorSD = SD) %>%
        dplyr::mutate(rep = r, scen_name = scen_names[j],
                      consistency = scn$consistency,
                      n_BC = scn$n_BC, n_BD = scn$n_BD, n_CD = scn$n_CD,
                      N_AB = scn$N_per_arm_AB, N_BC = scn$N_per_arm_BC, N_BD = scn$N_per_arm_BD,
                      N_CD = scn$N_per_arm_CD,
                      effect_class = scn$effect_class, kappa_BC = scn$kappa_BC, kappa_BD = scn$kappa_BD,
                      kappa_CD = scn$kappa_CD, em_scale = scn$em_scale,
                      p_A = scn$p_A, p_B_AB = scn$p_B_AB, p_B_CB = scn$p_B_CB, p_C = scn$p_C, 
                      p_B_DB = scn$p_B_DB, p_D = scn$p_D, p_C_DC = scn$p_C_DC, p_D_DC = scn$p_D_DC,
                      true_AC = true_AC,
                      reject = reject_null(LCL, UCL))%>%
        bind_cols(event_rates_rep)%>%
        bind_cols(cov_summary)%>%
        bind_cols(beta_0_info)
      
      out_list$bayes <- bay
    }
    
    rows[[j]] <- dplyr::bind_rows(out_list)
  }
  dplyr::bind_rows(rows)
}

### 2. Full Monte Carlo simulation ###
run_simulations_min <- function(
    scenarios,
    R = 1000,
    include_bucher = TRUE,
    include_bayes  = TRUE,  # JAGS 베이지안 NMA
    jags_ctrl = list(n.chains=3, n.adapt=3000, n.iter=8000, thin=8, seed=1000),
    base_seed = 1000,
    covdefs = define_covariates(),
    parallel = TRUE, n_cores = max(1, parallel::detectCores()-1)
){
  stopifnot(include_bucher || include_bayes) 
  
  scen_names <- names(scenarios)
  if (is.null(scen_names)) scen_names <- paste0("scen_", seq_along(scenarios))
  if (!parallel) {
    res_list <- lapply(seq_len(R), function(rr)
      run_rep_core(rr, scenarios, scen_names, base_seed, jags_ctrl,
                   covdefs, include_bucher, include_bayes))
  } else {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    RNGkind("L'Ecuyer-CMRG")
    if (is.null(base_seed)) base_seed <- 10L
    base_seed <- as.integer(base_seed)
    parallel::clusterSetRNGStream(cl, base_seed)
    
    export_objs <- c(
      "scenarios","scen_names","base_seed","jags_ctrl","covdefs",
      "include_bucher","include_bayes",
      # 데이터 생성/유틸
      "define_covariates","blend_schema_sigma",
      "generate_covariates","calibrate_intercept","simulate_one_study",
      "ipd_to_agd","generate_dataset_AB_BC",
      # 분석
      "do_bucher_agd_agd","run_bayes_nma_gemtc",
      "true_logOR_AB", "true_logOR_DB","true_logOR_DC",
      # 실행 유닛
      "run_rep_core"
    )
    parallel::clusterExport(cl, varlist = export_objs, envir = environment())
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(dplyr); library(tibble); library(mvtnorm); library(Matrix)
        library(rjags); library(broom); library(rlang); library(tidyr); library(gemtc)
      })
      NULL
    })
    
    res_list <- parallel::parLapply(
      cl, seq_len(R),
      function(rr) run_rep_core(rr, scenarios, scen_names, base_seed, jags_ctrl,
                                covdefs, include_bucher, include_bayes)
    )
  }
  dplyr::bind_rows(res_list)  # <- results
}

### 3. Performance summary ###
summarize_performance_min <- function(results){
  results %>%
    dplyr::group_by(scen_name, method) %>%
    dplyr::summarise(
      N_sim        = dplyr::n(),
      bias     = mean(logOR - true_AC),
      rmse     = sqrt(mean((logOR - true_AC)^2)),
      cover95  = mean(LCL <= true_AC & true_AC <= UCL),
      power    = dplyr::if_else(any(true_AC != 0),
                                mean(reject[true_AC != 0]), NA_real_),
      coverage_width  = mean(UCL - LCL),
      #Observed event rates
      mean_pA = mean(pA),
      mean_pB_AB = mean(pB_AB), 
      mean_pB_BC = mean(pB_BC),  
      mean_pC = mean(pC),
      mean_pD = mean(pD),
      mean_pB_BD = mean(pB_BD),
      mean_pC_DC = mean(pC_DC, na.rm = TRUE),
      mean_pD_DC = mean(pD_DC, na.rm = TRUE),
      
      mean_overall_rate = mean(overall_rate, na.rm=TRUE),
      
      #Covariate means
      mean_AB_C1 = mean(AB_mean_C1, na.rm = TRUE),
      mean_BC_C1 = mean(BC_mean_C1, na.rm = TRUE),
      mean_BD_C1 = mean(BD_mean_C1, na.rm = TRUE),
      mean_CD_C1 = mean(CD_mean_C1, na.rm = TRUE),
      
      mean_AB_C2 = mean(AB_mean_C2, na.rm = TRUE),
      mean_BC_C2 = mean(BC_mean_C2, na.rm = TRUE),
      mean_BD_C2 = mean(BD_mean_C2, na.rm = TRUE),
      mean_CD_C2 = mean(CD_mean_C2, na.rm = TRUE),
      
      mean_AB_C3 = mean(AB_mean_C3, na.rm = TRUE),
      mean_BC_C3 = mean(BC_mean_C3, na.rm = TRUE),
      mean_BD_C3 = mean(BD_mean_C3, na.rm = TRUE),
      mean_CD_C3 = mean(CD_mean_C3, na.rm = TRUE),
      
      #Covariate coefficient estimates
      mean_est_C1_coef_AB = mean(est_C1_coef_AB, na.rm = TRUE),
      mean_est_C1_coef_BC = mean(est_C1_coef_BC, na.rm = TRUE),
      mean_est_C1_coef_BD = mean(est_C1_coef_BD, na.rm = TRUE),
      mean_est_C1_coef_CD = mean(est_C1_coef_CD, na.rm = TRUE),
      
      mean_est_C2_coef_AB = mean(est_C2_coef_AB, na.rm = TRUE),
      mean_est_C2_coef_BC = mean(est_C2_coef_BC, na.rm = TRUE),
      mean_est_C2_coef_BD = mean(est_C2_coef_BD, na.rm = TRUE),
      mean_est_C2_coef_CD = mean(est_C2_coef_CD, na.rm = TRUE),
      
      mean_est_C3_coef_AB = mean(est_C3_coef_AB, na.rm = TRUE),
      mean_est_C3_coef_BC = mean(est_C3_coef_BC, na.rm = TRUE),
      mean_est_C3_coef_BD = mean(est_C3_coef_BD, na.rm = TRUE),
      mean_est_C3_coef_CD = mean(est_C3_coef_CD, na.rm = TRUE),
      
      .groups = "drop"
    )
}