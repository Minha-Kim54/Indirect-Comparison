
#methods

#### Bucher ####
do_bucher_agd_agd <- function(agd_AB, agd_BC){
  #AB(single study)
  ab <- agd_AB %>% mutate(trt=factor(trtn, levels=c("B","A")))
  ab_A <- ab %>% filter(trt == "A")
  ab_B <- ab %>% filter(trt == "B")
  p_A <- ab_A$r / ab_A$N
  p_B <- ab_B$r / ab_B$N
  logOR_AB <- log(p_A/(1-p_A)) - log(p_B/(1-p_B))
  a <- ab_A$r; b <- ab_A$N - a; c <- ab_B$r; d <- ab_B$N - c
  var_AB <- 1/a + 1/b + 1/c + 1/d # standard error for the log odds ratio
  
  #BC(K studies)
  bc <- agd_BC %>% mutate(trt=factor(trtn, levels=c("B","C")))
  
  bc_studies <- unique(bc$studyn)
  K <- length(bc_studies)
  
  logOR_CB <- numeric(K)
  var_CB <- numeric(K)
  se_CB <- numeric(K)
  
  for(k in 1:K) {
    bc_k <- bc %>% filter(studyn == bc_studies[k])
    
    # Compute log OR from observed event rates
    bc_C <- bc_k %>% filter(trt == "C")
    bc_B <- bc_k %>% filter(trt == "B")
    
    p_C <- bc_C$r / bc_C$N
    p_B <- bc_B$r / bc_B$N
    
    logOR_CB[k] <- log(p_C/(1-p_C)) - log(p_B/(1-p_B))
    
    # Compute variance from event counts
    a <- bc_C$r
    b <- bc_C$N - a
    c <- bc_B$r
    d <- bc_B$N - c
    
    var_CB[k] <- 1/a + 1/b + 1/c + 1/d
  }
  
  #pooling by inverse-variance weighting
  w_CB <- 1 / var_CB #weight
  logOR_CB_pooled <- sum(w_CB * logOR_CB) / sum(w_CB) 
  var_CB_pooled <- 1 / sum(w_CB)
  
  # Bucher results: AC
  logOR_AC <- logOR_AB - logOR_CB_pooled
  var_AC <- var_AB + var_CB_pooled
  se_AC <- sqrt(var_AC)
  
  tibble::tibble(method="Bucher (AgD/AgD)", contrast="A vs C",
                 logOR=logOR_AC, SE=se_AC,
                 LCL=logOR_AC-1.96*se_AC, UCL=logOR_AC+1.96*se_AC)
}


#### Bayesian NMA via gemtc (fixed effects) ####

run_bayes_nma_gemtc <- function(agd_AB, agd_BC,agd_BD,agd_CD,
                                n.chains=2, n.adapt=2000, n.iter=5000, thin=6,
                                em_scale=0,seed=1000,
                                beta_0_AB = NA, beta_0_BC = NA, beta_0_BD = NA,beta_0_CD = NA){
  
  set.seed(seed)
  
  agd_data<- dplyr::bind_rows(agd_AB,agd_BC,agd_BD,agd_CD)%>%
    dplyr::select(studyn,trtn,r,N)%>%
    dplyr::rename(study = studyn,
                  treatment = trtn,
                  responders = r,
                  sampleSize = N)%>%
    dplyr::mutate(treatment=as.character(treatment))
  
  network<-mtc.network(data.ab=agd_data)
  
  model<-mtc.model(network, type="consistency",n.chain=n.chains, linearModel="fixed",re.prior.sd=1000)
  
  #MCMC
  result<-mtc.run(model, n.adapt=n.adapt, n.iter=n.iter, thin=thin)
  
  # Extract posterior samples for d_CA (= log OR_AC)
  effect_AC <- relative.effect(result, t1 = "C")
  samples <- as.matrix(as.mcmc.list(effect_AC))
  dCA_samples <- samples[, "d.C.A"]
  
  result_summary <- tibble::tibble(
    method = "Bayes NMA",
    contrast = "A vs C",
    logOR = mean(dCA_samples),
    SD = sd(dCA_samples),
    LCL = quantile(dCA_samples, 0.025, names = FALSE),
    UCL = quantile(dCA_samples, 0.975, names = FALSE),
    beta_0_AB = beta_0_AB, 
    beta_0_BC = beta_0_BC, 
    beta_0_BD = beta_0_BD,
    beta_0_CD = beta_0_CD
  )
  return(result_summary)
}