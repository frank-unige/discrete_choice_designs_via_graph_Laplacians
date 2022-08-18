
sim_study <- function(d,k,sig, tol = 1e-8, rng = NULL){
  ## perform a simulation study to measure performance of the EMTP2 block descent algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - method: algorithm for optimization.
  ##    - tol: tolerance for algorithm.
  ## Returns:
  ##  a tibble with time, KKT condition values and duality gap
  
  # check arguments
  
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  # perform simulation

  design_time <- numeric(1)
  dirder_max <- numeric(1)  
  eff <- numeric(1)
  support <- numeric(1)
  
  b=rnorm(n=d,mean=0,sd=sig)
  p=exp(b)
  
  ptm <- proc.time()[1]
  result <- discrete_choice_design(p=p,k=k)
  design_time <- proc.time()[1] - ptm

  dirder_max <- result$dirder_max
  dual_gap <- (d-1)*(sum(result$w) - 1)
  support <- sum(abs(result$w)>1e-8)
  
  Q_complete= diag(L(p=p))%*%hyper_incidence(d,k)%*%solve(diag(R(p=p,k=k)))%*%rep(1/choose(d,k),choose(d,k))
  Theta_complete_d = Qvec2Theta(Q_complete)[-d,-d]
  eff = (det(Theta_complete_d)*CM(Gamma2vec(result$Gamma_opt)))^(1/(d-1))
  
  tbl <- tibble(type = paste0("design_time"), 
                value = design_time) %>% 
    bind_rows(tibble(type = paste0("dirder_max"),
                     value =  dirder_max))%>% 
    bind_rows(tibble(type = paste0("dual_gap"),
                     value =  dual_gap))%>% 
    bind_rows(tibble(type = paste0("efficiency"),
                     value =  eff))%>% 
    bind_rows(tibble(type = paste0("support"),
                     value =  support))
  return(tbl)
}


wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Args:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_study, args = fun_args[i, ]) %>% 
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}



#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



