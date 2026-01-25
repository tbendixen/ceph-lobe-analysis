# To ensure that warnings are not turned to errors when running MCMC,
# and to automatically build tools required for Stan
options(warn = 1,
        buildtools.check = function(action) TRUE )

### wrapper function to brms
brms_func <-   function(formula, constrain_imp = FALSE, vars, iter = 20000) {
  
  # Extract response variables from formula
  bterms <- brmsterms(formula)
  resp_vars <- names(bterms$terms)
  
  # Helper function to check if response has covariates
  has_coef <- function(resp_idx) {
    length(all.vars(bterms$dpars[[resp_idx]]$fe)) > 1
  }
  
  # Helper function to check if response has random effects
  has_re <-  function(resp_idx) {
    length(bterms$dpars[[resp_idx]]$re$form) > 0
  }
  
  # Build prior list dynamically
  prior_list <-   set_prior("normal(0, 0.5)", class = "Intercept", resp = resp_vars[1])
  if (has_coef(1)) prior_list <- prior_list + set_prior("normal(0, 0.25)", class = "b", resp = resp_vars[1])
  if (has_re(1)) prior_list <- prior_list + set_prior("exponential(1)", class = "sd", resp = resp_vars[1])
  prior_list <- prior_list + set_prior("exponential(2)", class = "sigma", resp = resp_vars[1])
  
  if(length(resp_vars) > 1) {
    for (i in 2:length(resp_vars)) {
      r <-   resp_vars[i]
      prior_list <- prior_list + set_prior("normal(0, 0.5)", class = "Intercept", resp = r)
      if (has_coef(i)) prior_list <- prior_list + set_prior("normal(0, 0.25)", class = "b", resp = r)
      if (has_re(i)) prior_list <- prior_list + set_prior("exponential(1)", class = "sd", resp = r)
      prior_list <- prior_list + set_prior("exponential(2)", class = "sigma", resp = r)
    }
  }
  
  if(constrain_imp) {
    model <- brm(formula,
                 prior = prior_list,
                 data = d,
                 data2 = list(cormat = cormat),
                 chains = 0,
                 seed = 1,
                 backend = "rstan",
                 cores = 4) 
    
    scode <- stancode(model)
    imp_code <- scode
    
    for (var in vars) {
      brms_var <- gsub("[\\._]", "", var)
      stan_var <- paste0("Ymi_", brms_var)
      
      lower_bound <- min(model$data[[var]], na.rm = TRUE)
      upper_bound <- max(model$data[[var]], na.rm = TRUE)
      
      imp_code <- gsub(paste0("vector[Nmi_", brms_var, "] ", stan_var, ";"),
                       paste0("vector<lower=", lower_bound, ", upper=", upper_bound, ">[Nmi_", brms_var, "] ", stan_var, ";"),
                       imp_code, fixed = TRUE)
    }
    
    stan_model_obj <-  stan_model(model_code = imp_code)
    attributes(model$fit)$stanmodel <- stan_model_obj
    
    model_constrain <- update(model,
                              cores = 4, chains = 4, iter = iter,
                              recompile = FALSE, 
                              control = list(adapt_delta = 0.99))
    
    return(model_constrain)
  }
  
  else {
    brm(formula,
        prior = prior_list,
        data = d,
        data2 = list(cormat = cormat),
        iter = iter, chains = 4,
        seed = 1,
        control = list(adapt_delta = 0.99),
        backend = "rstan",
        cores = 4) 
  }
}

### function to extract posterior median of imputed values
impute <- function(vars, model) {
  
  # Get data from brms model object
  data <- model$data
  
  # Match brms var names
  ymi_vars <- gsub("[\\._]", "", vars)
  
  for (i in seq_along(vars)) {
    # Extract imputed values from brms model
    imputed_values <- brms::as_draws_df(model, variable = paste0("^Ymi_", ymi_vars[i]), regex = TRUE)
    
    # Remove unnecessary columns
    imputed_values$.chain <- NULL
    imputed_values$.iteration <- NULL
    imputed_values$.draw <- NULL
    
    # Identify missing indices in the dataset
    miss <- which(is.na(data[[vars[i]]]))
    
    # Compute median of imputed values
    imputed_median <- apply(imputed_values[miss, ], 2, median)
    
    # Create column indicating missing values, useful for plotting
    data[[paste0(vars[i], "_imputed")]] <- ifelse(is.na(data[[vars[i]]]), 1, 0)
    
    # Impute missing values in the dataset
    data[miss, ][[vars[i]]] <- imputed_median
    
  }
  return(data)
}


### function for lineribbon post pred plot
line_colors <- c("#2C70FE","black")

lineribbon_plot <- function(e, x_var, y_var, g_var = benthic, data, interaction = TRUE) {
  
  if (interaction) {
    p <- ggplot() + 
      stat_lineribbon(data = e, aes(x = {{ x_var }}, y = .epred, 
                                    color = as.factor({{ g_var }}), fill = as.factor({{ g_var }})), 
                      .width = c(.95),
                      point_interval = "median_hdci", 
                      alpha = 0.5) + 
      scale_fill_manual(values = line_colors, guide = "none") +
      scale_color_manual(name = "Habitat",
                         values = c("0" = line_colors[1], "1" = line_colors[2]),
                         labels = c("Pelagic", "Benthic")) +
      theme(axis.text.x = element_text(size = 13)) +
      theme_classic()
    
    if (paste0(deparse(substitute(x_var)), "_imputed") %in% names(data)) {
      p <- p + geom_point(aes(x = {{ x_var }}, y = {{ y_var }}, 
                              color = as.factor({{ g_var }}),
                              shape = factor(!!sym(paste0(deparse(substitute(x_var)), "_imputed")))),
                          alpha = 0.75,
                          data = data) + 
        scale_shape_manual(values = c("0" = 16, "1" = 1), guide = "none") 
    } else {
      p <- p + geom_point(aes(x = {{ x_var }}, y = {{ y_var }},
                              color = as.factor({{ g_var }})),
                          alpha = 0.75,
                          data = data) + 
        scale_shape_manual(values = 16, guide = "none")
    }
    
  } else {
    p <- ggplot() + 
      stat_lineribbon(data = e, aes(x = {{ x_var }}, y = .epred),
                      color = "#2C70FE",
                      fill = "grey40",
                      .width = c(.95),
                      point_interval = "median_hdci", 
                      alpha = 0.5,
                      show.legend = FALSE) + 
      theme(axis.text.x = element_text(size = 13)) +
      theme_classic()
    
    if (paste0(deparse(substitute(x_var)), "_imputed") %in% names(data)) {
      p <- p + geom_point(aes(x = {{ x_var }}, y = {{ y_var }}, shape = factor(!!sym(paste0(deparse(substitute(x_var)), "_imputed")))),
                          alpha = 0.75,
                          data = data) + 
        scale_shape_manual(values = c("0" = 16, "1" = 1), guide = "none") 
    } else {
      p <- p + geom_point(aes(x = {{ x_var }}, y = {{ y_var }}),
                          alpha = 0.75,
                          data = data) + 
        scale_shape_manual(values = 16, guide = "none")
    }
    
  }
  
  return(p)
}

### function for categorical post pred plot
cat_colors <- c("#2C70FE","#feba2c")

cat_plot <- function(e, x_var, y_var, g_var, interaction = FALSE, data) {
  
  if (interaction) {
    ggplot() + 
      geom_point(aes(x = {{x_var}},
                     y = {{y_var}},
                     color = as.factor({{ g_var }}), 
                     fill = as.factor({{ g_var }})), 
                 data = data, position = position_jitter(width = 0.15, seed = 1), alpha = 0.4, show.legend = FALSE) +
      stat_pointinterval(data = e, aes(x = {{x_var}}, 
                                       y = .epred,
                                       color = as.factor({{ g_var }}), 
                                       fill = as.factor({{ g_var }})),
                         point_interval = median_hdci,
                         .width = 0.95,
                         linewidth = 2,
                         size = 5,
                         position = position_dodge(width = 0.2)) +
      theme_classic() + 
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
  
  } else {
    ggplot() + 
    geom_point(aes(x = {{x_var}}, y = {{y_var}}, color = {{x_var}}), 
              data = data, position = position_jitter(width = 0.15, seed = 1), alpha = 0.4) +
    stat_pointinterval(data = e, aes(x = {{x_var}}, y = .epred),
                      point_interval = median_hdci,
                      .width = 0.95,
                      linewidth = 2,
                      size = 5) +
    scale_color_manual(values = cat_colors) +
    scale_fill_manual(values = cat_colors) +
    theme_classic() + 
    theme(legend.position = "none")  + 
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) + 
    theme(legend.position = "none")
  }
}

### function for gradient-density post pred plot
gdensity_colors <- c("#feba2c", "#2C70FE")
cat4_colors <- c("#feba2c", "#6A5ACD", "#4CAF50", "#2C70FE")

gdensity_plot <- function(e, x_var, y_var, g_var, interaction = FALSE, data) {
  
  # Check if x_var should be reversed (only for benthic and sociality)
  x_var_name <- deparse(substitute(x_var))
  should_reverse <- x_var_name %in% c("benthic", "sociality.bin")
  
  if (interaction) {
    ggplot() + 
      stat_slab(aes(x = if(should_reverse) forcats::fct_rev(as.factor({{x_var}})) else as.factor({{x_var}}),
                    y = {{y_var}},
                    fill = as.factor({{ g_var }}),
                    slab_alpha = after_stat(f)), 
                data = data,
                width = 0.5,
                n = 2000, 
                show.legend = FALSE) +
      stat_pointinterval(data = e, aes(x = if(should_reverse) forcats::fct_rev(as.factor({{x_var}})) else as.factor({{x_var}}),
                                       y = .epred,
                                       color = as.factor({{ g_var }}), 
                                       fill = as.factor({{ g_var }})),
                         point_interval = median_hdci,
                         .width = 0.95,
                         linewidth = 2,
                         position = position_dodge(width = 0.2)) +
      theme_classic() + 
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
    
  } else {
    ggplot() + 
      stat_slab(aes(x = if(should_reverse) forcats::fct_rev(as.factor({{x_var}})) else as.factor({{x_var}}),
                    y = {{y_var}},
                    fill = as.factor({{ x_var }}),
                    slab_alpha = after_stat(f)), 
                data = data,
                width = 0.5,
                n = 2000, 
                show.legend = FALSE) +
      stat_pointinterval(data = e, aes(x = if(should_reverse) forcats::fct_rev(as.factor({{x_var}})) else as.factor({{x_var}}),
                                       y = .epred),
                         point_interval = median_hdci,
                         .width = 0.95,
                         linewidth = 2) +
      theme_classic() + 
      scale_color_manual(values = gdensity_colors) +
      scale_fill_manual(values = gdensity_colors) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) + 
      theme(legend.position = "none")
  }
}
