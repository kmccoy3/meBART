


library(meBART)
library(tidyverse)
library(latex2exp)
library(ggridges)
library(progress)


setwd("~/Documents/documents-main/research/rice_research/projects/BART/working_repo/meBART/vignettes")


params <- list(
    notes = "Regular combo.",
    short_desc = "_combo",
        
        
    # Data generation parameters
    n = 200, # number of TOTAL observations 
    p = 1,
    sigma_meas = 0.1, # sd of measurement error
    sigma_y = 0.1, # sd of y noise
    mu_x = matrix(0.0),
    sigma_x = matrix(0.3),
    
    # Tree / Model parameters
    nskip = 100, # (Default: 100L) 
    ndpost = 1000, # (Default: 1000L)
    ntree = 200, # (Default: 200L)
    sigdf = 3, # degrees of freedom for sigma prior (Default: 3L)
    sigquant = .90, # quantile for sigma prior (Default: 0.90)
    
    
    meas_error_sd = 0.01*diag(x=1, nrow=1, ncol=1),
    function_type = "combo",
    num_iters = 25,
    testing = FALSE
)


pb <- progress_bar$new(total = params$num_iters)









get_function <- function(type=params$function_type){
    switch(
        type,
        "indicator" = function(x){1 * (x >= 0.0)},
        "sin" = function(x){sin(2*pi*x)},
        "combo" = function(x){cos(pi*x)*(x<=0) - cos(pi*x)*(x>0)}
    )
}




generate_data <- function(MY_SEED){
    
    set.seed(MY_SEED)

    x_true <- rnorm(params$n, as.numeric(params$mu_x), as.numeric(params$sigma_x))
    x_meas <- x_true + rnorm(params$n, 0, params$sigma_meas)
    
    my_func <- get_function(params$function_type)
    y_true <- my_func(x_true)
    y_obs <- y_true + rnorm(params$n, 0, params$sigma_y)
    
    # Convert to df
    df <- data.frame(X=x_meas, y=y_obs, X_true=x_true, y_true=y_true)
    return(df)
}


ci_calculation <- function(mdl, model){
        
    values <- seq(from=as.numeric(params$mu_x - 2*params$sigma_x), to=as.numeric(params$mu_x + 2*params$sigma_x), length.out = 100)
    output <- predict(mdl, data.frame(values))
    means <- colMeans(output)
    
    lower_CI <- c()
    upper_CI <- c()
    
    for (i in 1:ncol(output)) {
        lower_CI <- c(lower_CI, quantile(output[,i], 0.025))
        upper_CI <- c(upper_CI, quantile(output[,i], 0.975))
    }
    
    ci_df <- data.frame(X=values, y=means, lower_CI=lower_CI, upper_CI=upper_CI)
    
    ci_df$model <- model

    
    return(ci_df)
        
}


get_crps <- function(mdl){
    
    values <- seq(from=as.numeric(params$mu_x - 2*params$sigma_x), to=as.numeric(params$mu_x + 2*params$sigma_x), length.out = 100)
    output <- predict(mdl, data.frame(values))
    n <- nrow(mdl$varcount)
    
    my_func <- get_function(params$function_type)
    
    ys <- my_func(values)
    
    crps <- loo::crps(x=output[1:(n/2),], x2=output[(n/2+1):n,], y=ys)
    
    return(-1 * crps$estimates[[1]])

}


get_file_ext <- function(){
    datetime <- str_replace_all(format(Sys.time()), c(" "="_", ":"="-"))
    return(paste0(datetime, params$short_desc))
}

curr_datetime <- get_file_ext()



dir.create(file.path(paste0("results/", curr_datetime)))







get_coverage <- function(ci_df){
    # Coverage calculation
    
    coverage <- c()
    
    my_func <- get_function(params$function_type)
    
    for (i in 1:nrow(ci_df)){
        x_val <- my_func(ci_df$X[i])
        
        # x_val <- sin(2*pi*ci_df$X[i]) # TODO: Remove hardcoding
        
        if (ci_df$lower_CI[i] <= x_val && ci_df$upper_CI[i] >= x_val){
            coverage <- c(coverage, 1)
        } else {
            coverage <- c(coverage, 0)
        }
    }
    
    return(mean(coverage))
}


plot_indicator <- function(ci_df_combined, df_me, MY_SEED){
    
    
    ci_plot <- ggplot(ci_df_combined) +
        facet_grid(cols = vars(model)) +
        geom_line(aes(x=X, y=y), linetype=3) +
        geom_ribbon(aes(x=X, ymin=lower_CI, ymax=upper_CI), alpha=0.2) +
        labs(title = "Vanilla BART vs. meBART: Mean Prediction with 95% CI", 
             x = "X", 
             y = "Y") +
        geom_point(data=df_me, aes(x = X, y = y), color="#619CFF") +
        annotate("segment", x = -0.50, xend = 0, y = 0, yend = 0, linetype = "solid") +
        annotate("segment", x = 0, xend = 0.5, y = 1, yend = 1, linetype = "solid") +
        annotate("point", x = 0, y = 0, size = 3, shape = 1) +
        annotate("point", x = 0, y = 1, size = 3, shape = 16) # +
        # scale_x_continuous(breaks = seq(-1.2, 1.2, 0.2)) +
        # scale_y_continuous(limits = c(-0.25, 1.25), breaks = seq(-0.2, 1.2, 0.2))
    
    filename <- paste0("./results/", curr_datetime, "/SEED-", MY_SEED, "_CIs.png")
    
    suppressMessages(ggsave(filename, ci_plot, width = 6, height = 4, dpi=600))
    
}



plot_sine <- function(ci_df_combined, df_me, MY_SEED){
    
    values <- seq(from=as.numeric(params$mu_x - 3*params$sigma_x), to=as.numeric(params$mu_x + 3*params$sigma_x), length.out = 100)
    my_func <- get_function(params$function_type)
    
    
    ci_plot <- ggplot(ci_df_combined) +
        facet_grid(cols = vars(model)) +
        geom_line(aes(x=X, y=y), linetype=3) +
        geom_ribbon(aes(x=X, ymin=lower_CI, ymax=upper_CI), alpha=0.2) +
        labs(title = "meBART: Mean Prediction with 95% CI", x = "X", y = "Y") +
        geom_point(data=df_me, aes(x = X, y = y), color="#619CFF") +
        # scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
        geom_line(data=data.frame(x=values, y=my_func(values)), aes(x = x, y = y), color="#F8766D")
    
    filename <- paste0("./results/", curr_datetime, "/SEED-", MY_SEED, "CIs.png")
    
    suppressMessages(ggsave(filename, ci_plot, width = 6, height = 4, dpi=600))
    
}




plot_x_draws <- function(mdl, df_me, MY_SEED){
    
    x_draws <- mdl$x_draws[1,,]
    x_draws_means <- rowMeans(x_draws)
    
    df_me <- df_me[, c("X", "X_true")]
    
    df_me$x_true_est <- x_draws_means
    
    # 
    # df <- data.frame(x_meas=x_meas, x_true = x_true, x_true_est = x_draws_means)
    df <- df_me[order(df_me$X_true), ]
    df$y <- 1:(params$n/2)
    
    df2 <- reshape2::melt(df, id="y")
    
    plot <- ggplot(df2) +
        geom_point(aes(x=value, y=y, color=variable), 
                   shape=19, 
                   size=2) + 
        theme(axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position = "bottom") +
        labs(title=TeX("Estimation of True Value of $X_i$"), 
             x="X",
             color="Variable") +
        geom_segment(data=df, 
                     aes(x=X, xend=x_true_est, y=y, yend=y), 
                     arrow=arrow(length=unit(0.2,"cm"))) +
        # scale_x_continuous(limits = c(-1, 1), breaks = seq(-1.2, 1.2, 0.2)) +
        scale_color_discrete(labels = unname(TeX(c("Measured Value, $X_i^*$", "True Value, $X_i$", "Estimated True Value, $\\hat{X_i}$"))))

    filename <- paste0("./results/", curr_datetime, "/SEED-", MY_SEED, "_xdraws.png")
    
    suppressMessages(ggsave(filename, plot, height=6, dpi=300, units="in"))
    
}

plot_sigmas <- function(mdl_BART, mdl_meBART, MY_SEED){
    
    
    sigma_df <- data.frame(iter_num = seq_along(mdl_meBART$sigma), sigma = mdl_meBART$sigma, method = "meBART")
    sigma_df2 <- data.frame(iter_num = seq_along(mdl_BART$sigma), sigma = mdl_BART$sigma, method = "Vanilla BART")
    
    sigma_df <- rbind(sigma_df, sigma_df2)
    
    sigma_df$method <- factor(sigma_df$method, levels = c("Vanilla BART", "meBART"))
    
    # Create the ggplot line plot
    sigma_plot <- ggplot(sigma_df, aes(x = iter_num, y = sigma, color=method)) +
        geom_line() +
        geom_point(alpha=0) +
        labs(x = "MCMC Iteration", y = TeX("sigma, $\\sigma$"), title = "\nPosterior Draws of Sigma") +
        geom_hline(yintercept=params$sigma_y, color="black", linetype="dashed") +
        theme(legend.position = "bottom")

    filename <- paste0("./results/", curr_datetime, "/SEED-", MY_SEED, "_sigma_plots.png")
    suppressMessages(ggsave(filename, sigma_plot, height=6, dpi=300, units="in"))
}


calculate_sigma_crps <- function(mdl_BART, mdl_meBART){
    n <- nrow(mdl_meBART$varcount)
    BART_sigma_crps <- loo::crps(x=mdl_BART$sigma[1:(n/2)], x2=mdl_BART$sigma[(n/2+1):n], y=params$sigma_y)
    
    meBART_sigma_crps <- loo::crps(x=mdl_meBART$sigma[1:(n/2)], x2=mdl_meBART$sigma[(n/2+1):n], y=params$sigma_y)

    return(c(-1*BART_sigma_crps$estimates[[1]], -1* meBART_sigma_crps$estimates[[1]]))
    
}




results_df <- data.frame(matrix(ncol=4, nrow=0))
colnames(results_df) <- c("SEED", "Method", "Metric", "Value")








for (MY_SEED in 1:params$num_iters){

    df_me <- generate_data(MY_SEED)

    df_train <- df_me[1:(params$n/2),]
    df_test <- df_me[((params$n/2) + 1):params$n,]

    
    f = file()
    sink(f)
    
    # Vanilla BART
    set.seed(MY_SEED)
    mdl_BART <- BART::wbart(df_train$X, df_train$y, df_test$X,
                            nskip=params$nskip, ndpost=params$ndpost, ntree=params$ntree, 
                            sigdf=params$sigdf, sigquant=params$sigquant)
    
    
    # meBART
    set.seed(MY_SEED)
    mdl_meBART <- meBART::mebart(df_train$X, df_train$y,  df_test$X,
                             nskip=params$nskip, ndpost=params$ndpost, ntree=params$ntree,
                             sigdf=params$sigdf, sigquant=params$sigquant,
                             meas_error_sigma=params$meas_error_sd,                         
                             x_mu = params$mu_x,
                             x_sigma = params$sigma_x)
    
    
    

    
    ci_BART <- ci_calculation(mdl_BART, "BART")
    ci_meBART <- ci_calculation(mdl_meBART, "meBART")
    
    BART_coverage <- get_coverage(ci_BART)
    meBART_coverage <- get_coverage(ci_meBART)
    

    
    mse_BART <- mean((df_test$y - mdl_BART$yhat.test.mean)^2)
    mse_meBART <- mean((df_test$y - mdl_meBART$yhat.test.mean)^2)
    
    my_func <- get_function(params$function_type)
    
    mse_func_BART <- mean((ci_BART$y - my_func(ci_BART$X))^2)
    mse_func_meBART <- mean((ci_meBART$y - my_func(ci_meBART$X))^2)
    
    mae_BART <- mean(abs(df_test$y - mdl_BART$yhat.test.mean))
    mae_meBART <- mean(abs(df_test$y - mdl_meBART$yhat.test.mean))
    
    ar_meBART <- mean(mdl_meBART$acceptances)
    
    crps.BART <- get_crps(mdl_BART)
    crps.meBART <- get_crps(mdl_meBART)
    
    if (!params$testing) plot_x_draws(mdl_meBART, df_train, MY_SEED)
    
    a <- params$nskip + 2
    b <- params$nskip + params$ndpost + 1
    new_error <- sd(df_train$X_true - rowMeans(mdl_meBART$x_draws[,,a:b]))
    
    sigma_crps <- calculate_sigma_crps(mdl_BART, mdl_meBART)
    BART_sigma_crps <- sigma_crps[1]
    meBART_sigma_crps <- sigma_crps[2]

    results_df <- rbind(results_df, data.frame("SEED"=MY_SEED,
                                          "Method"= "meBART",
                                          "Metric"= "New Error SD",
                                          "Value" = new_error))

    
    results_df <- rbind(results_df, data.frame("SEED"=MY_SEED,
                                          "Method"= "meBART",
                                          "Metric"= "Acceptance Ratio",
                                          "Value" = ar_meBART))
    
    
    
    
    
    results_df <- rbind(results_df, data.frame("SEED"=rep(MY_SEED, 12),
                                          "Method"= rep(c("Vanilla BART", "meBART"), each=6),
                                          "Metric"= rep(c("Coverage", "MSE", "MAE", "MSE_func", "CRPS", "Sigma CRPS"), times=2),
                                        "Value" = c(BART_coverage, mse_BART, mae_BART, mse_func_BART, crps.BART, BART_sigma_crps,
                                                    meBART_coverage, mse_meBART, mae_meBART, mse_func_meBART, crps.meBART, meBART_sigma_crps)))
    
    if (!params$testing) plot_sigmas(mdl_BART, mdl_meBART, MY_SEED)

    
    if (!params$testing){
        if (params$function_type == "indicator"){
            plot_indicator(ci_df_combined = rbind(ci_BART, ci_meBART), 
                    df_me = df_train, 
                    MY_SEED = MY_SEED)
        } else if (params$function_type == "sin"){
            plot_sine(ci_df_combined = rbind(ci_BART, ci_meBART), 
                    df_me = df_train, 
                    MY_SEED = MY_SEED)
        } else if (params$function_type == "combo"){
            plot_sine(ci_df_combined = rbind(ci_BART, ci_meBART), 
                      df_me = df_train, 
                      MY_SEED = MY_SEED)  
        }
    }
    
    
    sink()

    pb$tick()
}







filename <- paste0("./results/", curr_datetime, "/full_results_indicator.csv")
if (!params$testing) write_csv(results_df, filename)


# Get metric names
metric_names <- c("Coverage", "MSE", "MAE", "MSE_func", "CRPS", "Sigma CRPS")

# Step 1: Summarize with lb, median, ub
results_summary <- results_df %>% 
  filter(Metric %in% metric_names) %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  group_by(Method) %>%
  summarise(across(everything(), list(
    lb = ~ quantile(., 0.025),
    median = median,
    ub = ~ quantile(., 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  ungroup()

# Step 2: Combine median, lb, ub into single string columns
format_estimate <- function(median, lb, ub) {
  sprintf("%.3f (%.3f, %.3f)", median, lb, ub)
}


# Step 3: Create formatted columns
for (metric in metric_names) {
  results_summary[[metric]] <- format_estimate(
    results_summary[[paste0(metric, "_median")]],
    results_summary[[paste0(metric, "_lb")]],
    results_summary[[paste0(metric, "_ub")]]
  )
}

# Step 4: Keep only Method and formatted metrics
results_final <- results_summary %>%
  dplyr::select(Method, all_of(metric_names))


filename <- paste0("./results/", curr_datetime, "/summary_results_indicator.csv")

if (!params$testing){
    write_csv(results_final, filename)


    filename <- paste0("./results/", curr_datetime, "/params.txt")
    sink(filename)
    
    for (param in names(params)) {
        cat(paste0(param, ": ", params[[param]], "\n"))
    }
    
    sink()
}



