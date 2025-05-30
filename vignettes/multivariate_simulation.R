


library(meBART)
library(tidyverse)
library(latex2exp)
library(ggridges)

setwd("~/Documents/documents-main/research/rice_research/projects/BART/working_repo/meBART/vignettes")


p = 5

params <- list(
    notes = "Test Multivariate Run.",
        
        
    # Data generation parameters
    n = 200, # number of TOTAL observations 
    p = 5,
    sigma_meas = 0.1, # sd of measurement error
    sigma_y = 0.1, # sd of y noise
    mu_x = matrix(0.0, nrow=1, ncol=p),
    sigma_x = 0.3^2 * diag(x=1, nrow=p, ncol=p),
    
    # Tree / Model parameters
    nskip = 100, # (Default: 100L) 
    ndpost = 1000, # (Default: 1000L)
    ntree = 200, # (Default: 200L)
    meas_error_sd = 0.01*diag(x=1, nrow=p, ncol=p),
    function_type = "Friedman",
    num_iters = 25
)










get_function <- function(type=params$function_type){
    switch(
        type,
        # "indicator" = function(x){1 * (x >= 0.0)},
        # "sin" = function(x){sin(2*pi*x)},
        # "sigmoid" = function(x){1/(1 + exp(-6*x))},
        "Friedman" = function(x){sin(pi * x[,1] * x[,2]) + 2*(x[,3]-0.5)^2 + x[,4] + 0.5*x[,5]},
    )
}




generate_data <- function(MY_SEED){
    
    set.seed(MY_SEED)

    x_true <- MASS::mvrnorm(params$n, params$mu_x, params$sigma_x)
    x_meas <- x_true + MASS::mvrnorm(params$n, rep(0, p), params$meas_error_sd)
    
    my_func <- get_function(params$function_type)
    y_true <- my_func(x_true)
    y_obs <- y_true + rnorm(params$n, 0, params$sigma_y)
    
    # Convert to df
    df <- data.frame(X=x_meas, y=y_obs, X_true=x_true, y_true=y_true)
    return(df)
}


ci_calculation <- function(mdl, model){

    values <- mdl$treedraws$cutpoints[[1]]
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
}

curr_datetime <- get_file_ext()


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
    
    filename <- paste0("./results/", curr_datetime, "_bakeoff_indicator_SEED-", MY_SEED, ".png")
    
    ggsave(filename, ci_plot, width = 6, height = 4, dpi=600)
    
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
    
    filename <- paste0("./results/", curr_datetime, "_bakeoff_sine_SEED-", MY_SEED, ".png")
    
    ggsave(filename, ci_plot, width = 6, height = 4, dpi=600)
    
}





plot_x_draws <- function(mdl, df_me, MY_SEED){
    
    my_list <- list()
    
    for (k in 1:params$p){
        # k <- 1
        
        
        x_draws <- mdl$x_draws[k,,]
        x_draws_means <- rowMeans(x_draws)
        
        new_df <- data.frame(x_true_est=x_draws_means, X = df_me[1:100, k], X_true=df_me[1:100, k+6])
        
        
        # 
        # df <- data.frame(x_meas=x_meas, x_true = x_true, x_true_est = x_draws_means)
        df <- new_df[order(new_df$X_true), ]
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
            labs(title=TeX(paste0("Estimation of True Value of $X_", k, "$")), 
                 x="X",
                 color="Variable") +
            geom_segment(data=df, 
                         aes(x=X, xend=x_true_est, y=y, yend=y), 
                         arrow=arrow(length=unit(0.2,"cm"))) # +
        # scale_x_continuous(limits = c(-1, 1), breaks = seq(-1.2, 1.2, 0.2)) +
        # scale_color_discrete(labels = unname(TeX(c("Measured Value, $X_i^*$", "True Value, $X_i$", "Estimated True Value, $\\hat{X_i}$"))))
        my_list[[k]] <- plot
        
    }
    
    
    grid.arrange(grobs=my_list, ncol=3, nrow=2)
    
    # filename <- paste0("./results/", curr_datetime, "_bakeoff_x_draws_SEED-", MY_SEED, ".png")
    
    # ggsave(filename, plot, height=6, dpi=300, units="in")
    
}




results_df <- data.frame(matrix(ncol=4, nrow=0))
colnames(results_df) <- c("SEED", "Method", "Metric", "Value")

f = file()
sink(f)


for (MY_SEED in 1:params$num_iters){

    df_me <- generate_data(MY_SEED)

    df_train <- df_me[1:(params$n/2), 1:6]
    df_test <- df_me[((params$n/2) + 1):params$n, 1:6]

    

    
    # Vanilla BART
    set.seed(MY_SEED)
    mdl_BART <- BART::wbart(df_train[, -6], df_train$y, df_test[, -6],
                             nskip=params$nskip, ndpost=params$ndpost, ntree=params$ntree)
    
    
    # meBART
    set.seed(MY_SEED)
    mdl_meBART <- meBART::mebart(df_train[, -6], df_train$y,  df_test[, -6],
                             nskip=params$nskip, ndpost=params$ndpost, ntree=params$ntree, 
                             meas_error_sigma=params$meas_error_sd,                         
                             x_mu = params$mu_x,
                             x_sigma = params$sigma_x)
    

    
    # ci_BART <- ci_calculation(mdl_BART, "BART")
    # ci_meBART <- ci_calculation(mdl_meBART, "meBART")
    # 
    # BART_coverage <- get_coverage(ci_BART)
    # meBART_coverage <- get_coverage(ci_meBART)
    

    
    mse_BART <- mean((df_test$y - mdl_BART$yhat.test.mean)^2)
    mse_meBART <- mean((df_test$y - mdl_meBART$yhat.test.mean)^2)
    
    # my_func <- get_function(params$function_type)
    
    # mse_func_BART <- mean((ci_BART$y - my_func(ci_BART$X))^2)
    # mse_func_meBART <- mean((ci_meBART$y - my_func(ci_meBART$X))^2)
    
    mae_BART <- mean(abs(df_test$y - mdl_BART$yhat.test.mean))
    mae_meBART <- mean(abs(df_test$y - mdl_meBART$yhat.test.mean))
    
    ar_meBART <- mean(mdl_meBART$acceptances)
    
    # crps.BART <- get_crps(mdl_BART)
    # crps.meBART <- get_crps(mdl_meBART)
    
    # plot_x_draws(mdl_meBART, df_train, MY_SEED)
    
    # a <- params$nskip + 2
    # b <- params$nskip + params$ndpost + 1
    # new_error <- sd(df_train$X_true - rowMeans(mdl_meBART$x_draws[,,a:b]))
    
    # results_df <- rbind(results_df, data.frame("SEED"=MY_SEED,
    #                                       "Method"= "meBART",
    #                                       "Metric"= "New Error SD",
    #                                       "Value" = new_error))

    
    results_df <- rbind(results_df, data.frame("SEED"=MY_SEED,
                                          "Method"= "meBART",
                                          "Metric"= "Acceptance Ratio",
                                          "Value" = ar_meBART))
    
    
    
    
    
    # results_df <- rbind(results_df, data.frame("SEED"=rep(MY_SEED, 10),
    #                                       "Method"= rep(c("Vanilla BART", "meBART"), each=5),
    #                                       "Metric"= rep(c("Coverage", "MSE", "MAE", "MSE_func", "CRPS"), times=2),
    #                                     "Value" = c(BART_coverage, mse_BART, mae_BART, mse_func_BART, crps.BART, 
    #                                                 meBART_coverage, mse_meBART, mae_meBART, mse_func_meBART, crps.meBART)))
    
    results_df <- rbind(results_df, data.frame("SEED"=rep(MY_SEED, 4),
                                               "Method"= rep(c("Vanilla BART", "meBART"), each=2),
                                               "Metric"= rep(c("MSE", "MAE"), times=2),
                                               "Value" = c(mse_BART, mae_BART, mse_meBART, mae_meBART)))
    
    
    # 
    # if (params$function_type == "indicator"){
    #     plot_indicator(ci_df_combined = rbind(ci_BART, ci_meBART), 
    #             df_me = df_train, 
    #             MY_SEED = MY_SEED)
    # } else if (params$function_type == "sin"){
    #     plot_sine(ci_df_combined = rbind(ci_BART, ci_meBART), 
    #             df_me = df_train, 
    #             MY_SEED = MY_SEED)
    # } else if (params$function_type == "sigmoid"){
    #     plot_sine(ci_df_combined = rbind(ci_BART, ci_meBART), 
    #               df_me = df_train, 
    #               MY_SEED = MY_SEED)  
    # }
    
    



}

closeAllConnections()





filename <- paste0("./results/", curr_datetime, "_full_results_indicator.csv")
write_csv(results_df, filename)


# Get metric names
metric_names <- c("Coverage", "MSE", "MAE", "MSE_func", "CRPS")
metric_names <- c("MSE", "MAE")

# Step 1: Summarize with lb, mean, ub
results_summary <- results_df %>% 
  filter(Metric %in% metric_names) %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  group_by(Method) %>%
  summarise(across(everything(), list(
    lb = ~ quantile(., 0.025),
    mean = mean,
    ub = ~ quantile(., 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  ungroup()

# Step 2: Combine mean, lb, ub into single string columns
format_estimate <- function(mean, lb, ub) {
  sprintf("%.3f (%.3f, %.3f)", mean, lb, ub)
}


# Step 3: Create formatted columns
for (metric in metric_names) {
  results_summary[[metric]] <- format_estimate(
    results_summary[[paste0(metric, "_mean")]],
    results_summary[[paste0(metric, "_lb")]],
    results_summary[[paste0(metric, "_ub")]]
  )
}

# Step 4: Keep only Method and formatted metrics
results_final <- results_summary %>%
  select(Method, all_of(metric_names))


filename <- paste0("./results/", curr_datetime, "_summary_results_indicator.csv")
write_csv(results_final, filename)


filename <- paste0("./results/", curr_datetime, "_params.txt")
sink(filename)

for (param in names(params)) {
    cat(paste0(param, ": ", params[[param]], "\n"))
}

sink()



