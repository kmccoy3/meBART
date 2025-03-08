


library(tidyverse)

# Build BART model
ntree = 100
ndpost = 100
bartFit = wbart(x, y, ntree = ntree, ndpost = ndpost)



var_usage <- function(bartModel, depth=0){

    # Extract the tree structure
    tc <- textConnection(bartModel$treedraws$tree)
    full_df <- read.table(file = tc, fill = TRUE, row.names = NULL, header = FALSE, col.names = c("node_id", "variable", "cutpoint", "prediction"))
    
    # Remove first row with MCMC info
    all_trees_df <- full_df[-1,] %>% drop_na()
    
    num_tot_trees <- sum(!complete.cases(full_df)) - 1
    
    all_trees_df$depth <- floor(log2(all_trees_df$node_id))
    
    depth0_rules <- all_trees_df[all_trees_df$depth == depth, ] 
    
    var_count <- depth0_rules %>% count(variable)
    
    splitrule_counts <- depth0_rules %>% count(variable, cutpoint)
    
    splitrule_counts$variable <- as.numeric(splitrule_counts$variable)
    splitrule_counts$cutpoint <- as.numeric(splitrule_counts$cutpoint)
    
    for (i in 1:nrow(splitrule_counts)){
        splitrule_counts$cutpoint_value[i] <- bartFit$treedraws$cutpoints[[splitrule_counts$variable[i]+1]][splitrule_counts$cutpoint[i]+1]
    }


    result <- list(var_count=var_count, splitrule_counts=splitrule_counts, num_tot_trees=num_tot_trees)
    
    close(tc)
    
    return(result)
    
    
}

tmp <- var_usage(bartFit)




# Example

hist(tmp$splitrule_counts$cutpoint_value[tmp$splitrule_counts$variable == 0], breaks = 20, main = "Cutpoint values", xlab = "Cutpoint value")







# TODO: Add new / different datasets
# TODO: Add ggplot2 plots
# TODO: Apply to measurement error dataset (1D)

