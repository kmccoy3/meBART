#library(treeshap)

#setwd("/Users/Wootz/Desktop/Quality Aussurance Papers")
######################################################################################################
# Get shapely value functions




splitTrees <- function(s, num_trees) {
    tree_vector = c()
    
    # Use the appropriate regex pattern to find positions of tree starts
    matches = unlist(gregexpr('\n[0-9]{1,3}\n', s))
    
    if (num_trees > 1) {
        for (i in 1:(num_trees - 1)) {
            tree_vector[i] = substring(s, matches[i], matches[i + 1] - 1)
        }
        tree_vector[num_trees] = substring(s, matches[num_trees])
    } else {
        # FIXME: does not work as expected if num_trees does not account for number of posterior draws
        tree_vector = substring(s, matches)
    }
    
    return(tree_vector)
}



getBARTTree <- function(bartModel,
                        k = 1,
                        labelVar = FALSE,
                        ntree) {
    # Extract the tree string
    
    # Vector of all the trees
    tree_vector = splitTrees(bartModel$treedraws$trees, ntree)
    
    # gets chosen tree k
    tree_str <- tree_vector[k]
    #print(tree_str)
    cut_vals = bartModel$treedraws$cutpoints
    #print(paste("-----tree-----",  tree_str))
    # Split the string into lines
    lines <- unlist(strsplit(tree_str, "\n"))
    lines <- lines[-1]
    
    #print(paste("-----Lines2-----",  lines ))
    # Parse tree nodes
    nodes = list()
    for (line in lines[-1]) {
        # Exclude the first line
        parts <- as.numeric(unlist(strsplit(line, " ")))
        node_id <- parts[1]
        nodes[[as.character(node_id)]] = parts
    }
    
    all_node_ids <- as.numeric(names(nodes))
    max_node_id <- max(all_node_ids)
    
    # Initialize empty vectors to store the tree information
    left_daughter <- numeric(max_node_id)
    right_daughter <- numeric(max_node_id)
    split_var <- numeric(max_node_id)
    split_point <- numeric(max_node_id)
    prediction <- numeric(max_node_id)
    
    # Iterate over all possible nodes
    for (node_id in 1:max_node_id) {
        if (is.null(nodes[[as.character(node_id)]])) {
            # if the node doesn't exist, continue with zeros (already initialized)
            next
        }
        
        parts <- nodes[[as.character(node_id)]]
        #print(parts[2])
        
        #print(parts)
        
        # If it's a split node
        #Need to add 1 to account for C++ base 0 and R is base 1
        
        #print(paste("-----PART3-----",  parts[3] ))
        
        if (parts[3] != 0) {
            split_var[node_id] <- parts[2] + 1
            split_point[node_id] <- cut_vals[[parts[2] + 1]][parts[3] + 1] # Convert cutpoint index to actual value
            left_daughter[node_id] <- node_id * 2
            right_daughter[node_id] <- node_id * 2 + 1
        } else {
            # If it's a leaf node
            prediction[node_id] <- parts[4]
        }
    }
    
    # Convert vectors into a data frame
    tree_df <- data.frame(
        "left daughter" = left_daughter,
        "right daughter" = right_daughter,
        "split var" = split_var,
        "split point" = split_point,
        "prediction" = prediction
    )
    
    # Handle the labelVar option
    if (labelVar && any(tree_df$`split var` != 0)) {
        tree_df$`split var` <- factor(colnames(bartModel$train.data$x)[tree_df$`split var`])
    }
    
    return(tree_df)
}




renameSplitVar <- function(tree_df, dataset) {
    # Extract column names of the dataset
    col_names <- colnames(dataset)
    # Update the split.var column
    tree_df$split.var <- ifelse(tree_df$split.var == 0, "NA", col_names[tree_df$split.var])
    
    return(tree_df)
}


###########################################
# see if we can 'unify the BART model'
library(data.table)

BART.unify <- function(rf_model,
                       data,
                       numtree,
                       ndpost,
                       ndpost_draw) {
    #if(!inherits(rf_model,'randomForest')){stop('Object rf_model was not of class "randomForest"')}
    #if(any(attr(rf_model$terms, "dataClasses") != "numeric")) {
    # stop('Models built on data with categorical features are not supported - please encode them before training.')
    #}
    
    n <- numtree
    ret <- data.table()
    
    start_tree <- (ndpost_draw - 1) * n + 1
    end_tree <- ndpost_draw * n
    
    x <- lapply(start_tree:end_tree, function(tree) {
        tree_data <- getBARTTree(rf_model, k = tree, ntree = numtree * ndpost)
    })
    
    times_vec <- sapply(x, nrow)
    
    y <- rbindlist(x)
    
    ## Fix feature names on the split.var
    # Extract column names of the dataset
    col_names <- colnames(data)
    
    # Update the split.var column
    y$split.var <- ifelse(y$split.var == 0, NA, col_names[y$split.var])
    
    
    #print(y)
    y[, Tree := rep(0:(n - 1), times = times_vec)]
    y[, Node := unlist(lapply(times_vec, function(x)
        0:(x - 1)))]
    setnames(y,
             c(
                 "Yes",
                 "No",
                 "Feature",
                 "Split",
                 "Prediction",
                 "Tree",
                 "Node"
             ))
    y[, Feature := as.character(Feature)]
    y[, Yes := Yes - 1]
    y[, No := No - 1]
    y[y$Yes < 0, "Yes"] <- NA
    y[y$No < 0, "No"] <- NA
    y[, Missing := NA]
    y[, Missing := as.integer(Missing)] # seems not, but needed
    
    ID <- paste0(y$Node, "-", y$Tree)
    y$Yes <- match(paste0(y$Yes, "-", y$Tree), ID)
    y$No <- match(paste0(y$No, "-", y$Tree), ID)
    
    y$Cover <- 0
    
    y$Decision.type <- factor(x = rep("<", times = nrow(y)),
                              levels = c("<=", "<"))
    y[is.na(Feature), Decision.type := NA]
    
    # Here we lose "Quality" information
    y[!is.na(Feature), Prediction := NA]
    
    # treeSHAP assumes, that [prediction = sum of predictions of the trees]
    # in random forest [prediction = mean of predictions of the trees]
    # so here we correct it by adjusting leaf prediction values
    #y[is.na(Feature), Prediction := Prediction / n]
    
    
    setcolorder(
        y,
        c(
            "Tree",
            "Node",
            "Feature",
            "Decision.type",
            "Split",
            "Yes",
            "No",
            "Missing",
            "Prediction",
            "Cover"
        )
    )
    
    ret <- list(model = as.data.frame(y), data = as.data.frame(data))
    class(ret) <- "model_unified"
    attr(ret, "missing_support") <- FALSE
    attr(ret, "model") <- "randomForest"
    #return( ret)
    
    return(set_reference_dataset(ret, as.data.frame(data)))
}




getBARTprediction <- function(tree_df, observation) {
    # Starting at the root node
    current_node <- 1
    
    # Recursively traverse the tree until reaching a terminal node
    while (tree_df$left.daughter[current_node] != 0 |
           tree_df$right.daughter[current_node] != 0) {
        # If the observation's value for the current split var is less than the split point
        if (observation[tree_df$split.var[current_node]] < tree_df$split.point[current_node]) {
            # Go to the left child
            current_node <- tree_df$left.daughter[current_node]
        } else {
            # Otherwise, go to the right child
            current_node <- tree_df$right.daughter[current_node]
        }
    }
    
    # Return the prediction for the terminal node
    return(tree_df$prediction[current_node])
}


#######################################################################################
library(BART)

? wbart()



##simulate data (example from Friedman MARS paper)
f = function(x) {
    10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - .5)^2 + 10 * x[, 4] + 5 *
        x[, 5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x = matrix(runif(n * 10), n, 10) #10 variables, only first 5 matter
Ey = f(x)
y = Ey + sigma * rnorm(n)
lmFit = lm(y ~ ., data.frame(x, y)) #compare lm fit to BART later

##test BART with token run to ensure installation works
set.seed(99)
bartFit = wbart(x, y, nskip = 5, ndpost = 5)

## Not run:
##run BART

ntree = 10
set.seed(99)
bartFit = wbart(x, y, ntree = ntree, ndpost = 100)

##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y, Ey, lmFit$fitted, bartFit$yhat.train.mean)
colnames(fitmat) = c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))


# Pulls the tree structure as character/string
t1 = splitTrees(bartFit$treedraws$trees, ntree)
t1

# Formats into dataframe explaining left or right split,
# which splitting varaible, and the split value
tree1 = getBARTTree(bartFit, k = 1, labelVar = FALSE, ntree)
tree1


### My additions

bartFit_1tree = wbart(x, y, ntree = 1, ndpost = 100)


