
########## Test 1 ##########

x <- 2 + 3

expect_equal(x, 5)

########## Test 1 ##########

f = file()
sink(file=f)



# Data

x = MASS::Boston[,c(6,13)] #rm=number of rooms and lstat= percent lower status
y = MASS::Boston$medv # median value

set.seed(99)  ## MCMC posterior sampling: set seed for reproducibility
nd=200        ## number of draws to keep
burn=50      

# Run two different BARTS


 ## silence upcoming output using anonymous file connection

bf1 = invisible(meBART::wbart(x, y, nskip=burn, ndpost=nd))
bf1 <- bf1[-9]

set.seed(99)
bf2 = invisible(meBART::mebart(x, y, nskip=burn, ndpost=nd))
bf2 <- bf2[-9]
sink()

expect_equal(bf1, bf2)

