# source: https://www.youtube.com/watch?v=_5hXQPTW-wU

library(foreach)
library(doParallel)
library(ggplot2)

# set up parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores-1) # not to overload computer
registerDoParallel(cl)
n = 1000
big_list <- list()

big_list <- foreach(i = 1:n) %dopar% {
  big_list[i] <- i^2
}

stopCluster(cl)
