library(stats)

find_t <- function(data, k){
  N <- length(data)
  dist <- N %/% k
  start_times <- numeric(k)
  previous_times <- numeric(k + 1)
  for (i in 1:k){
    start_times[[i]] <- (i-1) * dist + 1
  }  
  start_times[[k + 1]] <- N
  while (get_dist(start_times, previous_times) > 1){
    previous_times <- start_times
    start_times <- update_times(data, previous_times)
  }
  return(start_times)
}

get_dist <- function(vec1, vec2){
  return(dist(rbind(vec1, vec2))[[1]])
}

update_times <- function(data, times){
  k <- length(times) - 1
  N <- length(data)
  i <- 1
  
  means <- numeric(k)
  sigmas <- numeric(k)
  sigma <- 0
  
  new_times <- numeric(k + 1)
  new_times[[1]] <- 1
  new_times[[k + 1]] <- N
  
  data_fragment <- data[times[[1]]:times[[2]]]
  means[[1]] <- get_mean(data_fragment, 1)
  sigmas[[1]] <- sd(data_fragment) 
  sigma <- sigma + (length(data_fragment) - 1) * sigmas[[1]]
  i <- i + length(data_fragment)
  
  for (j in 2:k){
    data_fragment <- data[times[[j]]:(times[[j + 1]]-1)]
    means[[j]] <- get_mean(data_fragment, i)
    sigmas[[j]] <- get_sd(data_fragment, i) 
    sigma <- sigma + (length(data_fragment) - 1) * sigmas[[j]]
    new_times[[j]] <- min(c((means[[j-1]] + means[[j]]) %/% 2, N))
    i <- i + length(data_fragment)
  }
  sigma <- sigma / (N - k)
  return(new_times)
}

get_mean <- function(data, starting_point){
  N <- length(data)
  if (N < 1) {
    return(starting_point)
  }
  if (sum(data) < 1){
    return(starting_point + N %/% 2)
  }
  data <- data / sum(data)
  mean_data <- 0
  for (i in starting_point:(starting_point + N - 1)){
    mean_data <- mean_data + data[[i - starting_point + 1]] * i
  }
  return(mean_data)
}

get_sd <- function(data, strt){
  N <- length(data)
  if (N < 1){
    return(0)
  }
  if (sum(data) < 1){
    return(0)
  }
  mean_data <- get_mean(data, starting_point = strt)
  data <- data / sum(data)
  sd_data <- 0
  for (i in strt:(strt + N - 1)){
    sd_data <- sd_data + data[[i - strt + 1]] * i * i
  }
  return(sd_data - mean_data * mean_data)
}

signal <- readRDS("./Farhangsara.Rds")[1:1441]
plot(signal, type = 'l', col = 'blue')
k <- 24
ans <- find_t(signal, k)
print((ans-1)/60)
abline(v = ans, col = 'red')