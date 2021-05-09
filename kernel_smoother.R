library(changepoint)

# signal processing functions
load_signals <- function(number_of_stations=22){
  station_num <- number_of_stations
  station_path <- "./out_data/"
  stations <- c ()
  for (num in 1: station_num){
    stations[[num]] <- read.csv(paste(station_path, "a", toString(num), ".csv", sep=""))[[2]]
  }
  return(stations)
}

find_steps <- function(signal, num_sections=3, plt=FALSE){
  if (sum(signal) == 0){
    return(signal[1:1425])
  }
  new_signal <- c()
  num <- 8
  num <- num - 1
  past_val <- 0
  for (node in 1:(length(signal) - num)){
    sequence <- seq(from=node, to=node + num)
    val <- (2*max(signal[sequence]) + mean(signal[sequence]))/3
    new_signal <- append(new_signal, val)
  }
  # plot(signal)
  # plot(new_signal)
  len <- 9
  no_zeroes <- numeric(len)
  for (i in (len + 1):(length(signal) - len)){
    sequence <- seq(from = i - len, to = i + len)
    if (signal[[i]] < max(signal) %/% 8){
      if (mean(signal[sequence]) > 1){
        no_zeroes[[i]] <- mean(signal[sequence])
      } else {
        no_zeroes[[i]] <- 0
      }
    } else {
      no_zeroes[[i]] <- signal[[i]]
    }
  }
  # plot(no_zeroes)
  new_new_signal <- c()
  num <- 8
  num <- num - 1
  past_val <- 0
  for (node in 1:(length(no_zeroes) - num)){
    sequence <- seq(from=node, to=node + num)
    val <- (2*max(no_zeroes[sequence]) + mean(no_zeroes[sequence]))/3
    new_new_signal <- append(new_new_signal, val)
  }
  smoothed_signal <- new_new_signal
  # plot(smoothed_signal)
  # plot(smoothed_signal)
  cpts <- changepoint::cpt.meanvar(new_new_signal, method="PELT", class=FALSE, minseglen = 30)
  # abline(v=cpts)
  cpts <- append(1, cpts)
  means <- c(0)
  for (i in 2:length(cpts)){
    means[[i]] <- mean(signal[seq(from=cpts[[i-1]] + 1, to=cpts[[i]])])
  }
  max_mean <- max(means)
  means <- ((num_sections - 0.0001) * means) %/% max(means)
  final_points <- c(cpts[[1]])
  for (i in 2:length(means)){
    if (means[[i]] == means[[i - 1]]){
      final_points[[length(final_points)]] <- cpts[[i]]
    } else{
      final_points <- append(final_points, cpts[[i]])
    }
  }
  final_signal <- c()
  for (i in 2:length(cpts)){
    for (j in cpts[[i-1]]:cpts[[i]]){
      final_signal <- append(final_signal, means[[i]])
    }
  }
  if (plt){
    plot(signal)
    plot(smoothed_signal)
    abline(v=final_points)
    print(final_points)
    plot(final_signal)
  }
  return(final_signal)
}


# profile model functions
create_profile <- function(step_signals, steps, past_stations=3, past_minutes=5, time_step=2){
  profile <- c()
  number_of_vars <- sum((past_minutes + 1 - past_stations):(past_minutes + 1))
  for (i in 1:number_of_vars){
    profile[[i]] <- c(0)
  }
  for (station in (1 + past_stations):length(step_signals)){
    for (min in (10 + time_step*past_minutes) : length(step_signals[[1]])){
      temp <- min_to_vec(station, min, step_signals, past_stations, past_minutes, time_step)
      for (i in 1:length(temp)){
        profile[[i]] <- append(profile[[i]], temp[[i]])
      }
    }
  }
  return(table(profile))
}

min_to_vec <- function(station, min, step_signals, past_stations=3, past_minutes=5, time_step=2){
  ans <- c()
  for (j in 0:past_stations){
    for (i in j:past_minutes){
      ans <- append(ans, step_signals[[station - j]][[min - (i * time_step)]])
    }
  }
  return(ans)
}

predict_state <- function(station_num, min, signals, profile, section_num=3, past_stations=3, past_mins=5, time_step=2){
  temp <- c()
  vals <- c()
  for (s in 0:past_stations){
    temp_station  <- signals[[station_num - s]]
    for (m in s:past_mins){
      if (m == 0){
        next
      } else {
        vals <- append(vals, temp_station[[min - (m * time_step)]] + 1)
      }
    }
  }
  for (i in 1:section_num){
    temp <- append(temp, profile[i, vals[[1]], vals[[2]], vals[[3]], vals[[4]], vals[[5]], vals[[6]], vals[[7]], vals[[8]]])
    #temp <- append(temp, profile[i, vals[[1]], vals[[2]], vals[[3]], vals[[4]], vals[[5]], vals[[6]], vals[[7]], vals[[8]], vals[[9]], vals[[10]], vals[[11]], vals[[12]], vals[[13]], vals[[14]], vals[[15]], vals[[16]], vals[[17]]])
  }
  return(which.max(temp) - 1)
}

calc_accuracy <- function(section_num=3, past_stations=3, past_mins=5, time_step=2, length_of_stations=22, num_of_tests=6){
  stations <- load_signals(length_of_stations)
  step_signals <- c()
  for (i in 1:length(stations)){
    step_signals[[i]] <- find_steps(stations[[i]], section_num)
  }
  profile <- create_profile(step_signals[1:(length_of_stations - num_of_tests)], steps=section_num, past_stations = past_stations, past_minutes =  past_mins, time_step =  time_step)
  print("profile")
  error <- 0
  total <- 0
#  for (i in (1 + past_stations):length(step_signals)){
  for (i in (length_of_stations + 1 - num_of_tests):length_of_stations){
    curr_signal <- step_signals[[i]]
    for (min in 30:(length(curr_signal) - 10)){
      total <- total + 1
      predicted <- predict_state(i, min, step_signals, profile, section_num, past_stations, past_mins, time_step) 
      if (predicted != curr_signal[[min]]){
        error <- error + 1
      }
      
    }
  }
  return(1.0 - (error/total))
}

# linear model functions

round_prediction <- function(prediction, section_num=3){
  if (prediction <= 0){
    return(0)
  }
  if (prediction >= section_num - 1){
    return(section_num - 1)
  }
  return(round(prediction))
}

create_data <- function(signals, past_stations=3, past_mins=5, time_step=2, num_of_tests, is_test_data=FALSE){
  ans <- list()
  n <- length(signals)
  cntr <- 1
  if (is_test_data){
    start_station <- n - num_of_tests + 1
    end_station <- n
  }
  else {
    start_station <- 1 + past_stations
    end_station <- n - num_of_tests
  }
  for (station in start_station : end_station){
    for (min in (10 + time_step*past_mins) : length(signals[[1]])){
      temp <- min_to_vec(station, min, signals, past_stations, past_mins, time_step)
      ans[[cntr]] <- temp
      cntr <- cntr + 1
    }
  }
  return(as.data.frame(do.call(rbind, ans)))
}

test_model <- function(testing_data, model, section_num){
  total <- 0
  error <- 0
  test_num <- nrow(testing_data)
  for (i in 1:test_num){
    total <- total + 1
    df <- testing_data[i, ]
    predicted_state <- predict(model, df[1:(length(testing_data) - 1)])
    predicted_state <- round_prediction(predicted_state, section_num = section_num)
    if (predicted_state != df[[length(testing_data)]]){
      error <- error + 1
    }
  }
  return(paste("Accuracy:", toString(1 - (error/total))))
}

create_formula <- function(past_stations=3, past_mins=5, time_step=2){
  ans <- ""
  number_of_stamps <- as.integer((past_stations + 1) * ((past_mins + 1) - (past_stations / 2)))
  ans <- paste(ans, "V", toString(number_of_stamps), " ~ ", sep="")
  for (i in 1:(number_of_stamps-2)){
    ans <- paste(ans, "V", toString(i), " + ", sep="")
  }
  ans <- paste(ans, "V", toString(number_of_stamps - 1), sep="")
  return(as.formula(ans))
}

linear_model <- function(section_num=3, past_stations=3, past_mins=5, time_step=2, length_of_stations=22, num_of_tests=6){
  stations <- load_signals(length_of_stations)
  step_signals <- c()
  for (i in 1:length(stations)){
    step_signals[[i]] <- find_steps(stations[[i]], section_num)
  }
  model_formula <- create_formula(past_stations = past_stations, past_mins = past_mins, time_step = time_step) 
  training_data <- create_data(step_signals, past_stations, past_mins, time_step, num_of_tests)
  testing_data <- create_data(step_signals, past_stations, past_mins, time_step, num_of_tests, TRUE)
  
  # training model
  model <- lm(model_formula, data = training_data)
  
  # testing model
  print(test_model(testing_data, model, section_num))
  
  return(model)
}