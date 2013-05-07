# Libraries ####
library(deSolve)
library(reshape2)
library(ggplot2)
source("phase_diagrams.R")

# Wrapper for ode simulation of EAB models ####

sim <- function(initial_df, times=seq(from=0, to=10, by=0.1), factory=eab_factory, min_rates="auto", max_tries=1, ...){
  
  # initial_df is a dataframe of starting values 
  # variables by column and site by row
  # e.g. initial_df <- data.frame (S=c(0.9,1,1), I=c(0.1, 0, 0), L=c(0,0,0))
  initial_list <- unlist(initial_df)
    
  # Create the reduced function to use
  ode_func <- factory(...)
  
  # Run the simulation until state variables stabilize
  STABLE <- FALSE 
  state <- data.frame(time=0, t(initial_list))
  tries <- 0
  
  while (!STABLE & tries < max_tries){    
    # Extract the most recent state of the system
    initial_state <- state [nrow(state), ]
    
    # Extend the length of time the model stretches
    max_time <- initial_state$time
    run_times <- times + max_time
  
    # Format the initial state for ode engine
    initial_state <- unlist(initial_state[names(initial_state)!="time"])
  
    # Run the model from the most recent state
    new_state <- as.data.frame(ode(y=initial_state, times=run_times, func=ode_func))
      
    # Combine the old and new results
    # Truncate the repeated row
    state <- rbind (state, new_state[-1,])  
    
    # Test for stability
    STABLE <- const_stability(state, min_rates)
    
    # Update counter for # of tries
    tries <- tries + 1
  }
  
  # Process the output for pretty graphing
  melt_out <- melt_ode_out(state)
  
  return (list(state=state, melt=melt_out))
  
}

# Testing for convergence ####

const_stability <- function (state, min_rates="auto"){
  if (min_rates[1]=="auto"){
    time_step <- state$time[2]-state$time[1]
    
    ranges <- sapply(state[names(state)!="time"], range)
    ranges <- ranges[2,] - ranges[1,]
    
    # Set the sensitivity in proportion to the range of the variable divided by the time step
    sens_const <- 0.0001
    
    min_rates <-  sens_const*ranges/time_step      
  }
  
  # Linearly interpolate first derivative
  diff_state <- sapply(state, diff)
  
  rates <- diff_state[, colnames(diff_state)!="time"]/diff_state[,"time"]
  
  # Check whether rates are below the minimum threshold
  rate_test <- t(apply(rates, MARGIN=1, FUN=function(x){abs(x)<=min_rates}))
  
  # Lazy rule of thumb: if final rate is subthreshold, we can presume system is in equilibrium
  stable <- rate_test[nrow(rate_test),]
  
  
  # Status updates
  print (paste0("Time ", max(state$time), ": stability testing for an approximate constant state"))
  print (stable)
  
  # If all states are stable, we can conclude the simulation
  return (all(stable))
  
}

# Plotting functions ####

# Process the output for pretty graphing
melt_ode_out <- function (state){
  melt_out <- melt (as.data.frame(state), id.vars="time")
  
  var_list <- melt_out$variable
  
  # Relies on [A-Z][0-9]* naming conventions
  # Variable (1 character) followed by site number
  melt_out$var <- substr(var_list, 1,1)
  melt_out$site <- substr(var_list, 2, max(nchar(as.character(var_list))))
  
  return (melt_out)
}



# A clean variable x site plot
sim_facet_plot <- function(melt_out){
  my_plot <- ggplot (melt_out, aes(x=time, y=value, colour=var))+geom_line(size=1)+facet_grid(site~var, scales="free_y")+guides(colour=FALSE)+xlab("Time")+ylab("Value")+theme_bw()
  return(my_plot)
}

