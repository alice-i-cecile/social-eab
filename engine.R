# Libraries ####
library(deSolve)
library(reshape2)
library(ggplot2)

# Wrapper for ode simulation of EAB models ####

sim <- function(initial_df, times=seq(from=0, to=10, by=0.1), factory=eab_factory, ...){
  
  # initial_df is a dataframe of starting values 
  # variables by column and site by row
  # e.g. initial_df <- data.frame (S=c(0.9,1,1), I=c(0.1, 0, 0), L=c(0,0,0))
  initial_list <- unlist(initial_df)
  
  # Create the reduced function to use
  ode_func <- factory(...)
  
  # Run the simulation
  raw_out <- ode(y=initial_list, times=times, func=ode_func)
  
  # Process the output for pretty graphing
  melt_out <- melt (as.data.frame(raw_out), id.vars="time")
  
  var_list <- melt_out$variable
  
  # Relies on [A-Z][0-9]* naming conventions
  # Variable (1 character) followed by site number
  melt_out$var <- substr(var_list, 1,1)
  melt_out$site <- substr(var_list, 2, max(nchar(as.character(var_list))))
  
  return (melt_out)
  
}

# Plotting functions ####

sim_facet_plot <- function(melt_out){
  my_plot <- ggplot (melt_out, aes(x=time, y=value, colour=var))+geom_line()+facet_wrap(var~site, scales="free_y")+theme_bw()
  return(my_plot)
}

# Testing ####
initial_test <- data.frame (S=c(1), I=c(0), L=c(1))
times_test <- seq(from=0, to=1, by=0.1)


foo <- sim (initial_test, times_test)
print(sim_facet_plot(foo))