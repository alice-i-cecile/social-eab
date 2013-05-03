# Libraries ####
library(plyr)

# Lee-Ann EAB model, extended to multipatch ####

# Make a prepackaged function for the simulation
eab1_factory <- function (r=1, K=1, allee_threshold=0, transmission_rate=0,   fatality_rate=0,  social_rate=0, transport_cost=0, local_cost=0, norm_strength=0, concern_level=0, dispersal_constant=0){
  
  
  # Define the differential equations. 
  # Static equations can be defined once
  # To avoid recomputing 
  
  # Allee effect
  allee_term <- function (infected){
    if (infected >= allee_threshold){
      return (1)
    } else {
      return (0)
    }
  }
  
  cross_patch_infection_func <- function (healthy_list, infected_list, local_list){
    
    # Infection travelling from patch j to patch i
    cross_infection_ij <- function (healthy_i, infected_j, local_j, dispersal_ij){
      transmission_rate*dispersal_ij*(1-local_j)*healthy_i*infected_j
    }
    
    # Make a list of all possible patch-patch links
    N <- length(healthy_list)
    i <- 1:N
    j <- 1:N
    
    dispersal_df <- expand.grid (i, j)
    names(dispersal_df) <- c("i", "j")
    
    # Dispersal between any two patches is a constant
    dispersal_df$dispersal <- dispersal_constant
    
    # Dispersal between a patch and itself is 0
    dispersal_df[dispersal_df["i"]==dispersal_df["j"], "dispersal"] <- 0
    
    # Find infection rate for each combination
    dispersal_df$healthy_i <- healthy_list[dispersal_df$i]
    dispersal_df$infected_j <- infected_list[dispersal_df$j]
    dispersal_df$local_j <- local_list[dispersal_df$j]
    
    dispersal_df$infection <- mapply(cross_infection_ij, dispersal_df$healthy_i, dispersal_df$infected_j, dispersal_df$local_j, dispersal_df$dispersal)
    
    # Find total infection rate for each patch
    cross_patch_infection <- ddply (dispersal_df, ~i, summarize, total_infection=sum(infection))
    
    infection <- unlist(cross_patch_infection$total_infection)
    names(infection) <- i
    
    return(infection)
    
  }    
  
  # Local vs. transport strategists
  d_local_func <- function (local, infected){
    
    social_influence <- norm_strength*(2*local-1)
    concern <- concern_level*infected
    
    local_payoff <- transport_cost-local_cost+social_influence+concern
    
    #d_local <- (1-local)*m*local*(local_payoff)-local*m*(1-local)*(-local_payoff)
    #d_local <- 2*local*(1-local)*m*local_payoff
    #d_local <- social_rate*local*(1-local)*(transport_cost-local_cost+norm_strength*(2*local-1)+concern_level*infected)
    
    d_local <- social_rate*local*(1-local)*local_payoff
    
    return (d_local)        
  }
  
  
  # Make a reduced function to use in ODE solver
  eab <- function(time, state, parms){
    
    # Extract state
    healthy_id <- substr(names(state), 1, 1)=="S"
    infected_id <- substr(names(state), 1, 1)=="I"
    local_id <- substr(names(state), 1, 1)=="L"
    
    healthy_list <- state[healthy_id]
    infected_list <- state[infected_id] 
    local_list <- state[local_id]
    
    i <- 1:length(healthy_list)
    
    # Find cross-patch infection first
    cross_patch_infection <- cross_patch_infection_func(healthy_list, infected_list, local_list)
    
    # Make new versions of infection functions
    
    # Infection rates
    infection_term <- function (i, healthy, infected){
      infection_rate <- transmission_rate*healthy*infected*allee_term (infected) + cross_patch_infection[i]
      return (infection_rate)
    }
    
    
    # Healthy population
    d_healthy_func <- function (i, healthy, infected){
      d_healthy <- r*healthy*(1-(healthy+infected)/K) - infection_term(i, healthy, infected)
      
      return (d_healthy)
    }
    
    # Infected population
    d_infected_func <- function (i, healthy, infected){
      d_infected <- -fatality_rate*infected + infection_term(i, healthy, infected)
      return (d_infected)
    }
    
    # Find rates of change in state variables
    d_healthy <- mapply(d_healthy_func, i, healthy_list, infected_list)
    d_infected <- mapply(d_infected_func, i, healthy_list, infected_list)
    d_local <- mapply(d_local_func, local_list, infected_list)
    
    # Rename the rates. Optional but adds clarity
    names (d_healthy) <- names(state)[healthy_id]
    names (d_infected) <- names(state)[infected_id]
    names (d_local) <- names(state)[local_id]
    
    # Compile and return rates of change
    out <- list(c(d_healthy, d_infected, d_local))
    
    return (out)            
  }
  
  # Return the reduced function for easy use in ode function
  return (eab)
}

# Lee-Ann EAB model, extended to multipatch ####

# Make a prepackaged function for the simulation
eab2_factory <- function (a=1, b=1, allee_threshold=0, transmission_rate=0,   fatality_rate=0,  social_rate=0, transport_cost=0, local_cost=0, norm_strength=0, concern_level=0, dispersal_constant=0){
  
  
  # Define the differential equations. 
  # Static equations can be defined once
  # To avoid recomputing 
  
  # Allee effect
  allee_term <- function (infected){
    if (infected >= allee_threshold){
      return (1)
    } else {
      return (0)
    }
  }
  
  cross_patch_infection_func <- function (healthy_list, infected_list, local_list){
    
    # Infection travelling from patch j to patch i
    cross_infection_ij <- function (healthy_i, infected_j, local_j, dispersal_ij){
      transmission_rate*dispersal_ij*(1-local_j)*healthy_i*infected_j
    }
    
    # Make a list of all possible patch-patch links
    N <- length(healthy_list)
    i <- 1:N
    j <- 1:N
    
    dispersal_df <- expand.grid (i, j)
    names(dispersal_df) <- c("i", "j")
    
    # Dispersal between any two patches is a constant
    dispersal_df$dispersal <- dispersal_constant
    
    # Dispersal between a patch and itself is 0
    dispersal_df[dispersal_df["i"]==dispersal_df["j"], "dispersal"] <- 0
    
    # Find infection rate for each combination
    dispersal_df$healthy_i <- healthy_list[dispersal_df$i]
    dispersal_df$infected_j <- infected_list[dispersal_df$j]
    dispersal_df$local_j <- local_list[dispersal_df$j]
    
    dispersal_df$infection <- mapply(cross_infection_ij, dispersal_df$healthy_i, dispersal_df$infected_j, dispersal_df$local_j, dispersal_df$dispersal)
    
    # Find total infection rate for each patch
    cross_patch_infection <- ddply (dispersal_df, ~i, summarize, total_infection=sum(infection))
    
    infection <- unlist(cross_patch_infection$total_infection)
    names(infection) <- i
    
    return(infection)
    
  }    
  
  # Local vs. transport strategists
  d_local_func <- function (local, infected){
    
    social_influence <- norm_strength*(2*local-1)
    concern <- concern_level*infected
    
    local_payoff <- transport_cost-local_cost+social_influence+concern
    
    #d_local <- (1-local)*m*local*(local_payoff)-local*m*(1-local)*(-local_payoff)
    #d_local <- 2*local*(1-local)*m*local_payoff
    #d_local <- social_rate*local*(1-local)*(transport_cost-local_cost+norm_strength*(2*local-1)+concern_level*infected)
    
    d_local <- social_rate*local*(1-local)*local_payoff
    
    return (d_local)        
  }
  
  
  # Make a reduced function to use in ODE solver
  eab <- function(time, state, parms){
    
    # Extract state
    healthy_id <- substr(names(state), 1, 1)=="S"
    infected_id <- substr(names(state), 1, 1)=="I"
    local_id <- substr(names(state), 1, 1)=="L"
    
    healthy_list <- state[healthy_id]
    infected_list <- state[infected_id] 
    local_list <- state[local_id]
    
    i <- 1:length(healthy_list)
    
    # Find cross-patch infection first
    cross_patch_infection <- cross_patch_infection_func(healthy_list, infected_list, local_list)
    
    # Make new versions of infection functions
    
    # Infection rates
    infection_term <- function (i, healthy, infected){
      infection_rate <- transmission_rate*healthy*infected*allee_term (infected) + cross_patch_infection[i]
      return (infection_rate)
    }
    
    
    # Healthy population
    d_healthy_func <- function (i, healthy, infected){
      d_healthy <- a*(healthy+infected) - b*(healthy+infected)*healthy- infection_term(i, healthy, infected)
      
      return (d_healthy)
    }
    
    # Infected population
    d_infected_func <- function (i, healthy, infected){
      d_infected <- -fatality_rate*infected - b*(healthy+infected)*infected+ infection_term(i, healthy, infected)
      return (d_infected)
    }
    
    # Find rates of change in state variables
    d_healthy <- mapply(d_healthy_func, i, healthy_list, infected_list)
    d_infected <- mapply(d_infected_func, i, healthy_list, infected_list)
    d_local <- mapply(d_local_func, local_list, infected_list)
    
    # Rename the rates. Optional but adds clarity
    names (d_healthy) <- names(state)[healthy_id]
    names (d_infected) <- names(state)[infected_id]
    names (d_local) <- names(state)[local_id]
    
    # Compile and return rates of change
    out <- list(c(d_healthy, d_infected, d_local))
    
    return (out)            
  }
  
  # Return the reduced function for easy use in ode function
  return (eab)
}


