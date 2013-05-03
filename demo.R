# Load source files ####
source ("models.R")
source ("engine.R")

# Testing EAB1 ####
# Parameters
initial_test <- data.frame (S=c(0.8, 1, 1), I=c(0.2,0, 0), L=c(0.5,0.5, 0.5))
times_test <- seq(from=0, to=10, by=0.1)

# Run the model
sim_results <- sim (initial_test, times_test, eab1_factory, max_tries=100, r=3, K=1, allee_threshold=0, transmission_rate=1, fatality_rate=0.1,  social_rate=0.5, transport_cost=0, local_cost=0, norm_strength=0, concern_level=1, dispersal_constant=0.1)

# Graphing
print(sim_facet_plot(sim_results$melt))

# Testing EAB2 ####
# Parameters
initial_test <- data.frame (S=c(0.5), I=c(0.5), L=c(0.5))
times_test <- seq(from=0, to=10, by=1)

# Run the model
sim_results <- sim (initial_test, times_test, eab2_factory,max_tries=100, a=3, b=1, allee_threshold=0, transmission_rate=1.5, fatality_rate=0,  social_rate=0.5, transport_cost=0, local_cost=0.5, norm_strength=0, concern_level=1, dispersal_constant=0.1)

# Graphing
print(sim_facet_plot(sim_results$melt))


# Modified phase diagram ####
# Parameters
initial_test <- data.frame (S=c(4985, 5000), I=c(15,0), L=c(0.1,0.1))
times_test <- seq(from=0, to=10, by=0.1)

# Run the model
sim_results <- sim (initial_test, times_test, eab1_factory, max_tries=2, r=1, K=5000, allee_threshold=0, transmission_rate=0.1, fatality_rate=1,  social_rate=0.1, transport_cost=5, local_cost=6.75, norm_strength=0.1, concern_level=0.1, dispersal_constant=0)

# Graphing
print(sim_facet_plot(sim_results$melt))


# Parameters
initial_test <- data.frame (S=c(0.9, 1), I=c(0.1,0), L=c(0.1,0.1))
times_test <- seq(from=0, to=10, by=0.1)

# Run the model
sim_results <- sim (initial_test, times_test, eab2_factory, max_tries=10, a=1, b=1, allee_threshold=0, transmission_rate=1.5, fatality_rate=1,  social_rate=0.1, transport_cost=5, local_cost=6.75, norm_strength=0.1, concern_level=0.1, dispersal_constant=0.1)

sim_results <- sim (initial_test, times_test, eab1_factory, max_tries=10, r=1, K=1, allee_threshold=0, transmission_rate=1.5, fatality_rate=1,  social_rate=0.1, transport_cost=5, local_cost=6.75, norm_strength=0.1, concern_level=0.1, dispersal_constant=0.1)


# Graphing
print(sim_facet_plot(sim_results$melt))
print(sim_facet_plot(sim_results$melt)+ylim(0,100))
