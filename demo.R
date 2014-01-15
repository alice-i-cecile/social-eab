# Load source files ####
source ("models.R")
source ("engine.R")

# Default EAB ####
# Parameters
initial_test <- data.frame (S=c(4985, 5000), I=c(15,0), L=c(0.1,0.1))
times_test <- seq(from=0, to=10, by=0.1)

# Run the model
sim_results <- sim (initial_test, times_test, eab_factory, max_tries=100, r=0.06, K=5000, allee_threshold=1, transmission_rate=6.5e-4, fatality_rate=3,  social_rate=0.1, transport_cost=5, local_cost=6.75, norm_strength=0.1, concern_level=0.1, dispersal_constant=0.05)

# Graphing
print(sim_facet_plot(sim_results$melt))

# Sandbox EAB ####
# Parameters
initial_test <- data.frame (S=c(0.85, 1), I=c(0.15,0), L=c(0.1,0.1))
times_test <- seq(from=0, to=10, by=0.1)

# Run the model
sim_results <- sim (initial_test, times_test, eab_factory, max_tries=100, r=0.06, K=1, allee_threshold=1, transmission_rate=6.5e-4, fatality_rate=1/3,  social_rate=0.1, transport_cost=5, local_cost=6.75, norm_strength=0.1, concern_level=0.1, dispersal_constant=0.05)

# Graphing
print(sim_facet_plot(sim_results$melt))


# Phase diagram for EAB ####

# Setup
search_space <- data.frame(concern_level=c(0, 1), local_cost=c(0, 6.75))

search_points <- search_p(search_space=search_space, npoints=64^2, method="grid")

static_parm <- list(r=0.06, K=5000, allee_threshold=1, transmission_rate=6.5e-4, fatality_rate=3,  social_rate=0.1, transport_cost=5, norm_strength=0.1,  dispersal_constant=0.05)

initial_test <- data.frame (S=c(4985, 5000), I=c(15,0), L=c(0.1,0.1))
times_test <- seq(from=0, to=10, by=0.1)

# Search the parameter space
phase_df <- phase_search_points(search_points=search_points, phase_func=phase_eab, static_parm=static_parm, initial_df=initial_test, times=times_test, factory=eab_factory, max_tries=100)

# Plotting  

ggplot(data=phase_df, aes(x=local_cost, y=concern_level, fill=phase)) + geom_tile()

phase_plot <- ggplot(data=phase_df, aes(x=local_cost, y=concern_level, fill=phase)) + geom_tile()+xlab("Cost of local firewood (Cl)") + ylab("Concern level (f)") + scale_fill_manual(values=c("grey75", "grey25"), name="Outcome")+theme_bw()
print(phase_plot)

phase_plot + guides(fill=FALSE)
