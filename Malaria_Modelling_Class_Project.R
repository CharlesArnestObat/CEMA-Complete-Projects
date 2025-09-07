###### The following exercise ####
##Network Analysis and SIS Epidemic Model##

####Project Description####

##This project introduces the intersection of network analysis and epidemic modeling through a hands-on exploration of disease transmission dynamics. Using a dataset of participants from various countries attending a Malaria Modelling modular program in Nairobi, Kenya, participants learn to construct social networks based on participants interaction patterns and simulate disease spread using a Susceptible-Infected-Susceptible (SIS) model.

##The project aims to demonstrate how network structure influences epidemic outcomes by modeling connections between individuals from the same countries, with Kenya serving as a hub connecting participants from different nations. The aim is to explore key epidemiological concepts including transmission rates (β), recovery rates (γ), and intervention strategies like wearing masks and identification of superspreaders.

##Through this practical exercise, participants should gain a hands-on experience with network visualization, epidemic simulation, and quantitative analysis of public health interventions. The goals is to emphasize applications by examining scenarios such as the impact of protective measures introduced mid-outbreak and the role of highly connected individuals in accelerating disease spread.

##Exercises are presented at the end of this tutorial

##############################################################################

#Packages for this work
library(igraph) # For graphs
library(ggplot2) # For plots
library(dplyr) # For manipulation
library(RColorBrewer) # For colours

##Create the dataset - participants, and their country of origin.
data <- data.frame( Participant = c("KAM", "SIN", "RAH", "DAW", "LUC", "TUH", "BAS", "MBA", "MUS", "DAV", "CHA", "EMM", "LYD", "VAL", "GOD", "AGN", "ROS", "FRA", "JAM", "DGB", "ELI", "WAN", "MUT", "GEO", "GEF", "SAH", "PAT", "TAB", "BRI", "STA", "OUM", "STE", "THU", "EMI", "CAM", "ZEN", "BIL", "JOH", "ALI", "PUR", "WAI", "MIL", "BRL", "TRI"),  Country = c("Cameroon", "Ethiopia", "Tanzania", "Ethiopia", "Zambia", "Uganda", "Nigeria", "Cameroon", "Zimbabwe", "Zambia", "Kenya", "Kenya", "Kenya", "Kenya", "Malawi", "Malawi", "Tanzania", "Tanzania", "Kenya", "Ethiopia", "Tanzania", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Switzerland", "Switzerland", "Ghana", "Switzerland", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya"))

#####PART 1. NETWORK CONSTRUCION ####

##Remember Networks are based on matrices - adjacency matrix defines who is connected to whom.
##Create an adjacency matrix based on connection rules: 
#a) People from same country are connected; 
#b) Non-Kenyans are also connected to Kenyans. 
create_network <- function(data) {
  n <- nrow(data)
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- data$Participant
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
##Same country connect
        if(data$Country[i] == data$Country[j]) {
          adj_matrix[i,j] <- 1
        }
##Non-Kenyans connect to Kenyans
        else if(data$Country[i] != "Kenya" & data$Country[j] == "Kenya") {
          adj_matrix[i,j] <- 1
        }
        else if(data$Country[i] == "Kenya" & data$Country[j] != "Kenya") {
          adj_matrix[i,j] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}

##Create the network and add attributes - what is an attribute?

adj_matrix <- create_network(data) ### Explore the adjacency 
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected") #What does undirected mean? How does it differ from a directed graph?

##Add attributes
V(g)$country <- data$Country
V(g)$name <- data$Participant

##Create and visualize the network
plot_network <- function(g, data) {
##Color nodes by country
countries <- unique(data$Country)
colors <- RColorBrewer::brewer.pal(min(length(countries), 11), "Spectral")
country_colors <- setNames(colors[1:length(countries)], countries)
  
##Set vertex colors using the country attribute
V(g)$color <- country_colors[V(g)$country]
  
##Set the graph layout - there are many layouts - please explore them
set.seed(123)
layout <- layout_with_fr(g)
  
##Plot
plot(g,
     layout = layout,
     vertex.size = 8,
     vertex.label = V(g)$name,  # Explicitly use vertex names - what if you don't want to have people's names on the vertices? How would you modify this?
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = "gray70",
     edge.width = 0.8,
     main = "Participant Interaction Network")
  
##Add legend - you can change position or text size
legend("topright", legend = countries, fill = country_colors, cex = 0.8, title = "Country")
}

###Plot the network
plot_network(g, data)


##############################################################################
###Now you have a network - what is it showing? Why circle? Are there other things you can explore from this?

###You can stop here - and explore the network parameters - density, degree (in/out degree), betweenness and other centrality metrics.

####What type of questions can one explore with network analysis####
##What is my outcome of analysis?
# • Node-level
# • Edge-level
# • Network structure and topology (as a whole)
# • Process (sequence of events
##Outcome variables
# • Traits of nodes (e.g., Disease state)– Do nodes of different classes have different attributes?– Do connectivity and network position determine infection risk?
#   • Edge values– What determines whether two nodes are linked?– Do linked nodes share certain attributes?
#   • Network topology– What determines the structure of networks?– Does infection influence network structure?
#   • Transmission processes– Does the spread of a pathogen depend on network structure?– Simulation modeling

####Part 2 -SIS Epidemic Modelling ####
sis_model <- function(adj_matrix, beta, gamma, initial_infected, days, superspreaders = NULL, mask_day = NULL, mask_efficacy = 0.3) {
  n <- nrow(adj_matrix)
  participants <- rownames(adj_matrix)
  
##Initialize states: 0 = Susceptible, 1 = Infected
  states <- matrix(0, nrow = days + 1, ncol = n)
  colnames(states) <- participants
  
##Set initial infected
  states[1, initial_infected] <- 1
  
##Track daily infections for incidence calculation
  daily_new_infections <- rep(0, days + 1)
  
##Superspreaders - who is a suer spreader?
  ss_multiplier <- rep(1, n)
  if(!is.null(superspreaders)) {
    ss_idx <- which(participants %in% superspreaders)
    ss_multiplier[ss_idx] <- 5  # 5x transmission rate
  }
  
  for(day in 1:days) {
    current_state <- states[day, ]
    new_state <- current_state
    new_infections <- 0
    
##Apply Intervention - masking effect if applicable
    current_beta <- beta
    if(!is.null(mask_day) && day >= mask_day) {
      current_beta <- beta * (1 - mask_efficacy)
    }
    
    for(i in 1:n) {
      if(current_state[i] == 0) {  # If susceptible
##Calculate infection probability from infected neighbors
    infected_neighbors <- which(current_state == 1 & adj_matrix[i, ] == 1)
    if(length(infected_neighbors) > 0) {
##Calculate total force of infection
      total_force <- 0
      for(neighbor in infected_neighbors) {
        total_force <- total_force + current_beta * ss_multiplier[neighbor]
          }
          
##Probability of infection (1 - probability of avoiding all infections)
      prob_infection <- 1 - exp(-total_force)
      if(runif(1) < prob_infection) {
        new_state[i] <- 1
        new_infections <- new_infections + 1
        }
        }
      } else {  ##If infected
##Recovery
        if(runif(1) < gamma) {
          new_state[i] <- 0
        }
      }
    }
    states[day + 1, ] <- new_state
    daily_new_infections[day + 1] <- new_infections
  }
  
  return(list(states = states, daily_new_infections = daily_new_infections))
}

###3 - SIMULATIONS

# Set parameters
beta <- 0.101
gamma <- 0.1
days <- 30

##Find participant indices
zeinabu_idx <- which(data$Participant == "ZEN")
agnes_idx <- which(data$Participant == "AGN")
initial_infected <- c("ZEN", "AGN")

##Scenario 1: Basic SIS model
set.seed(42)
result1 <- sis_model(adj_matrix, beta, gamma, initial_infected, days) #Explore

##Scenario 2: With masks from day 3
set.seed(42)
result2 <- sis_model(adj_matrix, beta, gamma, initial_infected, days, mask_day = 3)

##Scenario 3: With superspreaders
superspreaders <- c("BAS", "BIL", "MUT", "PUR") #Are these the superspreaders? Who is a superspreader? What criteria do we use to identify a superspreader? Hint: Network metrics. 
set.seed(42)
result3 <- sis_model(adj_matrix, beta, gamma, initial_infected, days, superspreaders = superspreaders)

####Part 4. ANALYSIS AND RESULTS ####

##Epidemic Results

# a) Number of infected on days 2, 4, 8
infected_day2 <- sum(result1$states[3, ])  # Day 2 (row 3, 0-indexed)
infected_day4 <- sum(result1$states[5, ])  # Day 4
infected_day5 <- sum(result1$states[5, ])  # Day 5
infected_day8 <- sum(result1$states[9, ])  # Day 8

cat("a) Number of infected individuals:\n")
cat("   Day 2:", infected_day2, "people\n")
cat("   Day 4:", infected_day4, "people\n")
cat("   Day 8:", infected_day8, "people\n\n")

# b) Incidence at day 8
incidence_day8 <- result1$daily_new_infections[5]
cat("b) Incidence on Day 5:", incidence_day8, "new infections\n\n")

# c) Infections averted with masks
infected_day5_masks <- sum(result2$states[5, ])
infections_averted <- infected_day5 - infected_day5_masks
cat("c) With masks from Day 5:\n")
cat("   Infections on Day 5 with masks:", infected_day5_masks, "\n")
cat("   Infections averted:", infections_averted, "\n\n")

# d) Network metrics
cat("d) Network Metrics:\n")

# i) Density
density_val <- edge_density(g)
cat("   i) Network density:", round(density_val, 3), "\n")

# ii) Superspreader analysis
cat("   ii) Superspreader analysis:\n")
infected_normal <- sum(result1$states[6, ])
infected_superspreaders <- sum(result3$states[6, ])
cat("       Without superspreaders: ", infected_normal, " infected\n")
cat("       With superspreaders: ", infected_superspreaders, " infected\n")
cat("       Additional infections from superspreaders: ", infected_superspreaders - infected_normal, "\n")

## iii) Effect of interaction rate changes
cat("   iii) Effect of changing interaction rates:\n")
betas <- c(0.01, 0.02, 0.03, 0.04, 0.05)
final_infected <- numeric(length(betas))

for(i in seq_along(betas)) {
  set.seed(42)
  temp_result <- sis_model(adj_matrix, betas[i], gamma, initial_infected, days)
  final_infected[i] <- sum(temp_result$states[6, ])
}

beta_results <- data.frame(Beta = betas, Final_Infected = final_infected)
print(beta_results)

####Part 5. VISUALIZATION OF EPIDEMIC PROGRESSION ####
plot_epidemic <- function(states) {
  days <- nrow(states) - 1
  infected_counts <- rowSums(states)
  
  df <- data.frame(
    Day = 0:days,
    Infected = infected_counts
  )
  
  ggplot(df, aes(x = Day, y = Infected)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(title = "Epidemic Progression",
         x = "Day",
         y = "Number of Infected Individuals") +
    theme_minimal() +
    scale_x_continuous(breaks = 0:days)
}

##Plot epidemic progression
print(plot_epidemic(result1$states))

##Comparison plot
compare_df <- data.frame(
  Day = rep(0:days, 3),
  Infected = c(rowSums(result1$states), rowSums(result2$states), rowSums(result3$states)),
  Scenario = rep(c("Basic", "With Masks", "With Superspreaders"), each = days + 1)
)

comparison_plot <- ggplot(compare_df, aes(x = Day, y = Infected, color = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Comparison of Epidemic Scenarios",
       x = "Day",
       y = "Number of Infected Individuals") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:days) +
  scale_color_manual(values = c("red", "blue", "orange"))

print(comparison_plot)

##Enhanced infections averted visualization as a shaded area
plot_infections_averted <- function(result_no_control, result_with_control, intervention_name = "Intervention") {
  days <- nrow(result_no_control$states) - 1
  
##Calculate infections for both scenarios
no_control <- rowSums(result_no_control$states)
with_control <- rowSums(result_with_control$states)
  
##Calculate total infections averted
total_averted <- sum(no_control) - sum(with_control)
percentage_averted <- round((total_averted / sum(no_control)) * 100, 1)
  
##Create data frame
df <- data.frame(
  Day = 0:days,
  No_Control = no_control,
  With_Control = with_control
  )
##Plot
p <- ggplot(df, aes(x = Day)) +
##Shaded area between curves
geom_ribbon(aes(ymin = With_Control, ymax = No_Control),
            fill = "lightblue", alpha = 0.7) +
##Lines for both scenarios
geom_line(aes(y = No_Control, color = "Without Control"), size = 1.2) +
geom_line(aes(y = With_Control, color = paste("With", intervention_name)), size = 1.2) +
##Add text annotation in shaded area
  annotate("text", x = days/2, y = max(no_control)/2, 
           label = paste("Infections Averted:\n", total_averted, "cases\n(", percentage_averted, "%)", sep = ""), size = 4, fontface = "bold", color = "darkblue") +
  labs(title = paste("Impact of", intervention_name, "on Disease Transmission"), x = "Day", y = "Number of Infected Individuals", color = "Scenario") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue")) +
  theme(legend.position = "bottom")
return(p)
}

##Show infections averted with masks
print(plot_infections_averted(result1, result2, "Masks"))

##Summary
cat("Network has", vcount(g), "participants and", ecount(g), "connections\n")
cat("Density:", round(density_val, 3), "- indicates a", 
    ifelse(density_val > 0.5, "highly", "moderately"), "connected network\n")

##############################################################################
##EXERCISE AND SOLUTIONS##

##Exercise 1: Network Exploration and Visualization

cat("\n=== EXERCISE 1: NETWORK EXPLORATION ===\n")

## a). Modify the network to show node sizes proportional to degree centrality
plot_network_degree <- function(g, data) {
##Calculate degree centrality
  degrees <- degree(g)
  
##Color nodes by country
  countries <- unique(data$Country)
  colors <- RColorBrewer::brewer.pal(min(length(countries), 11), "Spectral")
  country_colors <- setNames(colors[1:length(countries)], countries)
  V(g)$color <- country_colors[V(g)$country]
  
##Layout
  set.seed(123)
  layout <- layout_with_fr(g)
  
##Plot with node sizes proportional to degree
  plot(g, 
       layout = layout,
       vertex.size = degrees * 2,  # Scale degree for visibility
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 0.8,
       main = "Network with Node Sizes Proportional to Degree Centrality")
  
  legend("topright", legend = countries, fill = country_colors, cex = 0.8, title = "Country")
}

plot_network_degree(g, data)

###Exercise 1.1- What is betweenness centrality. Redraw the network with node sizes proportional to the node betweenness. How does this differ from degree centrality?

## b). Calculate and compare degree centrality for Kenya vs other countries
degrees <- degree(g)
kenya_degrees <- degrees[data$Country == "Kenya"]
other_degrees <- degrees[data$Country != "Kenya"]

cat("b) Degree Centrality Comparison:\n")
cat("   Kenya participants - Mean degree:", round(mean(kenya_degrees), 2), "\n")
cat("   Other countries - Mean degree:", round(mean(other_degrees), 2), "\n")
cat("   Kenya participants are", round(mean(kenya_degrees)/mean(other_degrees), 1), "times more connected\n\n")

##Exercise 1.2 - Calculate and compare betweenness centrality for Kenya vs other countries

## c). Identify top 3 most connected individuals
top_connected <- sort(degrees, decreasing = TRUE)[1:3]
cat("c) Top 3 most connected individuals:\n")
for(i in 1:3) {
  participant <- names(top_connected)[i]
  connections <- top_connected[i]
  country <- data$Country[data$Participant == participant]
  cat("   ", i, ".", participant, "(", country, ") -", connections, "connections\n")
}


##Exercise 1.3 - Identify least 3 connected individuals


##Exercise 2: Parameter Sensitivity Analysis

cat("\n=== EXERCISE 2: PARAMETER SENSITIVITY ===\n")

## a) & b). Run SIS model with different β values and plot
betas <- seq(0.001, 0.2, by = 0.01)
final_infected_beta <- numeric(length(betas))

for(i in seq_along(betas)) {
  set.seed(42)
  temp_result <- sis_model(adj_matrix, betas[i], gamma, initial_infected, days)
  final_infected_beta[i] <- sum(temp_result$states[days + 1, ])
}

# Plot β sensitivity
beta_df <- data.frame(Beta = betas, Final_Infected = final_infected_beta)
beta_plot <- ggplot(beta_df, aes(x = Beta, y = Final_Infected)) +
  geom_line(color = "red", size = 1.2) +
  geom_point(color = "red", size = 2) +
  labs(title = "Effect of Transmission Rate (β) on Final Infection Count",
       x = "Transmission Rate (β)",
       y = "Final Number of Infected") +
  theme_minimal()
print(beta_plot)

##Exercise 2.1 - Why does the infection plateau? What's the epidemiological phenomenon?

## c). Vary γ while keeping β constant
gammas <- seq(0.05, 0.5, by = 0.05)
final_infected_gamma <- numeric(length(gammas))

for(i in seq_along(gammas)) {
  set.seed(42)
  temp_result <- sis_model(adj_matrix, 0.4, gammas[i], initial_infected, days)
  final_infected_gamma[i] <- sum(temp_result$states[days + 1, ])
}

# Plot γ sensitivity
gamma_df <- data.frame(Gamma = gammas, Final_Infected = final_infected_gamma)
gamma_plot <- ggplot(gamma_df, aes(x = Gamma, y = Final_Infected)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "blue", size = 2) +
  labs(title = "Effect of Recovery Rate (γ) on Final Infection Count",
       x = "Recovery Rate (γ)",
       y = "Final Number of Infected") +
  theme_minimal()
print(gamma_plot)

## d). Calculate R₀ and plot
R0_values <- betas / gamma  # Using gamma = 0.1
R0_df <- data.frame(Beta = betas, R0 = R0_values, Final_Infected = final_infected_beta)

R0_plot <- ggplot(R0_df, aes(x = R0, y = Final_Infected)) +
  geom_line(color = "purple", size = 1.2) +
  geom_point(color = "purple", size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.2, y = max(final_infected_beta)/2, label = "R₀ = 1\n(Threshold)", color = "red") +
  labs(title = "Basic Reproduction Number (R₀) vs Final Infections",
       x = "Basic Reproduction Number (R₀)",
       y = "Final Number of Infected") +
  theme_minimal()
print(R0_plot)

cat("d) R₀ Analysis: Epidemic threshold occurs at R₀ = 1\n")
cat("   When R₀ < 1: Disease dies out\n")
cat("   When R₀ > 1: Disease persists and spreads\n\n")

##Exercise 2.2 - What is R₀? How would you use the knowledge of R₀ for disease control? What's another way tomcalculate R₀ as Billy showed in class? 

##Exercise 3: Intervention Timing Analysis

cat("\n=== EXERCISE 3: INTERVENTION TIMING ===\n")

## a) & b). Implement masks on different days and calculate infections averted
intervention_days <- 1:6
infections_averted_list <- numeric(length(intervention_days))
mask_efficacy <- 0.4

for(i in seq_along(intervention_days)) {
  set.seed(42)
  result_mask <- sis_model(adj_matrix, beta, gamma, initial_infected, days, 
                           mask_day = intervention_days[i], mask_efficacy = mask_efficacy)
  
  # Calculate total infections over entire period
  total_no_mask <- sum(rowSums(result1$states))
  total_with_mask <- sum(rowSums(result_mask$states))
  infections_averted_list[i] <- total_no_mask - total_with_mask
}

## c). Create bar chart
timing_df <- data.frame(
  Day = intervention_days,
  Infections_Averted = infections_averted_list
)

timing_plot <- ggplot(timing_df, aes(x = factor(Day), y = Infections_Averted)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = Infections_Averted), vjust = -0.3) +
  labs(title = "Infections Averted by Day of Mask Intervention",
       subtitle = "Earlier intervention = More infections prevented",
       x = "Day of Intervention",
       y = "Total Infections Averted") +
  theme_minimal()
print(timing_plot)

cat("c) Intervention Timing Results:\n")
for(i in seq_along(intervention_days)) {
  cat("   Day", intervention_days[i], "intervention:", infections_averted_list[i], "infections averted\n")
}

cat("\nd) Early intervention is more effective because:\n")
cat("   - Smaller infected population to control\n")
cat("   - Prevents exponential growth phase\n")
cat("   - Delays can lead to overwhelming healthcare systems\n\n")

##Exercise 4: Superspreader Impact Assessment

cat("\n=== EXERCISE 4: SUPERSPREADER ANALYSIS ===\n")

## a). Identify top 6 most connected individuals
top6_connected <- sort(degrees, decreasing = TRUE)[1:6]
top6_names <- names(top6_connected)

cat("a) Top 6 most connected individuals:\n")
for(i in 1:6) {
  cat("   ", i, ".", top6_names[i], "-", top6_connected[i], "connections\n")
}

##Exercise 31. - Who are the top 6 most connected individuals abased on betweenness? Do they differ from those based on degree centrality? Comment on whether they are different or the same - why?

## b) & c). Compare superspreader scenarios
# Scenario 1: Top 6 most connected as superspreaders
set.seed(42)
result_ss_connected <- sis_model(adj_matrix, beta, gamma, initial_infected, days, superspreaders = top6_names)

# Scenario 2: Random 6 people as superspreaders
set.seed(123)
random_ss <- sample(data$Participant, 6)
set.seed(42)
result_ss_random <- sis_model(adj_matrix, beta, gamma, initial_infected, days, superspreaders = random_ss)

## d). Calculate superspreader effects
baseline_infections <- sum(rowSums(result1$states))
connected_ss_infections <- sum(rowSums(result_ss_connected$states))
random_ss_infections <- sum(rowSums(result_ss_random$states))

connected_effect <- round(((connected_ss_infections - baseline_infections) / baseline_infections) * 100, 1)
random_effect <- round(((random_ss_infections - baseline_infections) / baseline_infections) * 100, 1)

cat("\nb-d) Superspreader Effect Analysis:\n")
cat("   Baseline (no superspreaders):", baseline_infections, "total infections\n")
cat("   Most connected as superspreaders:", connected_ss_infections, "total infections (", connected_effect, "% increase)\n")
cat("   Random people as superspreaders:", random_ss_infections, "total infections (", random_effect, "% increase)\n")

##Exercise 3.1 - Calculate the superspreader effect based on betweenness centrality, and compare most connected as superspreader and random people as superspreader vs baseline? Do they differ from when using degree centrality?

## e). Design targeted intervention
# Target the most connected individuals for intervention (e.g., isolation)
# Simple approach: remove top 3 most connected from network
adj_matrix_targeted <- adj_matrix
top3_indices <- which(rownames(adj_matrix) %in% names(top_connected)[1:3])

# Zero out connections for top 3 (simulate isolation)
adj_matrix_targeted[top3_indices, ] <- 0
adj_matrix_targeted[, top3_indices] <- 0

set.seed(42)
result_targeted <- sis_model(adj_matrix_targeted, beta, gamma, initial_infected, days)
targeted_infections <- sum(rowSums(result_targeted$states))

cat("\ne) Targeted Intervention (Isolate top 3 connected):\n")
cat("   With targeted intervention:", targeted_infections, "total infections\n")
cat("   Infections prevented:", baseline_infections - targeted_infections, "\n")
cat("   Effectiveness:", round(((baseline_infections - targeted_infections) / baseline_infections) * 100, 1), "% reduction\n\n")

##Exercise 3.2 - How many infections are prevented if we consider betweenness centrality fo isolation? 

##Exercise 5: Network Structure Modification

cat("\n=== EXERCISE 5: NETWORK STRUCTURE COMPARISON ===\n")

## a). Create small-world network
create_small_world <- function(data, prob = 0.1) {
  n <- nrow(data)
  adj_matrix_sw <- matrix(0, n, n)
  rownames(adj_matrix_sw) <- colnames(adj_matrix_sw) <- data$Participant
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j && runif(1) < prob) {
        adj_matrix_sw[i,j] <- 1
      }
    }
  }
  return(adj_matrix_sw)
}

# Create small-world network
set.seed(456)
adj_matrix_sw <- create_small_world(data, 0.1)
g_sw <- graph_from_adjacency_matrix(adj_matrix_sw, mode = "undirected")

# Create random network for comparison
set.seed(789)
g_random <- erdos.renyi.game(nrow(data), edge_density(g), directed = FALSE)
V(g_random)$name <- data$Participant
adj_matrix_random <- as_adjacency_matrix(g_random, sparse = FALSE)

## b). Compare network structures
networks <- list(
  Original = list(graph = g, adj_matrix = adj_matrix),
  SmallWorld = list(graph = g_sw, adj_matrix = adj_matrix_sw),
  Random = list(graph = g_random, adj_matrix = adj_matrix_random)
)

##Run simulations on all networks
results_comparison <- list()
for(net_name in names(networks)) {
  set.seed(42)
  results_comparison[[net_name]] <- sis_model(
    networks[[net_name]]$adj_matrix, 0.4, 0.1, initial_infected, days
  )
}

##Calculate metrics
cat("b) Network Structure Comparison:\n")
for(net_name in names(networks)) {
  net_graph <- networks[[net_name]]$graph
  density <- round(edge_density(net_graph), 3)
  avg_path <- round(mean_distance(net_graph), 2)
  final_infections <- sum(results_comparison[[net_name]]$states[days + 1, ])
  
##Find time to peak
  daily_infections <- rowSums(results_comparison[[net_name]]$states)
  peak_day <- which.max(daily_infections) - 1
  
  cat("   ", net_name, "Network:\n")
  cat("     i) Density:", density, "\n")
  cat("     ii) Average path length:", avg_path, "\n")
  cat("     iii) Final infections:", final_infections, "\n")
  cat("     iv) Time to peak:", peak_day, "days\n\n")
}

##Visualize network structure comparison
comparison_epidemic_plot <- function(results_list) {
  all_data <- data.frame()
  
  for(net_name in names(results_list)) {
    infected_counts <- rowSums(results_list[[net_name]]$states)
    temp_df <- data.frame(
      Day = 0:days,
      Infected = infected_counts,
      Network = net_name
    )
    all_data <- rbind(all_data, temp_df)
  }
  
  ggplot(all_data, aes(x = Day, y = Infected, color = Network)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(title = "Epidemic Progression Across Different Network Structures",
         subtitle = "Impact of network topology on disease spread",
         x = "Day",
         y = "Number of Infected Individuals") +
    theme_minimal() +
    scale_color_manual(values = c("red", "blue", "green")) +
    theme(legend.position = "bottom")
}

print(comparison_epidemic_plot(results_comparison))

##############################################################################
##FINAL EPIDEMIOLOGICAL QUESTION##

##Exercise 4 - General Epidemiological Modelling Questions


cat("\n=== EPIDEMIOLOGICAL RESEARCH QUESTION ===\n")

##Exercise 4.1 -  Read the following 
cat("Based on your analysis, answer the following research question:\n\n")
cat("QUESTION: You are a public health official planning intervention strategies\n")
cat("for a respiratory disease outbreak in a diverse international community.\n")
cat("Given limited resources, you must choose ONE of the following strategies:\n\n")

cat("Strategy A: Implement universal masking for everyone starting on Day 3\n")
cat("Strategy B: Identify and isolate the top 3 most connected individuals immediately\n")
cat("Strategy C: Implement random testing and isolation of 6 individuals on Day 1\n\n")

###############################################################################
####END OF Questiions - BONUS IMPLEMENTATION BELOW
###############################################################################


###Calculate effectiveness of each strategy
cat("ANALYSIS TO SUPPORT YOUR DECISION:\n")

###Strategy A: Universal masking from Day 3
set.seed(42)
strategy_A <- sis_model(adj_matrix, beta, gamma, initial_infected, days, mask_day = 3, mask_efficacy = 0.4)
infections_A <- sum(rowSums(strategy_A$states))

##Strategy B: Isolate top 3 connected (already calculated above)
infections_B <- targeted_infections

##Strategy C: Random isolation of 6 people
adj_matrix_random_isolation <- adj_matrix
set.seed(999)
random_isolated <- sample(1:nrow(adj_matrix), 6)
adj_matrix_random_isolation[random_isolated, ] <- 0
adj_matrix_random_isolation[, random_isolated] <- 0

set.seed(42)
strategy_C <- sis_model(adj_matrix_random_isolation, beta, gamma, initial_infected, days)
infections_C <- sum(rowSums(strategy_C$states))

baseline_total <- sum(rowSums(result1$states))

cat("Baseline (no intervention):", baseline_total, "total infections\n")
cat("Strategy A (Universal masks Day 3):", infections_A, "total infections -", 
    round(((baseline_total - infections_A)/baseline_total)*100, 1), "% reduction\n")
cat("Strategy B (Isolate top 3 connected):", infections_B, "total infections -", 
    round(((baseline_total - infections_B)/baseline_total)*100, 1), "% reduction\n")
cat("Strategy C (Random isolation of 6):", infections_C, "total infections -", 
    round(((baseline_total - infections_C)/baseline_total)*100, 1), "% reduction\n\n")

# Visual comparison of strategies
strategy_comparison <- data.frame(
  Day = rep(0:days, 4),
  Infected = c(rowSums(result1$states), rowSums(strategy_A$states), 
               rowSums(result_targeted$states), rowSums(strategy_C$states)),
  Strategy = rep(c("No Intervention", "Universal Masks (Day 3)", 
                   "Isolate Top 3 Connected", "Random Isolation (6 people)"), 
                 each = days + 1)
)

strategy_plot <- ggplot(strategy_comparison, aes(x = Day, y = Infected, color = Strategy)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  labs(title = "Comparison of Public Health Intervention Strategies",
       subtitle = "Which strategy would you choose and why?",
       x = "Day",
       y = "Number of Infected Individuals") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green", "orange")) +
  theme(legend.position = "bottom")

print(strategy_plot)

cat("DISCUSSION POINTS FOR PARTICPANTS:\n")
cat("1. Which strategy is most effective in reducing total infections?\n")
cat("2. Consider the practical feasibility of each strategy:\n")
cat("   - Cost and resource requirements\n")
cat("   - Community acceptance and compliance\n")
cat("   - Ethical considerations\n")
cat("   - Speed of implementation\n\n")

cat("3. How might the effectiveness change if:\n")
cat("   - The outbreak started with more initial cases?\n")
cat("   - The transmission rate (β) was higher?\n")
cat("   - The recovery rate (γ) was lower?\n\n")

cat("4. What additional information would help you make a better decision?\n")
cat("5. How does network structure influence the success of different strategies?\n\n")

cat("BONUS ANALYSIS: Cost-Effectiveness\n")
cat("Assume costs: Universal masking = $1000, Isolation = $500/person\n")
infections_prevented_A <- baseline_total - infections_A
infections_prevented_B <- baseline_total - infections_B
infections_prevented_C <- baseline_total - infections_C

cost_A <- 1000
cost_B <- 3 * 500  # 3 people isolated
cost_C <- 6 * 500  # 6 people isolated

cost_effectiveness_A <- round(cost_A / infections_prevented_A, 0)
cost_effectiveness_B <- round(cost_B / infections_prevented_B, 0)
cost_effectiveness_C <- round(cost_C / infections_prevented_C, 0)

cat("Cost per infection prevented:\n")
cat("Strategy A: $", cost_effectiveness_A, " per infection prevented\n")
cat("Strategy B: $", cost_effectiveness_B, " per infection prevented\n")
cat("Strategy C: $", cost_effectiveness_C, " per infection prevented\n\n")

cat("=== SUMMARY FOR REFLECTION ===\n")
cat("This exercise demonstrates key epidemiological principles:\n")
cat("• Network structure critically influences disease transmission\n")
cat("• Early intervention is more effective than delayed action\n")
cat("• Targeted interventions can be more efficient than broad measures\n")
cat("• Public health decisions require balancing effectiveness, cost, and feasibility\n")
cat("• Mathematical modeling helps evaluate intervention strategies before implementation\n\n")

##############################################################################
