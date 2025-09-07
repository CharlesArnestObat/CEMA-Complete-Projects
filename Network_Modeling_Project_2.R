###### The following exercise ####
##Network Analysis and SEIR Epidemic Model##

####Project Description####

##This project introduces the intersection of network analysis and epidemic modeling through a hands-on exploration of disease transmission dynamics. Using a dataset of participants from various countries attending a Malaria Modelling modular program in Nairobi, Kenya, participants learn to construct social networks based on participants interaction patterns and simulate disease spread using a Susceptible-Exposed-Infected-Recovered (SEIR) model.

##The project aims to demonstrate how network structure influences epidemic outcomes by modeling connections between individuals with a 10% probability, creating a more realistic interaction pattern. The aim is to explore key epidemiological concepts including transmission rates (β), incubation rates (σ), recovery rates (γ), and intervention strategies like wearing masks and identification of superspreaders.

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
#Modified: Everyone has a 10% chance of connecting to everybody else
create_network <- function(data, connection_prob = 0.1) {
  n <- nrow(data)
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- data$Participant
  
  set.seed(123)  # For reproducibility
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        # 10% chance of connection for any pair of individuals
        if(runif(1) < connection_prob) {
          adj_matrix[i,j] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}
#increasing connection/mixing/probability to 70%-connection_prob = 0.7
create_network2 <- function(data, connection_prob = 0.7) {
  n <- nrow(data)
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- data$Participant
  
  set.seed(123)  # For reproducibility
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        # 10% chance of connection for any pair of individuals
        if(runif(1) < connection_prob) {
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
degree(g)
betweenness(g)
closeness(g)


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

####Part 2 - SEIR Epidemic Modelling ####
seir_model <- function(adj_matrix, beta, sigma, gamma, initial_exposed, days, superspreaders = NULL, mask_day = NULL, mask_efficacy = 0.3) {
  n <- nrow(adj_matrix)
  participants <- rownames(adj_matrix)
  
  # Initialize states: 0 = Susceptible, 1 = Exposed, 2 = Infected, 3 = Recovered
  states <- array(0, dim = c(days + 1, n, 4))
  dimnames(states) <- list(NULL, participants, c("S", "E", "I", "R"))
  
  # Set initial states - everyone starts susceptible except initially exposed
  states[1, , "S"] <- 1
  exposed_idx <- which(participants %in% initial_exposed)
  states[1, exposed_idx, "S"] <- 0
  states[1, exposed_idx, "E"] <- 1
  
  # Track daily new infections and exposures for incidence calculation
  daily_new_exposures <- rep(0, days + 1)
  daily_new_infections <- rep(0, days + 1)
  
  # Superspreaders - who is a superspreader?
  ss_multiplier <- rep(1, n)
  if(!is.null(superspreaders)) {
    ss_idx <- which(participants %in% superspreaders)
    ss_multiplier[ss_idx] <- 5  # 5x transmission rate
  }
  
  for(day in 1:days) {
    # Current state compartments
    S_current <- states[day, , "S"]
    E_current <- states[day, , "E"]
    I_current <- states[day, , "I"]
    R_current <- states[day, , "R"]
    
    # Initialize next day's states
    S_next <- S_current
    E_next <- E_current
    I_next <- I_current
    R_next <- R_current
    
    new_exposures <- 0
    new_infections <- 0
    
    # Apply Intervention - masking effect if applicable
    current_beta <- beta
    if(!is.null(mask_day) && day >= mask_day) {
      current_beta <- beta * (1 - mask_efficacy)
    }
    
    # SEIR model equations implemented stochastically
    for(i in 1:n) {
      # S -> E transition (Susceptible to Exposed)
      if(S_current[i] == 1) {
        # Calculate force of infection from infected neighbors
        infected_neighbors <- which(I_current == 1 & adj_matrix[i, ] == 1)
        if(length(infected_neighbors) > 0) {
          # Calculate total force of infection
          total_force <- 0
          for(neighbor in infected_neighbors) {
            total_force <- total_force + current_beta * ss_multiplier[neighbor]
          }
          
          # Probability of exposure (1 - probability of avoiding all exposures)
          prob_exposure <- 1 - exp(-total_force)
          if(runif(1) < prob_exposure) {
            S_next[i] <- 0
            E_next[i] <- 1
            new_exposures <- new_exposures + 1
          }
        }
      }
      
      # E -> I transition (Exposed to Infected)
      else if(E_current[i] == 1) {
        # Probability of becoming infectious (incubation period ends)
        if(runif(1) < sigma) {
          E_next[i] <- 0
          I_next[i] <- 1
          new_infections <- new_infections + 1
        }
      }
      
      # I -> R transition (Infected to Recovered)
      else if(I_current[i] == 1) {
        # Probability of recovery
        if(runif(1) < gamma) {
          I_next[i] <- 0
          R_next[i] <- 1
        }
      }
    }
    
    # Update states for next day
    states[day + 1, , "S"] <- S_next
    states[day + 1, , "E"] <- E_next
    states[day + 1, , "I"] <- I_next
    states[day + 1, , "R"] <- R_next
    
    daily_new_exposures[day + 1] <- new_exposures
    daily_new_infections[day + 1] <- new_infections
  }
  
  return(list(
    states = states, 
    daily_new_exposures = daily_new_exposures,
    daily_new_infections = daily_new_infections
  ))
}
####Part 2 - SEIR Epidemic Modelling #### mask_efficacy reduced to 10% -mask_efficacy = 0.1
seir_model2 <- function(adj_matrix, beta, sigma, gamma, initial_exposed, days, superspreaders = NULL, mask_day = NULL, mask_efficacy = 0.3) {
  n <- nrow(adj_matrix)
  participants <- rownames(adj_matrix)
  
  # Initialize states: 0 = Susceptible, 1 = Exposed, 2 = Infected, 3 = Recovered
  states <- array(0, dim = c(days + 1, n, 4))
  dimnames(states) <- list(NULL, participants, c("S", "E", "I", "R"))
  
  # Set initial states - everyone starts susceptible except initially exposed
  states[1, , "S"] <- 1
  exposed_idx <- which(participants %in% initial_exposed)
  states[1, exposed_idx, "S"] <- 0
  states[1, exposed_idx, "E"] <- 1
  
  # Track daily new infections and exposures for incidence calculation
  daily_new_exposures <- rep(0, days + 1)
  daily_new_infections <- rep(0, days + 1)
  
  # Superspreaders - who is a superspreader?
  ss_multiplier <- rep(1, n)
  if(!is.null(superspreaders)) {
    ss_idx <- which(participants %in% superspreaders)
    ss_multiplier[ss_idx] <- 5  # 5x transmission rate
  }
  
  for(day in 1:days) {
    # Current state compartments
    S_current <- states[day, , "S"]
    E_current <- states[day, , "E"]
    I_current <- states[day, , "I"]
    R_current <- states[day, , "R"]
    
    # Initialize next day's states
    S_next <- S_current
    E_next <- E_current
    I_next <- I_current
    R_next <- R_current
    
    new_exposures <- 0
    new_infections <- 0
    
    # Apply Intervention - masking effect if applicable
    current_beta <- beta
    if(!is.null(mask_day) && day >= mask_day) {
      current_beta <- beta * (1 - mask_efficacy)
    }
    
    # SEIR model equations implemented stochastically
    for(i in 1:n) {
      # S -> E transition (Susceptible to Exposed)
      if(S_current[i] == 1) {
        # Calculate force of infection from infected neighbors
        infected_neighbors <- which(I_current == 1 & adj_matrix[i, ] == 1)
        if(length(infected_neighbors) > 0) {
          # Calculate total force of infection
          total_force <- 0
          for(neighbor in infected_neighbors) {
            total_force <- total_force + current_beta * ss_multiplier[neighbor]
          }
          
          # Probability of exposure (1 - probability of avoiding all exposures)
          prob_exposure <- 1 - exp(-total_force)
          if(runif(1) < prob_exposure) {
            S_next[i] <- 0
            E_next[i] <- 1
            new_exposures <- new_exposures + 1
          }
        }
      }
      
      # E -> I transition (Exposed to Infected)
      else if(E_current[i] == 1) {
        # Probability of becoming infectious (incubation period ends)
        if(runif(1) < sigma) {
          E_next[i] <- 0
          I_next[i] <- 1
          new_infections <- new_infections + 1
        }
      }
      
      # I -> R transition (Infected to Recovered)
      else if(I_current[i] == 1) {
        # Probability of recovery
        if(runif(1) < gamma) {
          I_next[i] <- 0
          R_next[i] <- 1
        }
      }
    }
    
    # Update states for next day
    states[day + 1, , "S"] <- S_next
    states[day + 1, , "E"] <- E_next
    states[day + 1, , "I"] <- I_next
    states[day + 1, , "R"] <- R_next
    
    daily_new_exposures[day + 1] <- new_exposures
    daily_new_infections[day + 1] <- new_infections
  }
  
  return(list(
    states = states, 
    daily_new_exposures = daily_new_exposures,
    daily_new_infections = daily_new_infections
  ))
}

###3 - SIMULATIONS

# Set parameters
beta <- 0.3      # Transmission rate
sigma <- 0.2     # Rate of progression from exposed to infected (1/incubation period)
gamma <- 0.1     # Recovery rate
days <- 30

##Find participant indices
zeinabu_idx <- which(data$Participant == "ZEN")
agnes_idx <- which(data$Participant == "AGN")
initial_exposed <- c("ZEN", "AGN")

##Scenario 1: Basic SEIR model
set.seed(42)
result1 <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days) #Explore

##Scenario 2: With masks from day 3
set.seed(42)
result2 <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, mask_day = 3)
##Scenario 2: With masks from day 15
set.seed(42)
result2b <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, mask_day = 15)
top6_connected <- sort(degree(g), decreasing = TRUE)[1:6]

##Scenario 3: With superspreaders
superspreaders <- c("DGB", "BAS", "WAN", "PUR") #Are these the superspreaders? Who is a superspreader? What criteria do we use to identify a superspreader? Hint: Network metrics. 
set.seed(42)
result3 <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, superspreaders = superspreaders)

####Part 4. ANALYSIS AND RESULTS ####

##Epidemic Results

# Helper function to get compartment totals
get_compartment_totals <- function(states, day, compartment) {
  sum(states[day + 1, , compartment])  # day + 1 because array is 1-indexed
}

# a) Number of infected on days 2, 4, 8
infected_day2 <- get_compartment_totals(result1$states, 2, "I")
infected_day4 <- get_compartment_totals(result1$states, 4, "I")
infected_day5 <- get_compartment_totals(result1$states, 5, "I")
infected_day8 <- get_compartment_totals(result1$states, 8, "I")
infected_day20 <- get_compartment_totals(result1$states, 20, "I")
cat("a) Number of infected individuals:\n")
cat("   Day 2:", infected_day2, "people\n")
cat("   Day 4:", infected_day4, "people\n")
cat("   Day 8:", infected_day8, "people\n\n")
cat("   Day 20:", infected_day20, "people\n\n")
# a) Number of infected on days 2, 4, 8, 20 when masks introduced on day 15
infected_day2 <- get_compartment_totals(result2b$states, 2, "I")
infected_day4 <- get_compartment_totals(result2b$states, 4, "I")
infected_day5 <- get_compartment_totals(result2b$states, 5, "I")
infected_day8 <- get_compartment_totals(result2b$states, 8, "I")
infected_day20 <- get_compartment_totals(result2b$states, 20, "I")
cat("a) Number of infected individuals:\n")
cat("   Day 2:", infected_day2, "people\n")
cat("   Day 4:", infected_day4, "people\n")
cat("   Day 8:", infected_day8, "people\n\n")
cat("   Day 20:", infected_day20, "people\n\n")
# b) Incidence at day 8
incidence_day8 <- result1$daily_new_infections[5]
cat("b) Incidence on Day 5:", incidence_day8, "new infections\n\n")

# c) Infections averted with masks
infected_day5_masks <- get_compartment_totals(result2$states, 5, "I")
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
infected_normal <- get_compartment_totals(result1$states, 6, "I")
infected_superspreaders <- get_compartment_totals(result3$states, 6, "I")
cat("       Without superspreaders: ", infected_normal, " infected\n")
cat("       With superspreaders: ", infected_superspreaders, " infected\n")
cat("       Additional infections from superspreaders: ", infected_superspreaders - infected_normal, "\n")

## iii) Effect of interaction rate changes
cat("   iii) Effect of changing interaction rates:\n")
betas <- c(0.1, 0.2, 0.3, 0.4, 0.5)
final_infected <- numeric(length(betas))

for(i in seq_along(betas)) {
  set.seed(42)
  temp_result <- seir_model(adj_matrix, betas[i], sigma, gamma, initial_exposed, days)
  final_infected[i] <- get_compartment_totals(temp_result$states, 6, "I")
}

beta_results <- data.frame(Beta = betas, Final_Infected = final_infected)
print(beta_results)

####Part 5. VISUALIZATION OF EPIDEMIC PROGRESSION ####
plot_seir_epidemic <- function(states) {
  days <- dim(states)[1] - 1
  
  # Extract compartment counts for each day
  S_counts <- rowSums(states[, , "S"])
  E_counts <- rowSums(states[, , "E"])
  I_counts <- rowSums(states[, , "I"])
  R_counts <- rowSums(states[, , "R"])
  
  df <- data.frame(
    Day = rep(0:days, 4),
    Count = c(S_counts, E_counts, I_counts, R_counts),
    Compartment = rep(c("Susceptible", "Exposed", "Infected", "Recovered"), each = days + 1)
  )
  
  ggplot(df, aes(x = Day, y = Count, color = Compartment)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.5) +
    labs(title = "SEIR Epidemic Progression",
         x = "Day",
         y = "Number of Individuals") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "orange", "red", "green")) +
    scale_x_continuous(breaks = 0:days)
}

##Plot epidemic progression
print(plot_seir_epidemic(result1$states))

##Comparison plot for infected individuals only
compare_infected_df <- data.frame(
  Day = rep(0:days, 3),
  Infected = c(rowSums(result1$states[, , "I"]), 
               rowSums(result2$states[, , "I"]), 
               rowSums(result3$states[, , "I"])),
  Scenario = rep(c("Basic", "With Masks", "With Superspreaders"), each = days + 1)
)

comparison_plot <- ggplot(compare_infected_df, aes(x = Day, y = Infected, color = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Comparison of SEIR Epidemic Scenarios (Infected)",
       x = "Day",
       y = "Number of Infected Individuals") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:days) +
  scale_color_manual(values = c("red", "blue", "orange"))

print(comparison_plot)

##Enhanced infections averted visualization as a shaded area
plot_infections_averted <- function(result_no_control, result_with_control, intervention_name = "Intervention") {
  days <- dim(result_no_control$states)[1] - 1
  
  ##Calculate infections for both scenarios
  no_control <- rowSums(result_no_control$states[, , "I"])
  with_control <- rowSums(result_with_control$states[, , "I"])
  
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

## a) & b). Run SEIR model with different β values and plot
betas <- seq(0.01, 0.6, by = 0.02)
final_infected_beta <- numeric(length(betas))

for(i in seq_along(betas)) {
  set.seed(42)
  temp_result <- seir_model(adj_matrix, betas[i], sigma, gamma, initial_exposed, days)
  final_infected_beta[i] <- get_compartment_totals(temp_result$states, days, "I")
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
  temp_result <- seir_model(adj_matrix, 0.4, sigma, gammas[i], initial_exposed, days)
  final_infected_gamma[i] <- get_compartment_totals(temp_result$states, days, "I")
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
R0_values <- betas / gamma  # Basic reproduction number for SEIR model
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
  result_mask <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, 
                            mask_day = intervention_days[i], mask_efficacy = mask_efficacy)
  
  # Calculate total infections over entire period
  total_no_mask <- sum(result1$states[, , "I"])
  total_with_mask <- sum(result_mask$states[, , "I"])
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
result_ss_connected <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, superspreaders = top6_names)

# Scenario 2: Random 6 people as superspreaders
set.seed(123)
random_ss <- sample(data$Participant, 6)
set.seed(42)
result_ss_random <- seir_model(adj_matrix, beta, sigma, gamma, initial_exposed, days, superspreaders = random_ss)

## d). Calculate superspreader effects
baseline_infections <- sum(result1$states[, , "I"])
connected_ss_infections <- sum(result_ss_connected$states[, , "I"])
random_ss_infections <- sum(result_ss_random$states[, , "I"])

connected_effect <- round(((connected_ss_infections - baseline_infections) / baseline_infections) * 100, 1)
random_effect <- round(((random_ss_infections - baseline_infections) / baseline_infections) * 100, 1)

cat("\nb-d) Superspreader Effect Analysis:\n")
cat("   Baseline (no superspreaders):", baseline_infections, "total infections\n")
cat("Most connected as superspreaders: ", connected_ss_infections, " total infections (", connected_effect, "% increase)\n")

