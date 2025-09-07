###### The following exercise ####
##Network Analysis and SIT Epidemic Model##

####Project Description####

##This project introduces the intersection of network analysis and epidemic modeling through a hands-on exploration of disease transmission dynamics. Using a dataset of participants from various countries attending a Malaria Modelling modular program in Nairobi, Kenya, participants learn to construct social networks based on participants interaction patterns and simulate disease spread using a Susceptible-Infected-Treatment (SIT) model with separate symptomatic and asymptomatic infected compartments.

##The project aims to demonstrate how network structure influences epidemic outcomes by modeling connections between individuals with a 10% probability, creating a more realistic interaction pattern. The aim is to explore key epidemiological concepts including transmission rates (β), treatment rates (τ), recovery rates (γ), and intervention strategies like wearing masks and identification of superspreaders.

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


# Network metrics
cat("d) Network Metrics:\n")

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
plot_network_degree <- function(g, data) {
  ##Calculate degree centrality
  degrees <- betweenness(g)
  
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
plot_network_betweeness <- function(g, data) {
  ##Calculate degree centrality
  degrees <- betweenness(g)
  
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
       main = "Network with Node Sizes Proportional to Betweenness Centrality")
  
  legend("topright", legend = countries, fill = country_colors, cex = 0.8, title = "Country")
}

plot_network_betweeness(g, data)

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

## a). Identify top 6 most connected individuals
top6_connected <- sort(degrees, decreasing = TRUE)[1:6]
top6_names <- names(top6_connected)

cat("a) Top 6 most connected individuals:\n")
for(i in 1:6) {
  cat("   ", i, ".", top6_names[i], "-", top6_connected[i], "connections\n")
}

##Exercise 3.1. - Who are the top 6 most connected individuals based on betweenness? Do they differ from those based on degree centrality? Comment on whether they are different or the same - why?

