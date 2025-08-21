library(tidyverse)
smallpox_table=read.csv("smallpox_table.csv")
smallpox_table <- smallpox_table %>%
  select(age, alive_P, alive_Pinoculated, immune_R) %>%
  summarize(
    life_exp_inoculated = sum(alive_Pinoculated) / 1300,
    life_exp            = sum(alive_P) / 1300
  )
smallpox_table
# Convert to long format for plotting
smallpox_table_long <- smallpox_table%>%
  pivot_longer(cols = everything(),
               names_to = "group",
               values_to = "life_expectancy")
smallpox_table_long

# Simple bar plot
ggplot(smallpox_table_long, aes(x = group, y = life_expectancy, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(life_expectancy, 2)), vjust = -0.5) +
  labs(
    title = "Life expectancy (years)",
    x = NULL, y = "Life expectancy"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

