# Script to plot prior predictive checks for depth

library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)

mm.depthsummary <- mm.data %>%
  mutate(Depth_round = round(Depth_m, digits = -2)) %>%
  # filtered these out bc there were only a couple of obs
  filter(Depth_round != 300, Depth_round != 400) %>% 
  group_by(Depth_round) %>%
  summarize(P = sum(delphinus_all)/length(delphinus_all)) %>%
  mutate(PPos = P/max(P))

beta_vals <- rnorm(100, 0, 1)
intercept_vals <- (rnorm(100, 0, 1))
depth_vals <- seq(0, 500, by = 10)
depth_vals_centered <- depth_vals - 132.9

df <- data.frame("Betas" = beta_vals,
           "Intercepts" = intercept_vals,
           "Group" = 1:100) %>%
  expand_grid(Depth_Vals = depth_vals)

df$Depth_Centered <- df$Depth_Vals - 132.9

df$Pred <- inv.logit(df$Intercepts + 
                       df$Betas * df$Depth_Centered)

ggplot(df) +
  geom_line(aes(x= Depth_Vals, y = Pred, group = Group), col = "gray") +
  geom_boxplot(data = mm.depthsummary, aes(x=Depth_round, y= PPos, group = Depth_round)) +
  scale_y_continuous(name = "Predicted Probability of Occupancy",
                     sec.axis = sec_axis(~.*max(mm.depthsummary$P), name="Observed Proportion Occupied")) +
  xlab("Depth (m)") +
  theme_minimal()
