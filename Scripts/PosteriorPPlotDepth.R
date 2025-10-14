
# Script to plot posterior predictive checks for depth effect
# NOTE requires running detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo_v2.R
# source("./detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo_v2.R)

mm.depthsummary <- mm.data %>%
  mutate(Depth_round = round(Depth_m, digits = -2)) %>%
  # filtered these out bc there were only a couple of obs
#  filter(Depth_round != 300, Depth_round != 400) %>% 
  group_by(Depth_round) %>%
  summarize(P = sum(delphinidae_all)/length(delphinidae_all), N = length(delphinidae_all)) %>%
  mutate(PPos = P/max(P))

# change this to 25 for the sake of pretty plotting
mm.depthsummary$Depth_round[1] <- 25

# Create plot for site x depth model with observations as points
p1 <- ggplot(plot_data_occupancy_depth, aes(x = depth, y = occupancy)) +
  geom_point(data = mm.data, aes(x=Depth_m, y= delphinidae_all), color = "gray")+
 
  stat_lineribbon(alpha = 0.25, fill = "#4CAF50", color = "#2E7D32", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth (m)", y = "Predicted Probability of Occupancy") +
  theme_minimal()

# Create plot for site x depth model with observations summarized in bins
p2 <- ggplot(plot_data_occupancy_depth, aes(x = depth, y = occupancy)) +
  scale_y_continuous(name = "Predicted Probability of Occupancy") +
  # scale_y_continuous(name = "Predicted Probability of Occupancy",
  #                     sec.axis = sec_axis(~.*max(mm.depthsummary$P), name="Observed Proportion Occupied")) +
  stat_lineribbon(alpha = 0.25, fill = "#4CAF50", color = "#2E7D32", 
                  .width = c(0.25, 0.5, 0.75)) +
  geom_point(data = mm.depthsummary, 
             aes(x=Depth_round, y= PPos, group = Depth_round, size = N), 
             color = "gray") +
  labs(x = "Depth (m)") +
  theme_minimal()

ggsave(plot = p2, file = "./Figures/DepthPosteriorPredCheck.png", 
       width = 8, height = 5, units = "in")
