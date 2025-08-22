
# NOTE requires running detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo.R

mm.volsummary <- mm.data %>%
  group_by(Volume_filt_mL) %>%
  summarize(P = sum(delphinus_all)/length(delphinus_all)) %>%
  mutate(PCorr = P/max(P))

# Create 4-panel plot for site x depth model
ggplot(plot_data_capture_vol, aes(x = volume, y = capture)) +
#  geom_boxplot(data = mm.depthsummary, aes(x=Depth_round, y= PCorr, group = Depth_round)) +
  geom_point(data = mm.data, aes(x=Volume_filt_mL, y= delphinus_all), color = "gray")+
  stat_lineribbon(alpha = 0.25, fill = "#EE7AE9", color = "#EE7AE9", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume (mL)", y = "Probability of Capture") +
  theme_minimal()
