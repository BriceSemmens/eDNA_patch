
# NOTE requires running detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo.R

mm.volsummary <- mm.data %>%
  group_by(Vol_round = round(Volume_filt_mL, digits = -3)) %>%
  summarize(P = sum(delphinus_all)/length(delphinus_all), N = length(delphinus_all)) %>%
  mutate(PCorr = P/max(P))

ggplot(plot_data_capture_vol, aes(x = volume, y = capture)) +
  
  scale_y_continuous(name = "Predicted Probability of Capture",
                     sec.axis = sec_axis(~.*max(mm.volsummary$PCorr), name="Observed Proportion Captured")) +
  #geom_point(data = mm.data, aes(x=Volume_filt_mL, y= delphinus_all), color = "gray")+
  stat_lineribbon(alpha = 0.25, fill = "#EE7AE9", color = "#EE7AE9", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume (mL)", y = "Probability of Capture") +
  geom_point(data = mm.volsummary, aes(x=Vol_round, y= PCorr, group = Vol_round, size = N), color = "gray") +
  theme_minimal()

ggsave(plot = last_plot(), file = "./Figures/VolumePosteriorPredCheck.png", 
       width = 5, height = 4, units = "in")
