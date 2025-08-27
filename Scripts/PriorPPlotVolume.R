# Script to plot prior predictive checks for volume
# NOTE requires running detection_decomp_MMs_real_data_with_miFish_and_NCOG_depthAndSiteCombo.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)

mm.volsummary <- mm.data %>%
  group_by(Volume_filt_mL) %>%
  summarize(P = sum(delphinus_all)/length(delphinus_all)) %>%
  mutate(PCorr = P/max(P))

vol_vals_centered <- seq(-3000, 3000, by = 250)

df <- expand_grid(Cap_prob_site = seq(0, 1, by = 0.1),
              Vol_Centered = vol_vals_centered, 
              Method = 1,
              Betas_vol = seq(-.005, .005, by = 0.0005)) %>%
  mutate(Pred = inv.logit(Cap_prob_site + Vol_Centered*Betas_vol)) %>%
  mutate(Group = paste0(Betas_vol, "_", Cap_prob_site))

ggplot(filter(df)) +
  geom_line(aes(x= Vol_Centered+3000, y = Pred, group = Group), col = "gray") +
  geom_point(data = mm.data, aes(x=Volume_filt_mL, y= delphinus_all)) +
  xlab("Volume (mL)") +
  ylab("Predicted Probability of Capture") +
  theme_minimal()

ggplot(filter(df)) +
  geom_line(aes(x= Vol_Centered+3000, y = Pred, group = Group), col = "gray") +
  geom_boxplot(data = mm.volsummary, aes(x=Volume_filt_mL, y= PCorr, group = Volume_filt_mL)) +
  scale_y_continuous(name = "Predicted Probability of Capture",
                     sec.axis = sec_axis(~.*1, name="Observed Proportion Captured")) +
  xlab("Volume (mL)") +
  theme_minimal()
