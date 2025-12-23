library(dplyr)
library(tidyr)

# Import the CSV file
mm.data <- read.csv("../Data/intercal_cetacean_detections.csv")

mm.data.dl <- mm.data %>% filter(Marker == "d-loop")

mm.dl.wide <- mm.data.dl %>% select(!c(eDNA.Sample.ID, Marker)) %>%
  pivot_wider(names_from = Sequence.Sample.ID, values_from = Reads) %>%
  replace(is.na(.), 0)

mm.data.mfu <- mm.data %>% filter(Marker == "MiFish")

mm.mfu.wide <- mm.data.mfu %>% select(!c(eDNA.Sample.ID, Marker)) %>%
  pivot_wider(names_from = Sequence.Sample.ID, values_from = Reads) %>%
  replace(is.na(.), 0)

write.csv(mm.dl.wide, "intercal_DL_cetaceans_detections_wide.csv")

write.csv(mm.mfu.wide, "intercal_MFU_cetaceans_detections_wide.csv")
