library(readxl)
library(stringr)

data <- read_excel("./Data/intercal_cetacean_detections.xls") %>%
  filter(Reads > 0)


metadata <- read.csv("./Data/intercal_ALL_metadata_072325.csv") %>%
  select(sample.id, latitude, longitude) %>%
  distinct() %>%
  rename("eDNA Sample ID" = sample.id)


data <- data %>%
  left_join(metadata, by = "eDNA Sample ID") %>%
  filter(BestTaxon %in% c("Delphinus capensis", 
                          "Delphinus delphis",
                          "Delphinus",
                          "Delphinus sp. 1 AN-2013"))


ggplot(data) +
  geom_point(aes(x=longitude, y= latitude, size = Reads, color = BestTaxon)) +
  facet_wrap(~BestTaxon) +
  theme_bw()
