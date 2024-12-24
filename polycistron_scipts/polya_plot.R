

library(ggplot2)
library(viridis)




prefix <- "sample"
input_file <- paste0(prefix, "_transcript_VS_polya_length.txt")

polya_data <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(polya_data) <- c("PolyA_Size", "Transcript_Size")



ggplot(polya_data, aes(x = PolyA_Size, y = Transcript_Size)) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 10)) +
  ylim(0, 2500) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) +
  scale_fill_viridis(option = "magma") +
  labs(
    title = "PolyA Size vs Transcript Size",
    x = "PolyA Size",
    y = "Transcript Size",
    fill = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black", size = 16, face = "bold"),
    axis.title = element_text(color = "black", size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )