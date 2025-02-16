library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpattern)
library(stringr)
library(ggsignif)
library(glue)
library(gridExtra)

data <-
  read.csv(
    'C:\\Users/Filippo/Desktop/October exp/Evo/PFUDipCoatedBeadsGlass1xvs10xbead_MScThesisTesani2025/PFUDipCoatedBeadsGlass1xvs10xbead_MScThesisTesani2025.csv'
  )

# define column ids and colors
colors <- c("gray85", "gray60")
unique_conditions <- c("1", "10")
color_map <- setNames(colors, unique_conditions)
labels <-
  c("1x bead",
    "10x bead")

# Extract 1x or 10x
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("bead_count", "bead_material"),
    names_pattern = ".*\\.(\\d+)\\.(\\w+).*"
  )

data_wide <- data_long %>%
  pivot_wider(names_from = bead_count,
              values_from = value)

df <- data_wide
list_cols <- names(df)[sapply(df, is.list)]
# Create new dataframe starting with non-list columns
df_clean <- df[, !names(df) %in% list_cols]
# Add mean columns
for (col in list_cols) {
  df_clean[paste0("mean_", col)] <-
    sapply(df[[col]], function(x)
      mean(x, na.rm = TRUE))
}
# Add standard deviation columns
for (col in list_cols) {
  df_clean[paste0("sd_", col)] <-
    sapply(df[[col]], function(x)
      sd(x, na.rm = TRUE))
}

summary_data <- df_clean
summary_data_long <- summary_data %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    values_to = "value",
    names_prefix = "mean_"
  )

# Pivot means
df_means <- summary_data %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "bead_count",
    names_prefix = "mean_",
    values_to = "pfu_count"
  )

# Pivot sds
df_sds <- summary_data %>%
  pivot_longer(
    cols = starts_with("sd_"),
    names_to = "bead_count_sd",
    names_prefix = "sd_",
    values_to = "sd"
  )

# Join by day + condition
summary_data_long <- df_means %>%
  left_join(df_sds,
            by = c("bead_material", "bead_count" = "bead_count_sd")) %>%
  mutate(is_1bead = !grepl("10", bead_count))  # TRUE if it is 1 bead count

# Capitalise bead material names and switch "plastic" for "DaniaPlast"
summary_data_long$bead_material <-
  gsub("plastic", "DaniaPlast", summary_data_long$bead_material)
summary_data_long$bead_material <-
  tools::toTitleCase(summary_data_long$bead_material)

#############################

plot <- ggplot(summary_data_long,
               aes(
                 x = factor(bead_material),
                 y = pfu_count,
                 fill = bead_count,
                 pattern = is_1bead
               )) +
  geom_bar(
    pattern_fill = "gray40",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    position = position_dodge(width = 0.9),
    stat = "identity",
    width = 0.9,
    color = "black",
    linewidth = 0.8
  ) +
    geom_text(
      aes(label = ifelse(
        bead_material == "DaniaPlast", round(pfu_count, 1), ""
      )),
      # Show values only for Plastic
      vjust = -1,
      # Position the text above the bar
      position = position_dodge(width = 0.9),
      size = 5,
      color = "black"
    ) +
  ##### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0, 3.5e+05),
    expand = c(0, 0),
    breaks = seq(0, 3.5e+05, by = 0.5e+05),
    labels = scales::scientific
  ) +
  geom_errorbar(
    aes(
      ymin = ifelse(pfu_count - sd < 0, 0, pfu_count - sd),
      ymax = pfu_count + sd
    ),
    position = position_dodge(width = 0.9),
    width = 0.4,
    linewidth = 0.8
  ) +
  scale_fill_manual(name = "",
                    values = color_map,
                    labels = labels) +
  scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(x = "\nBiofilter type",
       y = "PFU/bead\n",
       fill = "Condition",
       title = "this shows you the couples with significance!!!") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.y = element_text(margin = margin(r = 30)),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "lines"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 22)
  ) +
  guides(fill = guide_legend(nrow = 4,
                             override.aes = list(pattern = c("stripe", "none"))),
         pattern = "none")

# DEBUG: show column names on columns to make sure they're right
# plot <- plot +
#   geom_text(
#     data = summary_data_long,
#     aes(
#       x = factor(bead_material),
#       y = pfu_count,
#       fill = bead_count,
#       label = bead_count
#     ),
#     position = position_dodge(width = 0.9),
#     vjust = 0,
#     # Place just above the bar
#     size = 3,
#     inherit.aes = FALSE
#   )

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
