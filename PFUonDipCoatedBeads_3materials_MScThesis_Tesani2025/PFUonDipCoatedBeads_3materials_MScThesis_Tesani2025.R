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
    'C:\\Users/jasmi/Desktop/October exp/Evo/PFUonDipCoatedBeads_3materials_MScThesis_Tesani2025/PFUonDipCoatedBeads_3materials_MScThesis_Tesani2025.csv'
  )

# apply log10 to data means
data$PFU.DaniaPlast.bead <- log10(data$PFU.DaniaPlast.bead)
data$PFU.Glass.bead <- log10(data$PFU.Glass.bead)
data$PFU.BioChip.bead <- log10(data$PFU.BioChip.bead)

colors <- c("gray85", "gray60", "gray35")
labels <- c("BioChip beads", "DaniaPlast beads", "Glass beads")

# Calculate column means
means <- sapply(data, function(x)
  mean(x, na.rm = TRUE))
means_df <- data.frame(t(means))  # Transpose to make a row

colnames(means_df) <- paste0("mean_", colnames(data))
summary_data_long <- means_df

# Calculate column sds
sds <- sapply(data, function(x)
  sd(x, na.rm = TRUE))
sds_df <- data.frame(t(sds))  # Transpose to make a row

# apply log10 to data sds
sds_df$PFU.DaniaPlast.bead <- log10(sds_df$PFU.DaniaPlast.bead)
sds_df$PFU.Glass.bead <- log10(sds_df$PFU.Glass.bead)
sds_df$PFU.BioChip.bead <- log10(sds_df$PFU.BioChip.bead)

colnames(sds_df) <- paste0("sd_", colnames(data))

summary_data_long <- merge(summary_data_long, sds_df)
df <- summary_data_long

df_long <- df %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "treatment",
               values_to = "mean") %>%
  pivot_longer(cols = starts_with("sd"),
               names_to = "treatment_sd",
               values_to = "sd") %>%
  mutate(
    treatment = sub("mean_", "", treatment),
    # Remove 'mean_' from treatment names
    treatment_sd = sub("sd_", "", treatment_sd)
  ) %>%
  filter(treatment == treatment_sd)  # Keep only matching treatment and sd

df_long$treatment <-
  str_split(df_long$treatment, "\\.", simplify = TRUE)[, 2]

# fix material names
df_long$treatment <-
  sapply(df_long$treatment, function(x)
    paste(x, " beads"))

df_long$treatment[grepl("dania", df_long$treatment, ignore.case = TRUE)] <-
  "1DaniaPlast beads"
df_long$treatment[grepl("glass", df_long$treatment, ignore.case = TRUE)] <-
  "2Glass beads"
df_long$treatment[grepl("biochip", df_long$treatment, ignore.case = TRUE)] <-
  "3BioChip beads"

unique_conditions <- df_long$treatment
color_map <- setNames(colors, unique_conditions)

#############################

plot <- ggplot(df_long, aes(x = treatment,
                            y = mean,
                            fill = treatment)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat = "identity",
    width = 0.9,
    color = "black",
    linewidth = 0.8
  ) +
  scale_x_discrete(
    labels = function(x)
      gsub("\\d+", "", x)
  ) +
  ##### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0,  8),
    expand = c(0, 0),
    breaks = seq(0,  8, by = 1),
    # labels = scales::scientific
  ) +
  geom_errorbar(
    aes(ymin = ifelse(mean - sd < 0, 0, mean - sd),
        ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.4,
    linewidth = 0.8
  ) +
  scale_fill_manual(name = "",
                    values = color_map,
                    labels = labels) +
  # scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(
    x = "\nBiofilter type",
    y = expression(bold(log[10] * "(PFU/bead)")),
    fill = "Condition",
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
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
                             override.aes = list(pattern = c("stripe"))),
         pattern = "none")

star_shift_factors <- c("cb" = -0.3, "cbw" = -0.1)
star_data <- star_data %>%
  mutate(shift = star_shift_factors[col1])

# plot <- plot +
#   geom_text(
#     data = star_data,
#     aes(
#       x = factor(day),
#       y = y_star,
#       label = glue("* {col1}-{col2}")
#     ),
#     vjust = 0,   # Place just above the bar
#     size = 3,
#     inherit.aes = FALSE,
#     nudge_x = star_data$shift
#   )

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
