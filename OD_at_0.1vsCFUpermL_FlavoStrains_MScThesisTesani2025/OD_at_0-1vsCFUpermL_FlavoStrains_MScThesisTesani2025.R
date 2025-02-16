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
    'C:\\Users/Filippo/Desktop/October exp/Evo/OD_at_0.1vsCFUpermL_FlavoStrains_MScThesisTesani2025/OD_at_0.1vsCFUpermL_FlavoStrains_MScThesisTesani2025.csv'
  )

colors <-
  c("gray85",
    "gray75",
    "gray65",
    "gray55",
    "gray45",
    "gray35",
    "gray25")
unique_conditions <- df_long$treatment
color_map <- c(
  "d23" = "gray85",
  "d27" = "gray75",
  "d42" = "gray65",
  "d45" = "gray55",
  "d49" = "gray45",
  "s6" = "gray35",
  "X960106.1.1" = "gray25"
)

# Calculate column means
means <- sapply(data, function(x)
  mean(x, na.rm = TRUE))
means_df <- data.frame(t(means))  # Transpose to make a row
colnames(means_df) <- paste0("mean_", colnames(data))
summary_data_long <- means_df

# apply log10
for (col in names(summary_data_long)) {
  summary_data_long[col] <- sapply(summary_data_long[col], log10)
}

# Calculate column sds
sds <- sapply(data, function(x) {
      sd(log10(x), na.rm = TRUE)
    }
  )
sds_df <- data.frame(t(sds))  # Transpose to make a row
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

df_long$treatment <- sub("^d23$", "FPS-D23", df_long$treatment)
df_long$treatment <- sub("^d27$", "FPS-D27", df_long$treatment)
df_long$treatment <- sub("^d42$", "FPS-D42", df_long$treatment)
df_long$treatment <- sub("^d45$", "FPS-D45", df_long$treatment)
df_long$treatment <- sub("^d49$", "FPS-D49", df_long$treatment)
df_long$treatment <- sub("^s6$", "FPS-S6", df_long$treatment)
df_long$treatment <- sub("^X960106.1.1$", "960106-1/1", df_long$treatment)

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
  ##### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0, 9.1),
    expand = c(0, 0),
    breaks = seq(0, 9.1, by = 1)
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
  labs(
    x = "F. psychrophilum strains",
    y = expression(bold(log[10] * "(CFU/mL)")),
    fill = "Condition",
    title = expression(bold("Bacterial cultures at an" ~ OD[525 * " nm"] ~ "of 0.1"))
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 30)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18, colour = "transparent"),
    legend.key.size = unit(0, "lines"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    axis.text = element_text(size = 19)
  )

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
