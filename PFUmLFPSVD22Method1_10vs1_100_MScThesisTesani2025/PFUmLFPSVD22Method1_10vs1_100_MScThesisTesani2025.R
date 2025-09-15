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
    'C:\\Users/jasmi/Desktop/October exp/Evo/PFUmLFPSVD22Method1_10vs1_100_MScThesisTesani2025/PFUmLFPSVD22Method1_10vs1_100_MScThesisTesani2025.csv'
  )

colors <- c("gray85", "gray60", "gray35")
labels <-
  c("1:10 phage stock", "1:100 phage stock", "100% phage stock")
color_map = c("100%" = "gray85",
              "1:100" = "gray60",
              "1:10" = "gray35")

# rename data columns
data <- data %>% rename("100%" = X100..phage.stock,
                        "1:10" = X1.10.phage.stock,
                        "1:100" = X1.100.phage.stock)

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
unique_conditions <- df_long$treatment

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
    limits = c(0,  3.5e+09),
    expand = c(0, 0),
    breaks = seq(0,  3.5e+09, by = 0.5e+09),
    labels = scales::scientific
  ) +
  geom_text(
    aes(label = ifelse(
      mean < 0.5e+09, format(mean, scientific = TRUE, digits = 2), ""
    )),
    # Show values only for Plastic
    vjust = -1,
    # Position the text above the bar
    position = position_dodge(width = 0.9),
    size = 5,
    color = "black"
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
  labs(x = "\nPhage stock",
       y = "PFU/mL\n",
       fill = "Condition",
       title = "FPSV-D22") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.y = element_text(margin = margin(r = 30)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
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

# # DEBUG: show column names on columns to make sure they're right
# plot <- plot +
#   geom_text(
#     data = df_long,
#     aes(
#       x = treatment,
#       y = mean + 0.05e+09,
#       label = glue("{treatment}")
#     ),
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
