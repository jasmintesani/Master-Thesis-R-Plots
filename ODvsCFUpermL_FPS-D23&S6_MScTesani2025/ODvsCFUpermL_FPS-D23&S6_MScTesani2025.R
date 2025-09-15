library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpattern)
library(stringr)
library(ggsignif)
library(glue)
library(gridExtra)

split_name <- function(string) {
  # get day number from column and remove it for plotting
  matches <- gsub("^d(\\d)(.*)$", "\\1,\\2", string)
  split_matches <- strsplit(matches, ",")[[1]]
  return(split_matches)
}

data <-
  read.csv(
    'C:\\Users/jasmi/Desktop/October exp/Evo/ODvsCFUpermL_FPS-D23&S6_MScTesani2025/ODvsCFUpermL_FPS-D23&S6_MScTesani2025.csv'
  )

# apply log10 to the values
data$CFU.mL.D23 <- sapply(data$CFU.mL.D23, log10)
data$CFU.ml.S6 <- sapply(data$CFU.ml.S6, log10)

data_wide <- data %>%
  group_by(OD525) %>%
  summarise(
    CFU.mL.D23 = list(CFU.mL.D23),
    CFU.ml.S6 = list(CFU.ml.S6)
  )

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
    names_to = "condition",
    names_prefix = "mean_",
    values_to = "absorbance"
  )

# Pivot sds
df_sds <- summary_data %>%
  pivot_longer(
    cols = starts_with("sd_"),
    names_to = "condition_sd",
    names_prefix = "sd_",
    values_to = "sd"
  )

# Join by day + condition
summary_data_long <- df_means %>%
  left_join(df_sds,
            by = c("OD525", "condition" = "condition_sd")) %>%
  mutate(is_cb = condition == "tyb")

# rename columns to make them look nice
labels <-
  c("FPS-D23",
    "FPS-S6" )

color_map <- c(
  "CFU.mL.D23" = "gray75",
  "CFU.ml.S6" = "gray40"
)



write_xlsx(summary_data_long,
           path = glue(
             "C:\\Users/jasmi/Desktop/October exp/Evo/",
             "ODvsCFUpermL_FPS-D23&S6_MScTesani2025",
             "/ODvsCFUpermL_FPS-D23&S6_MScTesani2025.xlsx"
           ),
           col_names = TRUE)

#############################

plot <- ggplot(summary_data_long,
               aes(
                 x = factor(OD525),
                 y = absorbance,
                 fill = condition,
                 pattern = is_cb
               )) +
  geom_bar_pattern(
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
  ###### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0, 11),
    expand = c(0, 0),
    breaks = seq(0, 11, by = 1),
    #labels = scales::scientific
  ) +
  geom_errorbar(
    aes(ymin = absorbance - sd,
        ymax = absorbance + sd),
    position = position_dodge(width = 0.9),
    width = 0.4,
    linewidth = 0.8
  ) +
  scale_fill_manual(name = "",
                    values = color_map,
                    labels = labels) +
  scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(x = "\nOD at 525 nm",
       y = expression(bold(log[10] * "(CFU/mL)")),
       fill = "Condition",
       title = "") +
  theme_minimal() +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.y = element_text(margin = margin(r = 20)),
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
  guides(fill = guide_legend(nrow = 8,
                             override.aes = list(
                               pattern = c("none",
                                           "none")
                             )),
         pattern = "none")

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
