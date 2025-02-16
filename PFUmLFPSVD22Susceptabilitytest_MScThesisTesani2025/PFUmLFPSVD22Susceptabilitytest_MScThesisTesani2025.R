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
  matches <- gsub("(D[0-9]{2})(.*)", "\\1,\\2", string)
  split_matches <- strsplit(matches, ",")[[1]]
  return(split_matches)
}

data <-
  read.csv(
    'C:\\Users/Filippo/Desktop/October exp/Evo/PFUmLFPSVD22Susceptabilitytest_MScThesisTesani2025/PFUmLFPSVD22Susceptabilitytest_MScThesisTesani2025.csv'
  )

# clean up column names to easily identify groups and variables later
colnames(data) <-
  gsub("^fpsv\\.", "", colnames(data), ignore.case = TRUE)
colnames(data) <-
  gsub("(against|on)\\.", "", colnames(data), ignore.case = TRUE)
colnames(data) <-
  gsub("fps\\.", "", colnames(data), ignore.case = TRUE)
colnames(data) <-
  gsub("\\.", "", colnames(data), ignore.case = TRUE)

# convert to logarithmic data
data <- data %>%
  mutate(across(where(is.numeric), log10))

# filter data for the rest of the code
# define column ids and colors
unique_conditions <-
  c("S11A",
    "D23",
    "D27",
    "D42",
    "D45",
    "D49",
    "S6", )
color_map <- c(
  "D23" = "gray30",
  "D27" = "gray40",
  "D42" = "gray50",
  "D45" = "gray60",
  "D49" = "gray70",
  "S11A" = "gray80",
  "S6" = "gray90"
)
# rename columns to make them look nice
labels <-
  c("FPS-D23",
    "FPS-D27",
    "FPS-D42",
    "FPS-D45",
    "FPS-D49",
    "FPS-S11A",
    "FPS-S6")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("variable", "day_suffix"),
    names_pattern = "(D[0-9]{2})(.*)"
  )
data_wide <- data_long %>%
  pivot_wider(names_from = day_suffix,
              values_from = value)

df <- data_wide %>% rename(phage_type = variable)

list_cols <- names(df)[sapply(df, is.list)]
# Create new dataframe starting with non-list columns
df_clean <- df[, !names(df) %in% list_cols]
# Add mean columns
for (col in list_cols) {
  df_clean[paste0("mean_", col)] <-
    sapply(df[[col]], function(x)
      mean(x, na.rm = TRUE))
}

# collect p-values in a list
wilcox_p_values <- list()

for (j in seq(1, length(colnames(data)), by = 7)) {
  mb_index <- j
  pf43_index <- j + 9
  
  group_to_compare <- colnames(data)[c(mb_index, pf43_index)]
  # calculate the p-value between the 2 cols extracted
  first <- group_to_compare[1]
  second <- group_to_compare[2]
  wilcox_result = wilcox.test(data[[first]], data[[second]])
  
  first <- split_name(first)
  day <- first[1] # should be same day as "second"
  first_col <- first[2]
  second <- split_name(second)
  second_col <- second[2]
  
  # Store the column names and the p-value as a list
  wilcox_p_values[[length(wilcox_p_values) + 1]] <- list(
    day = day,
    col1 = first_col,
    col2 = second_col,
    p_value = wilcox_result$p.value
  )
}

# we compare to the controls 1 by 1
for (control_col_i in seq(1, length(colnames(data)), by = 11)) {
  for (i in seq(control_col_i + 1 + 1, control_col_i + 10, by = 1)) {
    # start after the day column (i goes from 2)
    # Get the current pair of columns (e.g. col_1 and col_2)
    group_to_compare <- colnames(data)[c(control_col_i, i)]
    # calculate the p-value between the 2 cols extracted
    first <- group_to_compare[1]
    second <- group_to_compare[2]
    wilcox_result = wilcox.test(data[[first]], data[[second]])
    
    first <- split_name(first)
    day <- first[1] # should be same day as "second"
    first_col <- first[2]
    second <- split_name(second)
    second_col <- second[2]
    
    # Store the column names and the p-value as a list
    wilcox_p_values[[length(wilcox_p_values) + 1]] <- list(
      day = day,
      col1 = first_col,
      col2 = second_col,
      p_value = wilcox_result$p.value
    )
  }
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
            by = c("phage_type", "condition" = "condition_sd")) %>%
  mutate(is_cb = (condition == "tyb" | condition == "mb"))

# prepare significance data for plot
# Convert the list-of-lists into a data frame
df_signif <- do.call(rbind,
                     lapply(wilcox_p_values, function(x) {
                       data.frame(
                         day = x$day,
                         col1 = x$col1,
                         col2 = x$col2,
                         p_value = x$p_value,
                         stringsAsFactors = FALSE
                       )
                     }))

star_data <- df_signif %>%
  # Keep only significant
  filter(p_value <= 0.05) %>%
  rowwise() %>%
  mutate(# Compute the taller absorbance of the two relevant bars
    y_star = {
      val1 <- summary_data_long$absorbance[summary_data_long$day == day &
                                             summary_data_long$condition == col1]
      val2 <-
        summary_data_long$absorbance[summary_data_long$day == day &
                                       summary_data_long$condition == col2]
      # Add a small offset so the star is above the bar
      max(val1, val2) + 0.05
    }) %>%
  ungroup()

patterns <-
  setNames(rep("none", length(unique(
    summary_data_long$condition
  ))), unique(summary_data_long$condition))
# patterns[c("tyb", "mb")] <- c("stripe", "crosshatch")

# change materials names for nice x axis
summary_data_long$phage_type[summary_data_long$phage_type == "D22"] <-
  "FPSV-D22"
summary_data_long$phage_type[summary_data_long$phage_type == "D45"] <-
  "FPSV-D45"

#############################

plot <- ggplot(
  summary_data_long,
  aes(
    x = interaction(phage_type, condition),
    # factor(phage_type),
    y = absorbance,
    fill = condition,
    pattern = condition
  )
) +
  geom_bar_pattern(
    pattern_fill = "gray50",
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
    limits = c(0, 12),
    expand = c(0, 0),
    breaks = seq(0, 12, by = 1),
    #labels = scales::scientific
  ) +
  geom_errorbar(
    aes(
      ymin = ifelse(absorbance - sd < 0, 0, absorbance - sd),
      ymax = absorbance + sd
    ),
    position = position_dodge(width = 0.9),
    width = 0.4,
    linewidth = 0.8
  ) +
  scale_fill_manual(name = "",
                    values = color_map,
                    labels = labels) +
  scale_pattern_manual(values = patterns) +
  labs(
    x = "\nF. psychrophilum phages",
    y = expression(bold(log[10] * "(PFU/mL)")),
    fill = "Condition",
    title = "Susceptability testing of F. psychrophilum"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.x = element_text(colour = "transparent"),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
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
    axis.text = element_text(size = 19)
  ) +
  scale_x_discrete(labels = sapply(levels(
    interaction(summary_data_long$phage_type, summary_data_long$condition)
  ), function(x) {
    parts <- strsplit(x, "\\.")[[1]]  # Split the interaction label
    if (parts[2] %in% c("FPSV-D22", "FPSV-D45")) {
      parts[2]
    } else {
      paste0("\n", parts[2])
    }
  })) +
  guides(fill = guide_legend(nrow = 7,
                             override.aes = list(
                               pattern = c("none",
                                           "none",
                                           "none",
                                           "none",
                                           "none",
                                           "none",
                                           "none")
                             )),
         pattern = "none")

# add significance to graph
# plot <- plot +
#   geom_text(
#     data = star_data,
#     aes(
#       x = factor(day),
#       y = ifelse(col1 == "mb", y_star + 0.5, y_star + 0.2),
#       fill = col2,
#       label = ifelse(col1 == "mb", glue("*{col1}-{col2}"), glue("*{col2}"))
#     ),
#     position = position_dodge(width = 0.9),
#     vjust = -0.5,
#     # Place just above the bar
#     size = 4,
#     color = "red",
#     inherit.aes = FALSE
#   )

# DEBUG: check column names and legend
plot <- plot +
  geom_text(
    data = summary_data_long,
    aes(
      x = factor(phage_type),
      y = absorbance + 0.1,
      fill = condition,
      label = ifelse(absorbance < 0,
                     format(
                       0, scientific = TRUE, digits = 1
                     ),
                     "")
    ),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    # Place just above the bar
    size = 3,
    inherit.aes = FALSE
  )

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
