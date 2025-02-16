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
    'C:\\Users/Filippo/Desktop/October exp/Evo/CFUperBead_differentMOIs_FPSD23_Resistance_MScThesisTesani2025/CFUperBead_differentMOIs_FPSD23_Resistance_MScThesisTesani2025.csv'
  )

# remove weird "cbn" columns
data <- data[, !grepl("cbn", colnames(data), ignore.case = TRUE)]

# filter data for the rest of the code
# define column ids and colors
unique_conditions <-
  c("cbMOI0",
    "cbMOI00",
    "ncbMOI0",
    "ncbMOI00")
color_map <- c("cbMOI0" = "gray90",
               "cbMOI00" = "gray70",
               "ncbMOI0" = "gray50",
               "ncbMOI00" = "gray30")
# rename columns to make them look nice
labels <-
  c("Coated MOI 0",
    "Coated MOI 100",
    "Noncoated MOI 0",
    "Noncoated MOI 100")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("variable", "day_suffix"),
    names_pattern = "(^d[0-9]+)(.*)"
  )
data_wide <- data_long %>%
  pivot_wider(names_from = day_suffix,
              values_from = value)

df <- data_wide %>% rename(day = variable)

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
            by = c("day", "condition" = "condition_sd")) %>%
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
patterns[c("cbMOI0", "cbMOI00")] <- c("stripe", "stripe")

# change materials names for nice x axis
summary_data_long$day <-
  gsub("d", "", summary_data_long$day, ignore.case = TRUE)


# write_xlsx(summary_data,
#            path = glue(
#              "C:\\Users/Filippo/Desktop/October exp/Evo",
#              "/CFUperBead_differentMOIs_FPSD23_Resistance_MScThesisTesani2025/",
#              "CFUperBead_differentMOIs_FPSD23_Resistance_MScThesisTesani2025.xlsx"
#            ),
#            col_names = TRUE)

#############################

plot <- ggplot(summary_data_long,
               aes(
                 x = factor(as.numeric(day)),
                 y = absorbance,
                 fill = condition,
                 pattern = condition
               )) +
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
    limits = c(0, 8.7e+06),
    expand = c(0, 0),
    breaks = seq(0, 8.7e+06, by = 1e+06),
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
  labs(x = "\nTime (days)",
       y = "CFU/bead\n",
       fill = "Condition",
       title = "") +
  theme_minimal() +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
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
  guides(fill = guide_legend(nrow = 4,
                             override.aes = list(pattern = c(
                               "stripe",
                               "stripe",
                               "none",
                               "none"
                             ))),
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
# plot <- plot +
#   geom_text(
#     data = summary_data_long,
#     aes(
#       x = factor(as.numeric(day)),
#       y = absorbance,
#       fill = condition,
#       label = condition
#     ),
#     position = position_dodge(width = 0.9),
#     vjust = -0.5,
#     # Place just above the bar
#     size = 3,
#     inherit.aes = FALSE
#   )

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(25, "cm"),
  heights = unit(26, "cm")
)
