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
  matches <- gsub("^d(\\d+)(.*)$", "\\1,\\2", string)
  split_matches <- strsplit(matches, ",")[[1]]
  return(split_matches)
}

data <-
  read.csv(
    'C:\\Users/Filippo/Desktop/October exp/Evo/CFUperBead_ShorttermEvoFPSD23_ManuallyCoatedFPSVD22BioChip_MScThesisTesani2025/CFUperBead_ShorttermEvoFPSD23_ManuallyCoatedFPSVD22BioChip_MScThesisTesani2025.csv'
  )

# day 5 is bad, we delete it
data <- data[, !grepl("5", colnames(data))]

# extract controls
ncnbs <- data[, grepl("ncnb", colnames(data))]
ncnbs_unrolled <- unlist(ncnbs)
min <- min(ncnbs_unrolled)
max <- max(ncnbs_unrolled)

# filter data for the rest of the code
# define column ids and colors
data <- data[, !grepl("ncnb", colnames(data))]
data <- data[, !grepl("w", colnames(data))]
unique_conditions <- c("cb", "ncb")
colors <- c("gray85", "gray60")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(cols = everything(), names_to = c("variable", "day_suffix"), names_pattern = "(d[0-9]+)(.*)")

# apply log10 to the data
data_long$value <- sapply(data_long$value, log10)

data_wide <- data_long %>%
  pivot_wider(
    names_from = day_suffix,
    values_from = value
  )
numeric_parts <- gsub("\\d+", "", names(data_wide))
numeric_parts[1] <- "day"
renamed_data <- data_wide %>% rename_with(~ numeric_parts, .cols = names(data_wide))

df <- renamed_data
list_cols <- names(df)[sapply(df, is.list)]
# Create new dataframe starting with non-list columns
df_clean <- df[, !names(df) %in% list_cols]
# Add mean columns
for (col in list_cols) {
  df_clean[paste0("mean_", col)] <- sapply(df[[col]], function(x) mean(x, na.rm = TRUE))
}
# remove "d" from day column values (d3 -> 3, d5 -> 5)
df_clean$day <- as.numeric(gsub("d", "", df_clean$day))

# collect p-values in a list
wilcox_p_values <- list()
# we loop through the columns, 2 by 2, to look at 1 sampling day at a time
for (i in seq(1, ncol(data), by = 2)) { # start after the day column (i goes from 2)
  # Get the current pair of columns (e.g. col_1 and col_2)
  day_group <- colnames(data)[i:(i + 1)]
  # calculate the p-value between the 2 cols extracted
  first <- day_group[1]
  second <- day_group[2]
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

# Add standard deviation columns
for (col in list_cols) {
  df_clean[paste0("sd_", col)] <- sapply(df[[col]], function(x) sd(x, na.rm = TRUE))
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
  mutate(is_cb = !grepl("ncb", condition))  # TRUE if ends with 'cb'

# rename columns to make them look nice
labels <- c( "Coated",
             "Noncoated" )

color_map <- setNames(colors, unique_conditions)

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

## 1) Make sure 'day' is a factor in both datasets with the same levels
summary_data_long$day <- factor(summary_data_long$day,
                                levels = sort(unique(summary_data_long$day)))
df_signif$day <- factor(df_signif$day,
                        levels = levels(summary_data_long$day))

## 2) Define how to shift conditions left or right
condition_positions <- c(
  "cb"   = -0.2,
  "ncb"  =  0.2
)

## 3) In star_data, compute numeric day + condition offset
##    (the factor codes) for xmin and xmax
star_data <- df_signif %>%
  filter(df_signif$p_value < 0.05) %>%
  mutate(
    # underlying numeric code of the factor day
    day_num = as.numeric(day),
    x1 = day_num + condition_positions[col1],
    x2 = day_num + condition_positions[col2],
    # adjust y so that the bracket is over the error bars
    y_star = {
      val1 <- summary_data_long$absorbance[summary_data_long$day == day &
                                             summary_data_long$condition == col1]
      val2 <-
        summary_data_long$absorbance[summary_data_long$day == day &
                                       summary_data_long$condition == col2]
      sd1 <- summary_data_long$sd[summary_data_long$day == day &
                                    summary_data_long$condition == col1]
      sd2 <- summary_data_long$sd[summary_data_long$day == day &
                                    summary_data_long$condition == col2]
      max(val1 + sd1, val2 + sd2) + 0.01
    },
    # choose bracket labels
    label_ast = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Compute max heights dynamically for placing asterisks and brackets
max_heights <- summary_data_long %>%
  group_by(day) %>%
  summarise(ypos = max(absorbance + sd) + 0.3)

offsets <- seq(0, 0.7, length.out = length(unique_conditions))

# Assign offsets based on the ordering in unique_conditions
star_data$yoffset <-
  offsets[match(star_data$col2, unique_conditions)]

star_data <-
  left_join(star_data, max_heights, by = "day")  # Merge max heights

#############################

the_plot <- ggplot(summary_data_long, aes(
  x = factor(day),
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
  geom_signif(
    data = star_data,
    aes(
      xmin = x1,
      xmax = x2,
      y_position = ypos,
      annotations = label_ast,
      group = day,
      label = label_ast,
      tip_length = 0.01
    ),
    manual = TRUE,
    inherit.aes = FALSE
  ) +
  ###### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0, 7.5),
    expand = c(0, 0),
    breaks = seq(0, 7.5, by = 1)
    #labels = scales::scientific
  ) +
  geom_errorbar(aes(
    ymin = ifelse(absorbance - sd < 0, 0, absorbance - sd),
    ymax = absorbance + sd),
    position = position_dodge(width = 0.9),
    width = 0.4,
    linewidth = 0.8) +
  scale_fill_manual(
    name = "",
    values = color_map,
    labels = labels
  ) +
  scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(x = "\nTime (days)",
       y = expression(bold(log[10] * "(CFU/bead)")),
       fill = "Condition",
       title = "") +
  theme_minimal() +
  theme(
    legend.position = c(0.72, 0.98),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.8, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 22)
  ) +
  guides(
    fill = guide_legend(
      nrow = 4,
      override.aes = list(
        pattern = c("stripe", "none")
      )
    ),
    pattern = "none"
  )

grid.arrange(the_plot, nrow = 1, widths = unit(30, "cm"), heights = unit(20, "cm"))
