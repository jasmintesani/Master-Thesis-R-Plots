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
    'C:\\Users/jasmi/Desktop/October exp/Evo/CFUperOD595_ShorttermEvo_FPSD23_CoatedBioChips_MScThesisTesani2025/CFUperOD595_ShorttermEvo_FPSD23_CoatedBioChips_MScThesisTesani2025.csv'
  )

# extract controls
ncnbs <- data[, grepl("ncnb", colnames(data))]
ncnbs_unrolled <- unlist(ncnbs)
min <- min(ncnbs_unrolled)
max <- max(ncnbs_unrolled)

# filter data for the rest of the code
# define column ids and colors
data <- data[,!grepl("ncnb", colnames(data))]
data <- data[,!grepl("w", colnames(data))]
unique_conditions <- c("cb", "ncb")
colors <- c("gray85", "gray60")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("variable", "day_suffix"),
    names_pattern = "(d[0-9]+)(.*)"
  )

# apply log10 to the data
data_long$value <- sapply(data_long$value, log10)

data_wide <- data_long %>%
  pivot_wider(names_from = day_suffix,
              values_from = value)
numeric_parts <- gsub("\\d+", "", names(data_wide))
numeric_parts[1] <- "day"
renamed_data <-
  data_wide %>% rename_with( ~ numeric_parts, .cols = names(data_wide))

df <- renamed_data
list_cols <- names(df)[sapply(df, is.list)]
# Create new dataframe starting with non-list columns
df_clean <- df[,!names(df) %in% list_cols]
# Add mean columns
for (col in list_cols) {
  df_clean[paste0("mean_", col)] <-
    sapply(df[[col]], function(x)
      mean(x, na.rm = TRUE))
}
# remove "d" from day column values (d3 -> 3, d5 -> 5)
df_clean$day <- as.numeric(gsub("d", "", df_clean$day))

# collect p-values in a list
wilcox_p_values <- list()
# we loop through the columns, 2 by 2, to look at 1 sampling day at a time
for (i in seq(1, ncol(data), by = 2)) {
  # start after the day column (i goes from 2)
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
  mutate(is_cb = !grepl("ncb", condition))  # TRUE if ends with 'cb'

# rename columns to make them look nice
labels <- c("Coated",
            "Noncoated")

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




write_xlsx(summary_data,
           path = glue(
             "C:\\Users/jasmi/Desktop/October exp/Evo/",
             "CFUperOD595_ShorttermEvo_FPSD23_CoatedBioChips_MScThesisTesani2025",
             "/means_and_sds.xlsx"
           ),
           col_names = TRUE)

#############################

plot <- ggplot(summary_data_long,
               aes(
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
  ###### Y-AXIS adjustments
  scale_y_continuous(
    limits = c(0, 7.5),
    expand = c(0, 0),
    breaks = seq(0, 7.5, by = 1),
    # labels = scales::scientific
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
  scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(
    x = "\nTime (days)",
    y = expression(bold(log[10] * "(CFU/bead/Biomass(" * OD[590] * "))")),
    fill = "Condition",
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.title.y = element_text(margin = margin(r = 20), face = "bold"),
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
