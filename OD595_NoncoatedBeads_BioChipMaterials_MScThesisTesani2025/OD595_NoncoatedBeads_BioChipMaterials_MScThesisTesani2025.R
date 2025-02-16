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
    'C:\\Users/Filippo/Desktop/October exp/Evo/OD595_NoncoatedBeads_BioChipMaterials_MScThesisTesani2025/OD595_NoncoatedBeads_BioChipMaterials_MScThesisTesani2025.csv'
  )

# extract controls
controls <- data[, grepl("tyb", colnames(data))]
controls_unrolled <- unlist(ncnbs)
min <- min(ncnbs_unrolled, na.rm = TRUE)
max <- max(ncnbs_unrolled, na.rm = TRUE)

# filter data for the rest of the code
# define column ids and colors
# data <- data[,!grepl("tyb", colnames(data))]
data <- data[,!grepl("w", colnames(data))]
unique_conditions <-
  c("950106.1.1", "d23", "d27", "d42", "d45", "d49", "s6", "tyb")
colors <-
  c("gray85",
    "gray75",
    "gray65",
    "gray55",
    "gray45",
    "gray35",
    "gray25",
    "gray95")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("variable", "day_suffix"),
    names_pattern = "(d[0-9])(.*)"
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
# remove "d" from day column values (d3 -> 3, d5 -> 5)
df_clean$day <- as.numeric(gsub("d", "", df_clean$day))

# collect p-values in a list
wilcox_p_values <- list()
# we compare to the controls 1 by 1
for (control_col_i in seq(1, length(colnames(data)), by = 8)) {
  for (i in seq(control_col_i + 1, control_col_i + 7, by = 1)) {
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
  mutate(is_cb = condition == "tyb")

# rename columns to make them look nice
labels <-
  c(
    "950106-1/1",
    "FPS-D23",
    "FPS-D27",
    "FPS-D42",
    "FPS-D45",
    "FPS-D49",
    "FPS-S6",
    "TYES broth"
  )

color_map <- c(
  "950106.1.1" = "gray20",
  "d23" = "gray30",
  "d27" = "gray40",
  "d42" = "gray50",
  "d45" = "gray60",
  "d49" = "gray70",
  "mb" = "white",
  "pf430.3" = "gray80",
  "s11a" = "gray90",
  "s6" = "gray100",
  "tyb" = "white"
)

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
  "ncb"  =  0.2,
  "cbw"  = -0.2,
  "ncbw" =  0.2
)

summary_data_long$absorbance[summary_data_long$day == 7 &
                               summary_data_long$condition == "ncb"]
summary_data_long$sd[summary_data_long$day == 7 &
                       summary_data_long$condition == "ncb"]

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
      max(val1 + sd1, val2 + sd2)
    },
    # choose bracket labels
    label_ast = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# save star_dat to CSV
csv_data <- star_data %>% select(day, col1, col2, p_value, label_ast)
write.csv(
  csv_data,
  paste0(
    "C:\\Users/Filippo/Desktop/October exp/Evo",
    "/OD595_NoncoatedBeads_BioChipMaterials_MScThesisTesani2025",
    "/significance_table.csv"
  ),
  row.names = FALSE
)



write_xlsx(summary_data_long,
           path = glue(
             "C:\\Users/Filippo/Desktop/October exp/Evo/",
             "OD595_NoncoatedBeads_BioChipMaterials_MScThesisTesani2025",
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
    limits = c(0, 2),
    expand = c(0, 0),
    breaks = seq(0, 2, by = 0.1),
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
  labs(x = "\nTime (days)",
       y = "OD at 590 nm\n",
       fill = "Condition",
       title = "BioChip beads") +
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
    axis.text = element_text(size = 22)
  ) +
  guides(fill = guide_legend(nrow = 8,
                             override.aes = list(
                               pattern = c("none",
                                           "none",
                                           "none",
                                           "none",
                                           "none",
                                           "none",
                                           "none",
                                           "stripe")
                             )),
         pattern = "none")

grid.arrange(
  plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
