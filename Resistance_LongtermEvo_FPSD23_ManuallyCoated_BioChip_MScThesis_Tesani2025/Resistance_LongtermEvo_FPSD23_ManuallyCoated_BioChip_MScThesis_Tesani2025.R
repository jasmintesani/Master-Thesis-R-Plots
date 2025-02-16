library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpattern)
library(stringr)
library(ggsignif)
library(glue)
library(gridExtra)
library(ggpubr)

data <-
  read.csv(
    'C:\\Users/Filippo/Desktop/October exp/Evo/Resistance_LongtermEvo_FPSD23_ManuallyCoated_BioChip_MScThesis_Tesani2025/Resistance_LongtermEvo_FPSD23_ManuallyCoated_BioChip_MScThesis_Tesani2025.csv'
  )

# remove either "cb"s or "ncb"s from data
data <- data[, !grepl("ncb", colnames(data), ignore.case = TRUE)]

unique_conditions <- c("cb1",
                       # "ncb1",
                       "cb2",
                       # "ncb2",
                       "cb3") # "ncb3")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = -Resistance.category,
    # Exclude Resistance.category from pivoting
    names_to = c("variable", "day_suffix"),
    names_pattern = "(d[0-9]+)(.*)"
  ) %>%
  pivot_longer(cols = Resistance.category,
               names_to = "category_col",
               values_to = "category_value")
data_wide <- data_long %>%
  pivot_wider(names_from = day_suffix,
              values_from = value)
# numeric_parts <- gsub("\\d", "", names(data_wide))
numeric_parts <- names(data_wide)
numeric_parts[1] <- "day"
renamed_data <-
  data_wide %>% rename_with( ~ numeric_parts, .cols = names(data_wide))

df_clean <- renamed_data
# remove "d" from day column values (d3 -> 3, d5 -> 5)
df_clean$day <- as.numeric(gsub("d", "", df_clean$day))

# rename columns to make them look nice
labels <- c("Coated",
            "Noncoated")

data_long <- df_clean %>%
  pivot_longer(cols = unique_conditions,
               names_to = "condition",
               values_to = "value")

data_long <- data_long %>%
  arrange(day)

#############################

order_by_day <- function() {
  days <- sort(unique(data_long$day))
  unlist(lapply(days, function(d)
    c(
      # paste0(d, "\ncb1"),
      paste0(d, "\ncb1"),
      # paste0(d, "\ncb2"),
      paste0(d, "\ncb2"),
      # paste0(d, "\ncb3")
      paste0(d, "\ncb3")
    )))
}


the_plot <- ggplot(data_long,
                   aes(
                     x = interaction(day, condition, sep = "\n"),
                     y = value,
                     group = condition,
                     fill = factor(category_value)
                   )) +
  scale_x_discrete(
    limits = order_by_day(),
    labels = function(x) {
      sapply(x, function(label) {
        if (grepl("cb2", label)) {
          return(sub("^(\\d+).*", "\\1", label))
        }
        return("")
      })
    }
  ) +  # This will group by day
  geom_col(
    position = "fill",
    stat = "identity",
    width = 0.8,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_brewer(name = "Category") +
  scale_fill_manual(
    values = c(
      "1" = "gray90",
      "2" = "gray65",
      "3" = "gray45",
      "4" = "gray20"
    ),
    labels = c(
      "Sensitive",
      "Strongly reduced growth in phage track",
      "Slightly reduced growth in phage track",
      "Resistant"
    )
  ) +
  labs(x = "\nTime (days)",
       y = "",
       fill = "",
       title = "Coated beads") +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),
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
  )

grid.arrange(
  the_plot,
  nrow = 1,
  widths = unit(45, "cm"),
  heights = unit(20, "cm")
)
