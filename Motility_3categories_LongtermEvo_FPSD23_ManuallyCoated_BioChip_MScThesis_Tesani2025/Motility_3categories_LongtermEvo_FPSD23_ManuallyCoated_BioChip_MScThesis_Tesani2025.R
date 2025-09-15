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

split_name <- function(string) {
  # get day number from column and remove it for plotting
  matches <- gsub("^d(\\d+)(.*)$", "\\1,\\2", string)
  split_matches <- strsplit(matches, ",")[[1]]
  return(split_matches)
}

data <-
  read.csv(
    'C:\\Users/jasmi/Desktop/October exp/Evo/Motility_3categories_LongtermEvo_FPSD23_ManuallyCoated_BioChip_MScThesis_Tesani2025/Motility_3categories_LongtermEvo_FPSD23_ManuallyCoated_BioChip_MScThesis_Tesani2025.csv'
  )

data <- data[, !grepl("ncb", colnames(data), ignore.case = TRUE)]

unique_conditions <- c("cb")

# Extract day and suffix
data_long <- data %>%
  tidyr::pivot_longer(
    cols = -Motility.category,
    # Exclude Motility.category from pivoting
    names_to = c("variable", "day_suffix"),
    names_pattern = "(d[0-9]+)(.*)"
  ) %>%
  pivot_longer(cols = Motility.category,
               names_to = "category_col",
               values_to = "category_value")
data_wide <- data_long %>%
  pivot_wider(names_from = day_suffix,
              values_from = value)
numeric_parts <- gsub("\\d+", "", names(data_wide))
numeric_parts[1] <- "day"
renamed_data <-
  data_wide %>% rename_with(~ numeric_parts, .cols = names(data_wide))

df_clean <- renamed_data
# remove "d" from day column values (d3 -> 3, d5 -> 5)
df_clean$day <- as.numeric(gsub("d", "", df_clean$day))

# rename columns to make them look nice
labels <- c("Coated",
            "Noncoated")

data_long <- df_clean %>%
  pivot_longer(cols = c("cb"),
               names_to = "condition",
               values_to = "value")

data_long <- data_long %>%
  arrange(day)

#############################

order_by_day <- function() {
  days <- sort(unique(data_long$day))
  unlist(lapply(days, function(d)
    c(paste0(d, "\ncb"))))
}


the_plot <- ggplot(data_long,
                   aes(
                     x = interaction(day, condition, sep = "\n"),
                     y = value,
                     group = condition,
                     fill = factor(category_value)
                   )) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent
  ) +
  scale_x_discrete(
    limits = order_by_day(),
    labels = function(x) {
      # Replace the second label for each day
      sapply(x, function(label) {
        # Check if the label contains the second condition (e.g., "cb" or "ncb")
        if (grepl("\ncb", label)) {
          l <- sub("[a-zA-Z]+", "", label)
          return(l)  # Custom label for the first condition
        } else {
          return("cb ncb       ")  # Custom label for the second condition
        }
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
    values = c("1" = "gray85",
               "2" = "gray65",
               "3" = "gray45"),
    labels = c("No motility", "Moderate motility", "Strong motility")
  ) +
  labs(x = "Time (days)",
       y = "",
       fill = "",
       title = "Coated beads") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
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
  )
# +
#   guides(fill = guide_legend(nrow = 2,
#                              override.aes = list(pattern = c("stripe", "none"))),
#          pattern = "none")

grid.arrange(
  the_plot,
  nrow = 1,
  widths = unit(30, "cm"),
  heights = unit(20, "cm")
)
