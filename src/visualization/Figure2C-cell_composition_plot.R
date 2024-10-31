# cell_composition_plot.R
# Create stacked bar plot of cell type composition
# Data from Qiu et al., 2021

library(tidyr)
library(ggplot2)
library(dplyr)

# Define color scheme for cell types
CELL_COLORS <- c(
  'B' = '#0868ac',           # blue
  'CD14+ Mono' = '#ff7f00',  # orange
  'CD4+ T' = '#33a02c',      # light green
  'CD8+ T' = '#969696',      # light grey
  'NK' = '#6a51a3',          # purple
  'DC' = '#006d2c',          # dark green
  'Erythroid precursor' = '#ef3b2c',  # red
  'FCGR3A+ Mono' = '#feb24c',  # light orange
  'Neutrophil' = '#fb9a99',    # pink
  'CMP' = '#525252',           # dark grey
  'Platelet' = '#ffeda0'       # yellow
)

#' Process percentage data and calculate condition averages
#' @param input_file Path to percentage data file
#' @return Processed dataframe with condition averages
process_percentages <- function(input_file) {
  # Read data
  percentages <- read.delim(input_file, check.names = FALSE)
  
  # Ensure numeric columns
  percentages <- percentages %>%
    mutate(across(-sample, as.numeric))
  
  # Calculate condition averages
  percentages %>%
    mutate(
      HC = (HC1 + HC2) / 2,
      S = (S1_T0 + S1_T6 + S2_T0 + S2_T6 + S3_T0 + S3_T6) / 6,
      NS = (NS1_T0 + NS1_T6 + NS2_T0 + NS2_T6) / 4
    ) %>%
    select(CellType = sample, HC, S, NS) %>%
    gather("Condition", "Percentage", -CellType)
}

#' Create stacked bar plot
#' @param data Processed percentage data
#' @return ggplot object
create_composition_plot <- function(data) {
  ggplot(data, aes(x = Condition, y = Percentage, fill = CellType)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = CELL_COLORS) +
    theme_bw() +
    theme(
      text = element_text(size = 16),
      axis.title = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Save plot with consistent parameters
#' @param plot ggplot object
#' @param output_file Output file path
save_plot <- function(plot, output_file) {
  png(
    file = output_file,
    width = 6,
    height = 6,
    units = 'in',
    res = 300
  )
  print(plot)
  dev.off()
}

#' Main execution function
#' @param input_file Path to percentage data file
#' @param output_dir Output directory
main <- function(input_file, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process data
  composition_data <- process_percentages(input_file)
  
  # Create plot
  plot <- create_composition_plot(composition_data)
  
  # Save plot
  save_plot(
    plot,
    file.path(output_dir, "cell_type_composition.png")
  )
  
  # Save processed data
  write.csv(
    composition_data,
    file.path(output_dir, "cell_type_composition.csv"),
    row.names = FALSE
  )
  
  message("Plot and data saved in: ", output_dir)
  
  return(list(
    data = composition_data,
    plot = plot
  ))
}

# Example usage:
# results <- main(
#   input_file = "results/cluster_composition/cluster_percentages.txt",
#   output_dir = "results/figures"
# )
