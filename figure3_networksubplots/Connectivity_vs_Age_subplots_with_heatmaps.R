# Plotting connectivity versus gestational age for Lanxin's birth transition paper. 
# Author: Tanya Bhatia
# Date: Aug 27 2024

# Libraries
library(R.matlab)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(tidyr)

# Load and prepare data as before
setwd('/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes')

matfile <- readMat("/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/network_FC_vec_new.mat")
df <- matfile$network.avg
df_new <- as.data.frame(df)
df_t <- t(df_new)
network_avg <- df_t
network_avg_df <- as.data.frame(network_avg)

network_avg_combat <- read.csv('/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/network_df_combat_output_edited_as_matrix.csv')
age_data <- read.csv('/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/Data_for_R_n184_final.csv')

conn_age_df <- cbind.data.frame(age_data, network_avg_combat)
names(conn_age_df) <- gsub("^V([0-9]+)\\.combat$", "\\1", names(conn_age_df))

conn_age_df <- conn_age_df[!is.na(conn_age_df$Preterm),]
filtered_data <- subset(conn_age_df, abs(mean_pry)<=0.5 & abs(mean_xyz)<=0.5 & abs(mean_max_pry)<=1 & abs(mean_max_xyz)<=1)


# Step 1: Define the original mapping of numbers to network pairs
# (You'll need to define this according to your specific network pairs)
pair_mapping <- c(
  "RS-RS", "RS-LS", "RS-LT", "RS-Occ", "RS-Sub", "RS-Inf", "RS-Sup", "RS-RT",
  "LS-LS", "LS-LT", "LS-Occ", "LS-Sub", "LS-Inf", "LS-Sup", "LS-RT",
  "LT-LT", "LT-Occ", "LT-Sub", "LT-Inf", "LT-Sup", "LT-RT",
  "Occ-Occ", "Occ-Sub", "Occ-Inf", "Occ-Sup", "Occ-RT",
  "Sub-Sub", "Sub-Inf", "Sub-Sup", "Sub-RT",
  "Inf-Inf", "Inf-Sup", "Inf-RT",
  "Sup-Sup", "Sup-RT",
  "RT-RT")
  
# Step 2: Define the new order (where RT is next to LT)
new_order <- c(
  "RS-RS", "RS-LS", "RS-LT", "RS-RT", "RS-Occ", "RS-Sub", "RS-Inf", "RS-Sup",
  "LS-LS", "LS-LT", "LS-RT", "LS-Occ", "LS-Sub", "LS-Inf", "LS-Sup",
  "LT-LT", "LT-RT", "LT-Occ", "LT-Sub", "LT-Inf", "LT-Sup",
  "RT-RT", "RT-Occ", "RT-Sub", "RT-Inf", "RT-Sup",
  "Occ-Occ", "Occ-Sub", "Occ-Inf", "Occ-Sup",
  "Sub-Sub", "Sub-Inf", "Sub-Sup",
  "Inf-Inf", "Inf-Sup",
  "Sup-Sup")

# Function to convert a pair to a set-like representation
to_set <- function(pair) {
  elements <- unlist(strsplit(pair, "-"))
  sorted_elements <- sort(elements)
  return(paste(sorted_elements, collapse = "-"))
}

# Convert each pair in pair_mapping to a set-like representation
set_pair_mapping <- sapply(pair_mapping, to_set)
set_new_order <- sapply(new_order, to_set)

# Step 3: Create a vector to map the new order to original column indices
index_mapping <- match(set_new_order, set_pair_mapping)

# Create lists to store the ggplot objects
line_plots <- list()
heat_maps <- list()

# Loop to generate line plots and heatmaps
for (i in 1:36) {
  
  col_index <- index_mapping[i]
  
  # Line plots
  # Iterate col_index for new order with index mapping (e.g. changed order of networks)
  correlation_result <- cor.test(conn_age_df[[as.character(col_index)]], conn_age_df$GA_scan, 
                                 method = 'pearson', use = 'complete.obs')
  # Iterate i for original mapping 
  #correlation_result <- cor.test(conn_age_df[[as.character(i)]], conn_age_df$GA_scan, 
  #                               method = 'pearson', use = 'complete.obs')
  r <- correlation_result$estimate
  p_value <- correlation_result$p.value
  
  title_text <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "")))
  
  line_plot <- ggplot(conn_age_df, aes(x = GA_scan, y = !!sym(as.character(col_index)))) +
    geom_point(aes(colour = factor(Fetal)), alpha = 0.8, size = 0.3) +
    scale_color_discrete(name = "", type = c("gray31", "gray71")) +
    geom_smooth(method = 'gam', size = 0.45, fill = 'palevioletred1', color = 'violetred3') +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
          legend.position = "none",
          plot.title = element_text(size = 7, face = "bold"),
          axis.text = element_text(size = 6.4),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    #scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_x_continuous(breaks = c(30,40,50)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))#, labels = scales::label_number(accuracy = 0.01))
  
  # Heatmaps
  col_name <- as.character(col_index)
  y <- filtered_data[[col_name]]
  x <- filtered_data$GA_scan
  loess_fit <- loess(y ~ x)
  
  x0 <- seq(25, 55, length.out = 150)
  y_pred <- predict(loess_fit, data.frame(x = x0))
  dy_dx <- diff(y_pred) / diff(x0)
  
  df <- data.frame(values = dy_dx, position = 1:149)
  df <- drop_na(df)
  df$values[df$values > 0.004] <- 0.004
  df$values[df$values < -0.004] <- -0.004
  
  heat_map <- ggplot(df, aes(x = position, y = 1)) +
    geom_raster(aes(fill = values), interpolate = TRUE) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, limits = c(-0.004, 0.004)) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(0, 0, 0, 17, unit = "pt"),
          aspect.ratio = 0.14 / 0.6) +  # Maintain aspect ratio for height and width # 0.14/0.69
    coord_fixed(ratio = 0.14 / 0.6)  # Fix the aspect ratio exactly
  
  # Combine the line plot and heat map
  combined_plot <- plot_grid(line_plot, heat_map, ncol = 1, rel_heights = c(0.86, 0.14))
  
  # Append combined plot to list
  line_plots[[i]] <- combined_plot
}

# Create a matrix for layout
layout_matrix <- matrix(NA, 8, 8)
layout_matrix[lower.tri(layout_matrix, diag = TRUE)] <- 1:36

# Combine all plots and save
best_plot <- grid.arrange(grobs = line_plots, layout_matrix = layout_matrix)
#print(best_plot)
# Save the plot
setwd('/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/figures')
ggsave("connectivity_graphs_heatplots_RESTRUCTURED.png", best_plot, width = 8, height = 8, limitsize = FALSE)
