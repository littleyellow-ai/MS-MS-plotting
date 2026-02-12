# Load required libraries
library(ggplot2)
library(scales)
library(dplyr)

# Set working directory
setwd("E:/SI-MSMS/MSMS")

# Helper function: Determine if two m/z values match within tolerance
is_match <- function(mz1, mz2, tolerance = 0.005) abs(mz1 - mz2) <= tolerance

# Function to extract numeric data from lines
extract_numeric_data <- function(lines, section_start) {
  # Find data start line
  data_start <- NA
  for (i in seq_along(lines)) {
    if (grepl("^\\d+[\t ]+\\d+", lines[i]) || grepl("^\\d+\\.\\d+[\t ]+\\d+", lines[i]) ||
        grepl("^-?\\d+\\.?\\d*[\t ]+-?\\d+\\.?\\d*", lines[i])) {
      data_start <- i
      break
    }
  }
  
  if (is.na(data_start)) return(NULL)
  
  # Extract valid numeric lines
  valid_lines <- character(0)
  for (line in lines[data_start:length(lines)]) {
    if (grepl("^\\d+[\t ]+\\d+", line) || grepl("^\\d+\\.\\d+[\t ]+\\d+", line) ||
        grepl("^-?\\d+\\.?\\d*[\t ]+-?\\d+\\.?\\d*", line)) {
      valid_lines <- c(valid_lines, line)
    } else {
      break
    }
  }
  
  if (length(valid_lines) == 0) return(NULL)
  
  data <- read.table(text = valid_lines, header = FALSE, col.names = c("mz", "intensity"))
  data <- data[, 1:2]
  data$mz <- as.numeric(as.character(data$mz))
  data$intensity <- as.numeric(as.character(data$intensity))
  na.omit(data)
}

# Create output folder
dir.create("MS2_plots", showWarnings = FALSE)

# Process each file
for (file in list.files(path = "problem5", pattern = "\\.txt$", full.names = TRUE)) {
  tryCatch({
    lines <- readLines(file)
    
    # Find section starts
    standard_start <- sample_start <- NA
    for (i in seq_along(lines)) {
      line_lower <- tolower(lines[i])
      if (grepl("^standard", line_lower) && is.na(standard_start)) standard_start <- i
      if (grepl("^sample", line_lower) && is.na(sample_start)) sample_start <- i
    }
    
    if (is.na(standard_start) || is.na(sample_start)) {
      warning(paste("Standard or sample marker not found:", basename(file)))
      next
    }
    
    # Extract data
    standard_data <- extract_numeric_data(lines[standard_start:ifelse(sample_start > standard_start, sample_start - 1, length(lines))], "standard")
    sample_data <- extract_numeric_data(lines[sample_start:length(lines)], "sample")
    
    if (is.null(standard_data) || is.null(sample_data) || nrow(standard_data) == 0 || nrow(sample_data) == 0) {
      warning(paste("No valid data:", basename(file)))
      next
    }
    
    # Process data
    standard_data$intensity <- -abs(standard_data$intensity)
    
    pos_data <- sample_data[sample_data$intensity > 0, ]
    neg_data <- standard_data[standard_data$intensity < 0, ]
    
    pos_data$rel_intensity <- pos_data$intensity / max(abs(pos_data$intensity))
    neg_data$rel_intensity <- neg_data$intensity / max(abs(neg_data$intensity))
    
    data_norm <- rbind(
      data.frame(pos_data[, c("mz", "rel_intensity")], source = "sample"),
      data.frame(neg_data[, c("mz", "rel_intensity")], source = "standard")
    )
    
    # Find matching fragments
    matched_pairs <- list()
    if (nrow(pos_data) > 0 && nrow(neg_data) > 0) {
      for (i in 1:nrow(pos_data)) {
        for (j in 1:nrow(neg_data)) {
          if (is_match(pos_data$mz[i], abs(neg_data$mz[j]))) {
            matched_pairs[[length(matched_pairs) + 1]] <- list(
              mz_sample = pos_data$mz[i],
              mz_standard = neg_data$mz[j],
              rel_intensity_sample = pos_data$rel_intensity[i],
              rel_intensity_standard = neg_data$rel_intensity[j]
            )
          }
        }
      }
    }
    
    # Select top 2 matched fragments
    top_matched <- list()
    if (length(matched_pairs) > 0) {
      pair_intensities <- sapply(matched_pairs, function(pair) {
        abs(pair$rel_intensity_sample) + abs(pair$rel_intensity_standard)
      })
      top_count <- min(2, length(pair_intensities))
      top_indices <- order(pair_intensities, decreasing = TRUE)[1:top_count]
      top_matched <- matched_pairs[top_indices]
    }
    annotation_data <- data.frame()
    label_data <- data.frame()
    if (length(matched_pairs) > 0) {
      circle_positions <- sapply(matched_pairs, function(pair) pair$mz_sample)
      annotation_data <- data.frame(
        mz = circle_positions,
        y = 0,
        type = "circle"
      )
    }
    if (length(top_matched) > 0) {
      for (pair in top_matched) {
        pos_base_y <- pair$rel_intensity_sample
        pos_y <- min(pos_base_y + 0.08, 1.4)
        label_data <- rbind(label_data, data.frame(
          mz = pair$mz_sample,
          y = pos_y,
          label = sprintf("%.5f", pair$mz_sample),
          side = "sample",
          vjust = -0.5,
          color = "black"
        ))
        neg_base_y <- pair$rel_intensity_standard
        neg_y <- max(neg_base_y - 0.08, -1.4)
        
        label_data <- rbind(label_data, data.frame(
          mz = pair$mz_standard,
          y = neg_y,
          label = sprintf("%.5f", abs(pair$mz_standard)),
          side = "standard",
          vjust = 1.5,
          color = "red"
        ))
      }
    }
    
    # Calculate plot parameters
    all_data <- rbind(sample_data, standard_data)
    min_mz <- min(all_data$mz, na.rm = TRUE)
    max_mz <- max(all_data$mz, na.rm = TRUE)
    
    has_high_mz <- any(all_data$mz > 400, na.rm = TRUE)
    has_low_mz <- any(all_data$mz < 80, na.rm = TRUE)
    
    # Determine x-axis ticks
    if (has_high_mz) {
      start <- floor(min_mz / 100) * 100
      x_ticks <- seq(start, max_mz, by = 100)
    } else if (max_mz > 200) {
      start <- floor(min_mz / 50) * 50
      x_ticks <- seq(start, max_mz, by = 50)
    } else {
      x_ticks <- seq(ifelse(has_low_mz, 25, floor(min_mz / 25) * 25), max_mz, by = 25)
    }
    
    if (min_mz < 25 && has_low_mz) x_ticks <- unique(c(0, x_ticks))
    if (length(x_ticks) > 15) x_ticks <- x_ticks[seq(1, length(x_ticks), by = 2)]
    
    # Calculate x-axis limits
    x_limit_left <- 0
    if (nrow(label_data) > 0) {
      max_label_length <- max(nchar(label_data$label))
      padding <- ifelse(has_high_mz, 0.10, 0.03 * max_label_length) * (max_mz - min_mz)
      x_limit_right <- max_mz + padding
    } else {
      x_limit_right <- max_mz
    }
    
    # Create plot
    bar_width <- ifelse(has_high_mz, 1.2, 0.7)
    plot <- ggplot(data_norm, aes(x = mz, y = rel_intensity, fill = source)) +
      geom_bar(stat = "identity", width = bar_width, size = ifelse(has_high_mz, 0.3, 0)) +
      scale_fill_manual(values = c(sample = "black", standard = "red"), guide = "none")
    
    # Add annotations
    if (nrow(annotation_data) > 0) {
      plot <- plot + geom_point(data = annotation_data, aes(x = mz, y = y), 
                                shape = 21, size = 1.0, color = "gray40", 
                                fill = "lightblue", alpha = 0.5, stroke = 0.3)
    }
    
    # Add labels
    if (nrow(label_data) > 0) {
      sample_labels <- label_data[label_data$side == "sample", ]
      standard_labels <- label_data[label_data$side == "standard", ]
      
      if (nrow(sample_labels) > 0) {
        plot <- plot + geom_text(data = sample_labels, aes(x = mz, y = y, label = label),
                                 size = 2.2, vjust = -0.5, hjust = 0.5,
                                 color = "black", fontface = "bold")
      }
      
      if (nrow(standard_labels) > 0) {
        plot <- plot + geom_text(data = standard_labels, aes(x = mz, y = y, label = label),
                                 size = 2.2, vjust = 1.5, hjust = 0.5,
                                 color = "red", fontface = "bold")
      }
    }
    
    # Add plot elements
    plot <- plot +
      labs(x = "m/z", y = "Relative intensity", title = gsub("\\.[^.]+$", "", basename(file))) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
        axis.title.x = element_text(size = 9, margin = margin(t = 5)),
        axis.title.y = element_text(size = 9, margin = margin(r = 5)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(size = 0.5, color = "black"),
        axis.ticks = element_line(size = 0.4, color = "black"),
        axis.ticks.length = unit(0.1, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
      ) +
      scale_x_continuous(
        breaks = x_ticks,
        limits = c(x_limit_left, x_limit_right),
        expand = expansion(mult = c(0, ifelse(has_high_mz, 0.02, 0.1)))
      ) + 
      scale_y_continuous(
        breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
        limits = c(-1.4, 1.4),
        expand = expansion(mult = c(0, 0)),
        labels = function(x) format(x, digits = 1)
      ) +
      geom_hline(yintercept = 0, color = "black", size = 0.5)
    
    # Save plot
    ggsave(paste0("MS2_plots/", gsub("\\.[^.]+$", "", basename(file)), ".png"),
           plot, width = 7.26/2.54, height = 6.6/2.54, dpi = 300, units = "in")
    
    # Output summary
    cat(paste("Processed:", basename(file), 
              "| Sample:", nrow(sample_data), 
              "| Standard:", nrow(standard_data), 
              "| Matches:", length(matched_pairs), 
              "| Labels:", nrow(label_data), 
              "| m/z range:", round(min_mz, 2), "-", round(max_mz, 2), "\n"))
    
  }, error = function(e) {
    warning(paste("Failed:", basename(file), "Error:", e$message))
  })
}
