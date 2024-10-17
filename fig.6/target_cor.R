# Set the working directory --------------------------------------------------
setwd("D:/bilibiliR/45_Immune_Infiltration_Correlation_Lollipop_Plot")

# Install required packages --------------------------------------------------
# Check if ggplot2 is installed
if (!require("ggplot2", quietly = TRUE)) {
  # If not installed, use install.packages function to install
  install.packages("ggplot2")
}

# Load the necessary packages ------------------------------------------------
library(ggplot2)

# Load example data ----------------------------------------------------------
data <- read.csv("correlation_result.csv", row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Set the target gene
target_gene <- "SHISA3"
mydata <- data[data$Gene == target_gene,]  # Select data for the target gene
head(mydata)

# Data processing ------------------------------------------------------------
# Use reorder function to rearrange the y-axis
mydata$im_cell <- reorder(mydata$im_cell, mydata$Cor)

# Plotting -------------------------------------------------------------------
plot <- ggplot(mydata, aes(x = Cor, y = im_cell)) +  # Create scatter plot
  geom_segment(aes(x = 0, xend = Cor, y = im_cell, yend = im_cell), color = "black") +  # Add segment layer
  geom_point(aes(size = abs(Cor), colour = p.value), alpha = 0.5) +  # Add point layer
  scale_colour_gradient(low = "#339D5A", high = "#EDB306") +  # Set color gradient
  # Alternative color options: #ED063F, #339D5A, #EDB306
  scale_size_continuous(range = c(2, 10)) +  # Set the range for point sizes
  theme_minimal() +  # Set the theme
  # theme(legend.position = "none") +  # Remove legend
  # Adjust the axis settings
  # scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +  # Set x-axis
  # scale_y_discrete(limits = rev(levels(mydata$im_cell))) +  # Set y-axis
  # Customize text style
  theme(axis.line = element_line(size = 1.0),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(hjust = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center the title
        axis.title = element_text(size = 14, face = "bold")) +
  # Add labels
  labs(title = "Correlation of SHISA3 with Immune Cells", 
       x = "Cor", y = "im_cell")

plot  # Display the plot

# Apply different themes -----------------------------------------------------
plot + theme_grey() + ggtitle("theme_grey()")
plot + theme_bw() + ggtitle("theme_bw()")
plot + theme_linedraw() + ggtitle("theme_linedraw()")
plot + theme_light() + ggtitle("theme_light()")
plot + theme_dark() + ggtitle("theme_dark()")
plot + theme_classic() + ggtitle("theme_classic()")

# Save the plot --------------------------------------------------------------
ggsave("target_cor.png", width = 8, height = 6, dpi = 300)

# Tips:
# To quickly learn ggplot2, check out the ggplotAssist tutorial released on 2023-08-16 
# Install ggplotAssist from GitHub
devtools::install_github("cardiomoon/ggplotAssist")
# After installation, restart RStudio
library(ggplotAssist)
