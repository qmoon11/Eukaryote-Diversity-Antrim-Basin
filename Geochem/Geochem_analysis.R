# This script processes the chemical and geochemical analyses of the Antrim gradient

#### Load Packages #### 
library(readxl)
library(purrr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(cowplot)
library(geometry)
setwd("/Users/quinnmoon/Downloads/Antrim_Microbiome")
#### Read in metadata ####
metadata <-read.csv(file = "metadata_AS.csv") 
#### Chemistry by depth ####
# Assign custom shapes and colors to the sites. 
#shapes: 21–25 filled + star (8) for the 6th site
custom_shapes <- c(
  "Mortensen" = 21,  # filled circle
  "Conant"   = 22,  # filled square
  "Greg"     = 23,  # filled diamond
  "Bull"     = 24,  # filled triangle up
  "Horn"     = 25,  # filled triangle down
  "West"     = 8    # star (*), open
)

# Fill colors: red for filled shapes, NA for star (open)
custom_fills <- c(
  "Mortensen" = "black",
  "Conant"   = "black",
  "Greg"     = "black",
  "Bull"     = "black",
  "Horn"     = "black",
  "West"     = NA
)


# ----- Salinity -----
# This maps Cl min/max to TDS min/max for secondary axis alignment
#only keep first 6 rows of metadaa
metadata <- metadata[1:6, ]

# Scale factor to align Cl to TDS axis
scale_factor <- max(metadata$TDS, na.rm = TRUE) / max(metadata$Cl, na.rm = TRUE)

# Make sure Sample is a factor with consistent order
metadata$Sample <- factor(metadata$Sample)

# Get 6 unique sample names
sample_names <- levels(metadata$Sample)

# Assign shapes 21–26 to the 6 samples
shape_values <- setNames(21:26, sample_names)
fill_values <- setNames(rep("red", 6), sample_names)

# Spearman correlations
tds_cor <- cor.test(metadata$Depth..m., metadata$TDS, method = "spearman")
cl_cor  <- cor.test(metadata$Depth..m., metadata$Cl, method = "spearman")

# Extract rho and p
tds_rho <- round(tds_cor$estimate, 2)
tds_p   <- signif(tds_cor$p.value, 2)

cl_rho  <- round(cl_cor$estimate, 2)
cl_p    <- signif(cl_cor$p.value, 2)

# Calculate max x and y positions for annotations
# Position for annotation: bottom right inside plot
depth_max_95pct <- max(metadata$Depth..m., na.rm = TRUE) * 0.95
TDS_min <- min(metadata$TDS, na.rm = TRUE)
TDS_range <- max(metadata$TDS, na.rm = TRUE) - TDS_min
TDS_annotation_y1 <- TDS_min + 0.05 * TDS_range    # 5% above min (near bottom)
TDS_annotation_y2 <- TDS_min + 0.10 * TDS_range    # 10% above min


x_ann <- depth_max_95pct + 0.05 * diff(range(metadata$Depth..m.))   # move right by 5% of depth range
y1_ann <- TDS_annotation_y1 - 0.05 * diff(range(metadata$TDS))      # move down by 5% of TDS range
y2_ann <- TDS_annotation_y2 - 0.05 * diff(range(metadata$TDS))      # move down by 5% of TDS range



salinity <- ggplot(metadata, aes(x = Depth..m.)) +
  # Linear fits
  geom_smooth(aes(y = TDS), method = "lm", se = FALSE, color = "blue", size = 1) +
  geom_smooth(aes(y = Cl * scale_factor), method = "lm", se = FALSE, color = "red", size = 1) +
  
  # Points
  geom_point(aes(y = TDS, shape = Sample), fill = "blue", size = 3, color = "blue", stroke = 1) +
  geom_point(aes(y = Cl * scale_factor, shape = Sample), fill = NA, size = 3, color = "red", stroke = 1) +
  
  # Custom shapes
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  
  # Axes
  scale_y_continuous(
    name = "TDS (mg/L)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cl (mg/L)")
  ) +
  labs(x = "Depth Below Surface (m)") +
  
  # Spearman stats annotation with intermediate font size, correct rho symbol, and moved
  annotate("text", x = x_ann, y = y1_ann,
           label = bquote("TDS: " ~ rho == .(tds_rho) * ", p = " * .(tds_p)),
           hjust = 1, size = 5.5, color = "blue") +
  annotate("text", x = x_ann, y = y2_ann,
           label = bquote("Cl: " ~ rho == .(cl_rho) * ", p = " * .(cl_p)),
           hjust = 1, size = 5.5, color = "red") +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y.left = element_text(color = "blue", size = 18),
    axis.title.y.right = element_text(color = "red", size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    
    plot.margin = margin(5, 5, 5, 5)
  )
#salinity


# ----- pH -----

# Spearman correlations
pH_cor <- cor.test(metadata$Depth..m., metadata$pH, method = "spearman")

# Extract rho and p
pH_rho <- round(pH_cor$estimate, 2)
pH_p   <- signif(pH_cor$p.value, 2)

pH <- ggplot(metadata, aes(x = Depth..m.)) +
  geom_smooth(aes(y = pH), method = "lm", se = FALSE, color = "black", size = 1) +
  geom_point(aes(y = pH, shape = Sample, fill = Sample), size = 3, color = "black", stroke = 1) +
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  scale_fill_manual(name = "Sample", values = custom_fills) +
  labs(x = "Depth Below Surface (m)", y = "pH") +
  annotate(
    "text", x = 556.26, y = 7.93,
    label = bquote(rho == .(pH_rho) * ",  p = " * .(pH_p)),
    hjust = 1, vjust = 1, size = 5.5, color = "black"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )
#pH
# ----- Temp -----
# Spearman correlation for temperature vs depth
temp_cor <- cor.test(metadata$Depth..m., metadata$temperature, method = "spearman")

# Extract rho and p values
temp_rho <- round(temp_cor$estimate, 2)
temp_p <- signif(temp_cor$p.value, 2)

# Set annotation location slightly inside plot range
x_max <- max(metadata$Depth..m., na.rm = TRUE) * 0.99
y_max <- quantile(metadata$temperature, 0.98, na.rm = TRUE)  # Avoid stretching y-axis

temp <- ggplot(metadata, aes(x = Depth..m.)) +
  # Trend line
  geom_smooth(aes(y = temperature), method = "lm", se = FALSE, color = "black", size = 1) +
  
  # Sample points
  geom_point(aes(y = temperature, shape = Sample, fill = Sample), size = 3, color = "black", stroke = 1) +
  
  # Custom shapes/fills
  scale_shape_manual(name = "BNG Well", values = custom_shapes) +
  scale_fill_manual(name = "BNG Well", values = custom_fills) +
  
  # Axis labels
  labs(x = "Depth Below Surface (m)", y = "Temperature (°C)") +
  
  # Spearman annotation (top-right) with intermediate font size and correct rho
  annotate(
    "text", x = x_max, y = y_max,
    label = bquote(rho == .(temp_rho) * ",  p = " * .(temp_p)),
    hjust = 1, vjust = 1, size = 5.5, color = "black"
  ) +
  
  # Theme adjustments for style consistency
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

#temp



# ----- TDN -----
# Spearman correlation for TDN vs depth
# Spearman correlation for TDN vs depth
TDN_cor <- cor.test(metadata$Depth..m., metadata$TDN, method = "spearman")

# Extract rho and p values rounded
TDN_rho <- round(TDN_cor$estimate, 2)
TDN_p <- signif(TDN_cor$p.value, 2)

# Max and min coordinates for annotation placement
x_max <- max(metadata$Depth..m., na.rm = TRUE)
y_min <- min(metadata$TDN, na.rm = TRUE)

TDN <- ggplot(metadata, aes(x = Depth..m.)) +
  geom_smooth(aes(y = TDN), method = "lm", se = FALSE, color = "black", size = 1) +
  geom_point(aes(y = TDN, shape = Sample, fill = Sample), size = 3, color = "black", stroke = 1) +
  scale_shape_manual(name = "BNG Well", values = custom_shapes) +
  scale_fill_manual(name = "BNG Well", values = custom_fills) +
  labs(x = "Depth Below Surface (m)", y = "TDN (mg/L)") +
  annotate(
    "text", x = x_max, y = y_min, 
    label = bquote(rho == .(TDN_rho) * ",  p = " * .(TDN_p)),
    hjust = 1, vjust = 0, size = 5.5, color = "black"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )
#TDN


#### Water Isotope Analysis ####


# Read Michigan GNIP-style reference water data
Anna_reference <- read.csv("Antrim_AntrimCty_Data.csv", stringsAsFactors = FALSE, header = TRUE)

# Prepare published data (no Sample names needed, no error bars)
Anna_reference_grouped <- Anna_reference %>%
  mutate(Sample = "Published Data", group = "Published Data") %>%
  dplyr::select(Sample, delta18O, delta2H, group)

# Prepare metadata (your study), only keep rows with non-missing delta18O
metadata_clean_water <- metadata %>%
  filter(
    !is.na(delta18O),
    !is.na(delta2H),
    !is.na(delta18O.precision..1_.),
    !is.na(delta2H.precision..1_.)
  ) %>%
  mutate(
    delta18O_lower = delta18O - abs(delta18O.precision..1_.),
    delta18O_upper = delta18O + abs(delta18O.precision..1_.),
    delta2H_lower  = delta2H - abs(delta2H.precision..1_.),
    delta2H_upper  = delta2H + abs(delta2H.precision..1_.),
    group = "This Study"
  )

# For combined plotting (includes points from both your study + published data)
metadata_grouped <- metadata_clean_water %>%
  dplyr::select(Sample, delta18O, delta2H, group)

combined_plot_data <- bind_rows(metadata_grouped, Anna_reference_grouped)

# Fill colors
custom_fills_water <- c(
  "Mortensen" = "red",
  "Conant"   = "red",
  "Greg"     = "red",
  "Bull"     = "red",
  "Horn"     = "red",
  "West"     = NA  # Unfilled
)
fill_values <- c(custom_fills_water, "Published Data" = NA)

# Shapes: yours + published data
shape_values <- c(custom_shapes, "Published Data" = 1)

# Plot
water_isotope <- ggplot() +
  # GMWL line
  geom_line(data = gmwl_data, aes(x = delta18O, y = delta2H),
            color = "blue", size = 1.2) +
  
  # Error bars for your data only
  geom_errorbarh(
    data = metadata_clean_water,
    aes(y = delta2H, xmin = delta18O_lower, xmax = delta18O_upper),
    height = 0.5, color = "black"
  ) +
  geom_errorbar(
    data = metadata_clean_water,
    aes(x = delta18O, ymin = delta2H_lower, ymax = delta2H_upper),
    width = 0.5, color = "black"
  ) +
  
  # All points (your data + published data)
  geom_point(data = combined_plot_data,
             aes(x = delta18O, y = delta2H, shape = Sample, fill = Sample),
             size = 3, color = "black", stroke = 1,
             show.legend = FALSE) +
  
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = fill_values) +
  
  labs(
    x = expression(delta^{18}*O~"(‰ VSMOW)"),
    y = expression(delta^{2}*H~"(‰ VSMOW)")
  ) +
  
  # Inline legend items with increased label size
  annotate("segment", x = -5, xend = -3, y = -130, yend = -130, color = "blue", size = 1.2) +
  annotate("text", x = -2.8, y = -130, label = "GMWL", hjust = 0, size = 5.5) +
  geom_point(aes(x = -5, y = -137), shape = 21, fill = "red", color = "black", size = 3, stroke = 1) +
  annotate("text", x = -4.2, y = -137, label = "This Study", hjust = 0, size = 5.5) +
  geom_point(aes(x = -5, y = -144), shape = 21, fill = NA, color = "black", size = 3, stroke = 1) +
  annotate("text", x = -4.2, y = -144, label = "Published Antrim Basin Values", hjust = 0, size = 5.5) +
  
  # Set axis limits so max is 0
  scale_x_continuous(limits = c(NA, 0)) +
  scale_y_continuous(limits = c(NA, 0)) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

water_isotope


#### Carbon Isotope ####

# Assuming metadata has 'Site' column for grouping samples:
carbon_isotope <- ggplot(metadata, aes(x = d13CH4, y = d13C02, shape = Sample, fill = Sample)) +
  # Microbial gas field (Whiticar 1999)
  annotate("rect",
           xmin = -55, xmax = -49,
           ymin = -15, ymax = 30,
           fill = "forestgreen", alpha = 0.3) +
  
  # Thermogenic gas field (Whiticar 1999)
  annotate("rect",
           xmin = -50, xmax = -45,
           ymin = -15, ymax = -5,
           fill = "gray40", alpha = 0.3) +
  
  # Your sample points
  geom_point(size = 3, color = "black", stroke = 1) +
  
  # Custom shapes and fills
  scale_shape_manual(
    name = "Sample Site",
    values = custom_shapes
  ) +
  scale_fill_manual(
    name = "Sample Site",
    values = custom_fills
  ) +
  
  labs(
    x = expression(delta^{13}*CH[4]~"(‰)"),
    y = expression(delta^{13}*CO[2]~"(‰)")
  ) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

carbon_isotope
#The isotope boundaries are based on (Whiticar 1999, Carbon and hydrogen isotope systematics of bacterial formation and oxidation of methane)


#### Tapp ####
metadata_tapp <- metadata %>%
  filter(!is.na(Tapp), !is.na(Tapp_PosError), !is.na(Tapp_NegError)) %>%
  mutate(
    Tapp_upper = Tapp + abs(Tapp_PosError),
    Tapp_lower = Tapp - abs(Tapp_NegError)
  )


# Get depth range from full metadata
depth_min <- min(metadata$Depth..m., na.rm = TRUE)
depth_max <- max(metadata$Depth..m., na.rm = TRUE)

tapp_plot <- ggplot(metadata_tapp, aes(x = Depth..m., y = Tapp)) +
  # Shaded origin regions
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 80,
           fill = "forestgreen", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 80, ymax = 120,
           fill = "gray40", alpha = 0.3) +
  
  # Vertical error bars
  geom_linerange(aes(ymin = Tapp_lower, ymax = Tapp_upper),
                 color = "black", width = 0.5) +
  
  # Points
  geom_point(aes(shape = Sample, fill = Sample),
             size = 3, color = "black", stroke = 1) +
  
  # Custom styling, no legend for shape/fill
  scale_shape_manual(values = custom_shapes, guide = "none") +
  scale_fill_manual(values = custom_fills, guide = "none") +
  
  labs(
    x = "Depth Below Surface (m)",
    y = "Apparent Formation Temperature (°C)"
  ) +
  coord_cartesian(xlim = c(depth_min, depth_max), ylim = c(0, 120)) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )
#tapp_plot




#### Combine all plots into a single layout ####

# Remove legends from plots
salinity_nolegend <- salinity + theme(legend.position = "none")
pH_nolegend <- pH + theme(legend.position = "none")
TDN_nolegend <- TDN + theme(legend.position = "none")
temp_nolegend <- temp + theme(legend.position = "none")
tapp_nolegend <- tapp_plot + theme(legend.position = "none")
carbon_isotope_nolegend <- carbon_isotope + theme(legend.position = "none")

# Extract Sample Site (shape/fill) legend
site_legend <- cowplot::get_legend(
  TDN +
    guides(shape = guide_legend(ncol = 1)) +  # Set to 1 column (vertical)
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "cm")
    )
)

# Build a dummy plot for Gas Origin legend
site_legend <- cowplot::get_legend(
  TDN +
    guides(shape = guide_legend(ncol = 1)) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(1.2, "cm"),
      legend.spacing.y = unit(0.4, "cm")
    )
)

gas_origin_legend_plot <- ggplot(data.frame(x = 1, y = 1, Origin = c("Biogenic", "Thermogenic")),
                                 aes(x = x, y = y, fill = Origin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Gas Origin",
    values = c(
      "Biogenic" = alpha("forestgreen", 0.6),
      "Thermogenic" = alpha("gray40", 0.6)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.y = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

gas_legend <- cowplot::get_legend(gas_origin_legend_plot)

combined_legends <- plot_grid(site_legend, gas_legend, ncol = 2, rel_heights = c(1, 1))

combined_plot_6panel <- (
  (salinity_nolegend | pH_nolegend) /
    (TDN_nolegend | temp_nolegend) /
    (tapp_plot | carbon_isotope_nolegend) /
    ((water_isotope | combined_legends) + plot_layout(widths = c(3, 1)))
) +
  plot_annotation(
    tag_levels = list(c("(A)  1.", "(A)  2.", "(A)  3.", "(A)  4.", "(B)", "(C)", "(D)"))
  ) &
  theme(plot.tag = element_text(size = 25, face = "bold")) # This sets the tag size!

#set dimensions
combined_plot_6panel <- combined_plot_6panel + plot_layout(heights = c(1, 1, 1, 1))

combined_plot_6panel
#save plot as pdf
ggsave("combined_geochem.pdf",
       plot = combined_plot_6panel,
       width = 18, height = 24, units = "in", device = "pdf", limitsize = FALSE)


