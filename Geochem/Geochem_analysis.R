# This script processes the isotope and geochemical analyses to produce figure 2 and supplementary figure 7 of the Antrim Basin manuscript
# Script includes code to plot depth on both x and y axis

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
reference_geo <- read.csv(file = "ReferenceGeochemistryIsotopes.csv")
#### Chemistry by depth (on x) ####
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
#only keep first 6 rows of metadate
metadata <- metadata[1:6, ]

# Scale factor to align Cl to TDS axis
scale_factor <- max(metadata$TDS, na.rm = TRUE) / max(metadata$Cl, na.rm = TRUE)

# Make sure Sample is a factor with consistent order
metadata$Sample <- factor(metadata$Sample)

# Get 6 unique sample names
sample_names <- levels(metadata$Sample)

# Assign shapes 21–26 to the 6 samples
shape_values <- setNames(21:26, sample_names)

# Spearman correlations
tds_cor <- cor.test(metadata$Depth..m., metadata$TDS, method = "spearman")
cl_cor  <- cor.test(metadata$Depth..m., metadata$Cl, method = "spearman")

# Extract rho and p
tds_rho <- round(tds_cor$estimate, 2)
tds_p   <- signif(tds_cor$p.value, 2)

cl_rho  <- round(cl_cor$estimate, 2)
cl_p    <- signif(cl_cor$p.value, 2)

# Linear regression for TDS
tds_lm <- lm(TDS ~ Depth..m., data = metadata)
summary(tds_lm)$r.squared  # 0.82

# Linear regression for Cl
cl_lm <- lm(Cl ~ Depth..m., data = metadata)
summary(cl_lm)$r.squared #0.81

# Positions for annotation
depth_max_95pct <- max(metadata$Depth..m., na.rm = TRUE) * 0.95
TDS_min <- min(metadata$TDS, na.rm = TRUE)
TDS_range <- max(metadata$TDS, na.rm = TRUE) - TDS_min
TDS_annotation_y1 <- TDS_min + 0.05 * TDS_range
TDS_annotation_y2 <- TDS_min + 0.10 * TDS_range
x_ann <- depth_max_95pct + 0.05 * diff(range(metadata$Depth..m.))
y1_ann <- TDS_annotation_y1 - 0.05 * diff(range(metadata$TDS))
y2_ann <- TDS_annotation_y2 - 0.05 * diff(range(metadata$TDS))

# Salinity band limits
Cl_max_plot <- max(metadata$Cl, na.rm = TRUE)
TDS_max_plot <- Cl_max_plot * scale_factor
brine_min_Cl <- 50000         # or your threshold for brine
brine_min_TDS <- brine_min_Cl * scale_factor

#Blues for Water Regions
band_data <- data.frame(
  ymin  = c(0,        1000,       35000,        brine_min_TDS),
  ymax  = c(1000,     35000,      brine_min_TDS, TDS_max_plot),
  label = c("Fresh",  "Brackish", "Ocean",      "Brine"),
  fill  = c("aliceblue", "lightblue", "deepskyblue", "dodgerblue3")
)

salinity <- ggplot(metadata, aes(x = Depth..m.)) +
  # Water blocks (bands)
  geom_rect(
    data = band_data,
    inherit.aes = FALSE,
    aes(
      xmin = min(metadata$Depth..m., na.rm = TRUE),
      xmax = max(metadata$Depth..m., na.rm = TRUE),
      ymin = ymin, ymax = ymax,
      fill = label
    ),
    alpha = 0.42
  ) +
  scale_fill_manual(
    name = "Salinity Class",
    values = setNames(band_data$fill, band_data$label)
  ) +
  # Linear fits: TDS black, Cl red
  geom_smooth(aes(y = TDS), method = "lm", se = FALSE, color = "black", size = 1) +
  geom_smooth(aes(y = Cl * scale_factor), method = "lm", se = FALSE, color = "red", size = 1) +
  # Points: TDS black, Cl red outline
  geom_point(aes(y = TDS, shape = Sample), fill = "black", size = 3, color = "black", stroke = 1) +
  geom_point(aes(y = Cl * scale_factor, shape = Sample), fill = NA, size = 3, color = "red", stroke = 1) +
  scale_shape_manual(name = "Sample", values = shape_values) +
  # Axes
  scale_y_continuous(
    name = "TDS (mg/L)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cl (mg/L)")
  ) +
  labs(x = "Depth Below Surface (m)") +
  # Spearman stats annotation
  annotate("text", x = x_ann, y = y1_ann,
           label = bquote("TDS: " ~ rho == .(tds_rho) * ", p = " * .(tds_p)),
           hjust = 1, size = 5.5, color = "black") +
  annotate("text", x = x_ann, y = y2_ann,
           label = bquote("Cl: " ~ rho == .(cl_rho) * ", p = " * .(cl_p)),
           hjust = 1, size = 5.5, color = "red") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y.left = element_text(color = "black", size = 18), # TDS axis now black
    axis.title.y.right = element_text(color = "red", size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

salinity
# ----- pH -----

# Spearman correlations
pH_cor <- cor.test(metadata$Depth..m., metadata$pH, method = "spearman")


# Extract rho and p
pH_rho <- round(pH_cor$estimate, 2)
pH_p   <- signif(pH_cor$p.value, 2)

pH_lm <- lm(pH ~ Depth..m., data = metadata)
summary(pH_lm)$r.squared #[1] 0.621143

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

TDN_lm <- lm(TDN ~ Depth..m., data = metadata)
summary(TDN_lm)$r.squared #0.86

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
#ggsave("combined_geochem.pdf",
       plot = combined_plot_6panel,
       width = 18, height = 24, units = "in", device = "pdf", limitsize = FALSE)































#-------- Inverted Plots (depth on y) --------

# -------- Salinity Plot: TDS, Cl vs Depth (y, reversed) --------
# Get min/max for TDS, Cl, Depth
# 1. Calculate plotting limits and scale factor
tds_limits <- range(metadata$TDS, na.rm = TRUE)
cl_limits  <- range(metadata$Cl, na.rm = TRUE)
depth_limits <- range(metadata$Depth..m., na.rm = TRUE)
scale_factor <- diff(tds_limits) / diff(cl_limits)

metadata$Cl_plot <- (metadata$Cl - cl_limits[1]) * scale_factor + tds_limits[1]

# 2. Define chloride band thresholds and convert to TDS axis
cl_band_limits <- c(0, 250, 18000, 22000, cl_limits[2]) # Cl thresholds
tds_band_limits <- (cl_band_limits - cl_limits[1]) * scale_factor + tds_limits[1]

band_data <- data.frame(
  xmin = tds_band_limits[1:4],
  xmax = tds_band_limits[2:5],
  fill = c("Fresh", "Brackish", "Ocean", "Brine"),
  color = c("aliceblue", "lightblue", "deepskyblue", "dodgerblue3")
) %>%
  mutate(
    ymin = depth_limits[1],
    ymax = depth_limits[2],
    xmin = pmax(xmin, tds_limits[1]),
    xmax = pmin(xmax, tds_limits[2])
  ) %>%
  filter(xmin < xmax)

# 3. Filter reference data for complete cases, then transform
ref_points <- reference_geo %>%
  filter(!is.na(Depth..m.), !is.na(Cl..mg.L.))
ref_points$Cl_plot <- (ref_points$Cl..mg.L. - cl_limits[1]) * scale_factor + tds_limits[1]

# 4. Annotation variables
x_ann <- mean(depth_limits)
y1_ann <- mean(tds_limits)
y2_ann <- mean(tds_limits) + 0.1 * diff(tds_limits)
tds_rho <- 0.85
tds_p <- 0.001
cl_rho <- 0.77
cl_p <- 0.002

# 5. Plot -- salinity bands FIRST!
salinity2 <- ggplot(metadata, aes(y = Depth..m.)) +
  # Salinity bands (background first)
  geom_rect(
    data = band_data,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.42
  ) +
  # Reference Cl points (gray), on rect background
  geom_point(
    data = ref_points,
    inherit.aes = FALSE,
    aes(x = Cl_plot, y = Depth..m.),
    color = "gray50",
    size = 4,
    shape = 21,
    fill = NA,
    stroke = 1,
    alpha = 0.7
  ) +
  # Metadata sample points TDS (black, filled)
  geom_point(aes(x = TDS, shape = Sample), fill = "black", size = 5, color = "black", stroke = 1) +
  # Metadata sample points Cl (red, unfilled)
  geom_point(aes(x = Cl_plot, shape = Sample), fill = NA, size = 5, color = "red", stroke = 1) +
  # Smooth lines for TDS/Cl
  geom_smooth(aes(x = TDS), method = "lm", se = FALSE, color = "black", size = 1) +
  geom_smooth(aes(x = Cl_plot), method = "lm", se = FALSE, color = "red", size = 1) +
  # Manual sample shapes
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  # Salinity band colors
  scale_fill_manual(
    name = "Salinity Class",
    values = setNames(band_data$color, band_data$fill)
  ) +
  # X axis = TDS (top), Cl as secondary (bottom)
  scale_x_continuous(
    name = "Total Dissolved Solids (mg/L)",
    position = "top",
    limits = tds_limits,
    sec.axis = sec_axis(
      ~ (.-tds_limits[1]) / scale_factor + cl_limits[1],
      name = "Chloride (mg/L)",
      breaks = pretty(cl_limits)
    )
  ) +
  # Y axis = Depth, reversed
  scale_y_reverse(name = "Depth Below Surface (m)") +

  # Minimal theme and aesthetic tweaks
  theme_minimal() +
  theme(
    axis.title.x.top = element_text(color = "black", size = 18, vjust = 0),
    axis.text.x.top  = element_text(color = "black", size = 18),
    axis.title.x.bottom = element_text(color = "red", size = 18, vjust = 1),
    axis.text.x.bottom  = element_text(color = "red", size = 18),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

salinity2


# For TDS vs. Depth
lm_tds <- lm(Depth..m. ~ TDS, data = metadata)
summary_tds <- summary(lm_tds)
r2_tds <- summary_tds$r.squared
cat("R² for TDS fit:", r2_tds, "\n")
#TDS R2=0.819

# For Cl_plot vs. Depth
lm_cl <- lm(Depth..m. ~ Cl_plot, data = metadata)
summary_cl <- summary(lm_cl)
r2_cl <- summary_cl$r.squared
cat("R² for Cl_plot fit:", r2_cl, "\n")
#Cl R2=0.810

#Grad axis marks for remainder of plots
salinity_build <- ggplot_build(salinity2)
salinity_ticks <- salinity_build$layout$panel_params[[1]]$y$breaks
salinity_limits <- salinity_build$layout$panel_params[[1]]$y$range

# Also, get the reversed limits for scale_y_reverse:
salinity_limits_rev <- rev(salinity_limits)

# -------- pH Plot: pH vs Depth (y, reversed) --------
# Filter reference_geo for valid pH and depth
ref_pH_points <- reference_geo %>%
  filter(!is.na(pH), !is.na(Depth..m.))

pH2 <- ggplot(metadata, aes(y = Depth..m.)) +
  # Reference pH points—background layer
  geom_point(
    data = ref_pH_points,
    inherit.aes = FALSE,
    aes(x = pH, y = Depth..m.),
    color = "gray50",
    size = 4,
    shape = 21,
    fill = NA,
    stroke = 1,
    alpha = 0.7
  ) +
  # Your data—smoothing line + points
  geom_smooth(aes(x = pH), method = "lm", se = FALSE, color = "black", linetype="dashed", size = 1) +
  geom_point(aes(x = pH, shape = Sample, fill = Sample), size = 5, color = "black", stroke = 1) +
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  scale_fill_manual(name = "Sample", values = custom_fills) +
  scale_x_continuous(name = "pH", position = "top") +
  scale_y_reverse(name = "Depth Below Surface (m)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

pH2
# Fit the linear model
lm_pH <- lm(Depth..m. ~ pH, data = metadata)
r2_pH <- summary(lm_pH)$r.squared
#pH r2=0.621

# -------- Temp Plot: Temperature vs Depth (y, reversed) --------
depth_limits_rev <- c(707.425, 64.02)          # Deepest at bottom, shallowest at top
desired_depth_breaks <- c(200, 400, 600) # The exact y-axis ticks and labels you want

temp2 <- ggplot(metadata_clean, aes(y = Depth..m.)) +
  # Reference temperature points (gray, hollow)
  geom_point(
    data = ref_temp_points,
    inherit.aes = FALSE,
    aes(x = Temp..c., y = Depth..m.),
    color = "gray50",
    size = 4,
    shape = 21,
    fill = NA,
    stroke = 1,
    alpha = 0.7
  ) +
  # Metadata sample points: temperature
  geom_point(aes(x = temperature, shape = Sample), fill = "black", size = 5, color = "black", stroke = 1) +
  # Smooth line for temperature
  geom_smooth(aes(x = temperature), method = "lm", se = FALSE, color = "black", size = 1) +
  # Manual sample shapes
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  # X axis: temp only (top)
  scale_x_continuous(
    name = "Temperature (°C)",
    position = "top",
    limits = c(13.6, NA)
    # You can set limits here if needed: limits = c(min_temp, max_temp)
  ) +
  # Y axis: depth, reversed, with exact breaks
  scale_y_reverse(
    name = "Depth Below Surface (m)",
    limits = depth_limits_rev,
    breaks = desired_depth_breaks,
    labels = desired_depth_breaks
  ) +
  # Minimal theme and exact aesthetic tweaks
  theme_minimal() +
  theme(
    axis.title.x.top = element_text(color = "black", size = 18, vjust = 0),
    axis.text.x.top  = element_text(color = "black", size = 16),
    axis.title.y     = element_text(size = 18, color = "black"),
    axis.text.y      = element_text(size = 18, color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 18),
    legend.key.size  = grid::unit(1.2, "cm"),
    plot.margin      = margin(5, 5, 5, 5)
  )

print(temp2)
# -------- TDN Plot: TDN vs Depth (y, reversed) --------
# Calculate correlation stats first
TDN_cor <- cor.test(metadata$Depth..m., metadata$TDN, method = "spearman")
TDN_rho <- round(TDN_cor$estimate, 3) #0.886
TDN_p   <- round(TDN_cor$p.value, 3) #0.033

# Get axis limits for annotation placement
x_max <- max(metadata$Depth..m., na.rm = TRUE)
y_min <- min(metadata$TDN, na.rm = TRUE)

# Main plot, no reference layer
TDN2 <- ggplot(metadata, aes(y = Depth..m.)) +
  geom_smooth(aes(x = TDN), method = "lm", se = FALSE, color = "black", size = 1) +
  geom_point(aes(x = TDN, shape = Sample, fill = Sample), size = 5, color = "black", stroke = 1) +
  scale_shape_manual(name = "Sample", values = custom_shapes) +
  scale_fill_manual(name = "Sample", values = custom_fills) +
  scale_x_continuous(name = "TDN (mg/L)", position = "top") +
  scale_y_reverse(
    name = "Depth Below Surface (m)",
    limits = depth_limits_rev,
    breaks = desired_depth_breaks,
    labels = desired_depth_breaks
  ) +
  theme_minimal() +
  theme(
    axis.title.x.top  = element_text(color = "black", size = 18, vjust = 0),
    axis.text.x.top   = element_text(color = "black", size = 16),
    axis.title.y      = element_text(size = 18, color = "black"),
    axis.text.y       = element_text(size = 18, color = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title      = element_text(size = 22),
    legend.text       = element_text(size = 20),
    legend.key.size   = grid::unit(1.2, "cm"),
    plot.margin       = margin(5, 5, 5, 5)
  )

print(TDN2)



#### Water Isotope Analysis ####
# Reference isotope data (hollow gray circles)
ref_iso_points <- reference_geo %>%
  filter(!is.na(d18O.H2O), !is.na(d2H.H2O)) %>%
  mutate(delta18O = d18O.H2O, delta2H = d2H.H2O)

# My 6 sites: "This Study"
my_sites_data <- combined_plot_data %>%
  filter(group == "This Study" & !is.na(delta18O) & !is.na(delta2H))

# MWL line (Craig 1961)
mwl_range <- range(c(my_sites_data$delta18O, ref_iso_points$delta18O), na.rm = TRUE)
mwl_data <- data.frame(
  delta18O = seq(mwl_range[1], mwl_range[2], length.out = 100)
)
mwl_data$delta2H <- 8 * mwl_data$delta18O + 10

# Error bars for my sites (±1 for delta2H)
my_sites_data$delta2H_lower <- my_sites_data$delta2H - 1
my_sites_data$delta2H_upper <- my_sites_data$delta2H + 1

water_isotope <- ggplot() +
  # Reference isotope points (hollow gray)
  geom_point(
    data = ref_iso_points,
    aes(x = delta18O, y = delta2H),
    color = "gray50",
    shape = 21,
    fill = NA,
    alpha = 0.7,
    size = 4,
    stroke = 1
  ) +
  # MWL line (Craig 1961)
  geom_line(
    data = mwl_data,
    aes(x = delta18O, y = delta2H),
    linetype = "solid",
    color = "blue",
    size = 1.2
  ) +
  # Error bars for sites
  geom_errorbar(
    data = my_sites_data,
    aes(x = delta18O, ymin = delta2H_lower, ymax = delta2H_upper),
    width = 0.5,
    color = "black"
  ) +
  #custom shapes, all filled/surrounded in black
  geom_point(
    data = my_sites_data,
    aes(x = delta18O, y = delta2H, shape = Sample),
    color = "black",   # border/line
    fill = "black",    # inside fill
    size = 5,
    stroke = 1
  ) +
  scale_shape_manual(values = shape_values) + # You keep custom shapes!
  labs(
    x = expression(delta^{18}*O~"(‰ VSMOW)"),
    y = expression(delta^{2}*H~"(‰ VSMOW)")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

water_isotope



#### Carbon Isotope ####
# 1. Prepare reference isotope points
ref_carbon_points <- reference_geo %>%
  filter(!is.na(d13C.CH4), !is.na(d13C.CO2)) %>%
  mutate(d13CH4 = d13C.CH4, d13CO2 = d13C.CO2)

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
  
  # Reference isotope points (hollow gray, behind your samples)
  geom_point(
    data = ref_carbon_points,
    inherit.aes = FALSE,
    aes(x = d13CH4, y = d13CO2),
    color = "gray50",
    shape = 21,
    fill = NA,
    alpha = 0.7,
    size = 4,
    stroke = 1
  ) +
  
  #sample points (custom shapes/fills, solid)
  geom_point(size = 5, color = "black", stroke = 1) +
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
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

carbon_isotope
#The isotope boundaries are based on (Whiticar 1999, Carbon and hydrogen isotope systematics of bacterial formation and oxidation of methane)



# -------- Tapp Plot: Tapp vs Depth (y, reversed) --------
tapp_plot2 <- ggplot(metadata_tapp, aes(y = Depth..m., x = Tapp)) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 0, xmax = 80,
           fill = "forestgreen", alpha = 0.3) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 80, xmax = 120,
           fill = "gray40", alpha = 0.3) +
  geom_linerange(aes(xmin = Tapp_lower, xmax = Tapp_upper),
                 color = "black", width = 0.5) +
  geom_point(aes(shape = Sample, fill = Sample),
             size = 5, color = "black", stroke = 1) +
  scale_shape_manual(values = custom_shapes, guide = "none") +
  scale_fill_manual(values = custom_fills, guide = "none") +
  scale_x_continuous(name = "Apparent Formation Temperature (°C)", position = "top") +
  scale_y_reverse(name = "Depth Below Surface (m)") +
  coord_cartesian(ylim = c(depth_max, depth_min), xlim = c(0, 120)) +
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
tapp_plot2

# -------- NPDOC Plot, y reversed --------
# Calculate correlation stats first
NPDOC_cor <- cor.test(metadata$Depth..m., metadata$NPDOC, method = "spearman")
NPDOC_rho <- round(NPDOC_cor$estimate, 3) #0.143
NPDOC_p   <- round(NPDOC_cor$p.value, 3) #0.803

# Get axis limits for annotation placement
x_max <- max(metadata$Depth..m., na.rm = TRUE)
y_min <- min(metadata$NPDOC, na.rm = TRUE)

# Main plot, no reference layer
NPDOC2 <- ggplot(metadata, aes(y = Depth..m.)) +
  geom_smooth(aes(x = NPDOC), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 1) +
  geom_point(aes(x = NPDOC, shape = Sample, fill = Sample), size = 5, color = "black", stroke = 1) +
  scale_shape_manual(name = "BNG Well", values = custom_shapes) +
  scale_fill_manual(name = "BNG Well", values = custom_fills) +
  scale_x_continuous(name = "Non-Purgeable Dissolved Organic Carbon (mg/L)", position = "top") +
  scale_y_reverse(
    name = "Depth Below Surface (m)",
    limits = depth_limits_rev,
    breaks = desired_depth_breaks,
    labels = desired_depth_breaks
  ) +
  theme_minimal() +
  theme(
    axis.title.x.top  = element_text(color = "black", size = 18, vjust = 0),
    axis.text.x.top   = element_text(color = "black", size = 16),
    axis.title.y      = element_text(size = 18, color = "black"),
    axis.text.y       = element_text(size = 18, color = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title      = element_text(size = 18),
    legend.text       = element_text(size = 16),
    legend.key.size   = grid::unit(1.2, "cm"),
    plot.margin       = margin(5, 5, 5, 5)
  )

print(NPDOC2)


lm_NPDOC <- lm(Depth..m. ~ NPDOC, data = metadata)
r2_NPDOC <- summary(lm_NPDOC)$r.squared #0.08
# -------- CH4 Dissolved Plot, y reversed --------
CH4_cor <- cor.test(metadata$Depth..m., metadata$Dissolved.CH4, method = "spearman")
CH4_rho <- round(CH4_cor$estimate, 3) #-0.80
CH4_p   <- round(CH4_cor$p.value, 3) #0.13

x_max <- max(metadata$Depth..m., na.rm = TRUE)
y_min <- min(metadata$Dissolved.CH4, na.rm = TRUE)

CH4_plot <- ggplot(metadata, aes(y = Depth..m.)) +
  geom_smooth(aes(x = Dissolved.CH4), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 1) +
  geom_point(aes(x = Dissolved.CH4, shape = Sample, fill = Sample), size = 5, color = "black", stroke = 1) +
  scale_shape_manual(name = "BNG Well", values = custom_shapes) +
  scale_fill_manual(name = "BNG Well", values = custom_fills) +
  scale_x_continuous(name = "Dissolved Methane (mg/L)", position = "top") +
  scale_y_reverse(
    name = "Depth Below Surface (m)",
    limits = depth_limits_rev,
    breaks = desired_depth_breaks,
    labels = desired_depth_breaks
  ) +
  theme_minimal() +
  theme(
    axis.title.x.top  = element_text(color = "black", size = 18, vjust = 0),
    axis.text.x.top   = element_text(color = "black", size = 16),
    axis.title.y      = element_text(size = 18, color = "black"),
    axis.text.y       = element_text(size = 18, color = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title      = element_text(size = 18),
    legend.text       = element_text(size = 16),
    legend.key.size   = grid::unit(1.2, "cm"),
    plot.margin       = margin(5, 5, 5, 5)
  )

print(CH4_plot)

lm_CH4<- lm(Depth..m. ~ Dissolved.CH4, data = metadata)
r2_CH4 <- summary(lm_CH4)$r.squared #0.68


#### Combine all inverted plots into a single layout ####
# Remove legends from updated panels (do not overwrite originals)
salinity2_panel             <- salinity2      + theme(legend.position = "none")
pH2_panel                   <- pH2            + theme(legend.position = "none")
TDN2_panel                  <- TDN2           + theme(legend.position = "none", axis.title.y = element_blank())
temp2_panel                 <- temp2          + theme(legend.position = "none", axis.title.y = element_blank())
tapp_plot2_panel            <- tapp_plot2     + theme(legend.position = "none")
carbon_isotope_panel        <- carbon_isotope + theme(legend.position = "none")
CH4_plot_panel             <- CH4_plot       + theme(legend.position = "none")
NPDOC2_panel                <- NPDOC2         + theme(legend.position = "none")
water_isotope_panel          <- water_isotope  + theme(legend.position = "none")
SO4_plot_panel              <- SO4_plot     +  theme(legend.position = "none", axis.title.y = element_blank())

# Extract Sample Site (shape/fill) legend from TDN2
site_legend_panel <- cowplot::get_legend(
  TDN2 +
    guides(shape = guide_legend(ncol = 1)) +
    theme(
      legend.position   = "right",
      legend.title      = element_text(size = 22),
      legend.text       = element_text(size = 18),
      legend.key.width  = unit(1.2, "cm"),
      legend.key.height = unit(1.2, "cm"),
      legend.spacing.y  = unit(0.4, "cm")
    )
)

# Build a dummy plot for Gas Origin legend
gas_origin_legend_plot_panel <- ggplot(data.frame(x = 1, y = 1, Origin = c("Biogenic", "Thermogenic")),
                                       aes(x = x, y = y, fill = Origin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Gas Origin",
    values = c(
      "Biogenic"    = alpha("forestgreen", 0.6),
      "Thermogenic" = alpha("gray40", 0.6)
    )
  ) +
  theme_void() +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 18),
    legend.text        = element_text(size = 16),
    legend.key.width   = unit(1.2, "cm"),
    legend.key.height  = unit(1.2, "cm"),
    legend.spacing.y   = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

gas_legend_panel <- cowplot::get_legend(gas_origin_legend_plot_panel)

# Combine both legends
combined_legends_panel <- plot_grid(site_legend_panel, gas_legend_panel, ncol = 2, rel_heights = c(1, 1))

#Combined geochem and isotope plot
# Assemble final combined plot with new name
combined_plot_6panel_v2 <- (
  (salinity2_panel | pH2_panel) /
    (TDN2_panel      | temp2_panel) /
    (tapp_plot2_panel| carbon_isotope_panel) /
    ((water_isotope  | combined_legends_panel) + patchwork::plot_layout(widths = c(3, 1)))
) +
  patchwork::plot_annotation(
    tag_levels = list(c("(A)  1.", "(A)  2.", "(A)  3.", "(A)  4.", "(B)", "(C)", "(D)"))
  ) &
  theme(plot.tag = element_text(size = 25, face = "bold"))

combined_plot_6panel_v2 <- combined_plot_6panel_v2 + patchwork::plot_layout(heights = c(1, 1, 1, 1))

# Show the plot
print(combined_plot_6panel_v2)

# Save as PDF
ggsave("combined_geochem3.pdf",
       plot = combined_plot_6panel_v2,
       width = 18, height = 24, units = "in", device = "pdf", limitsize = FALSE)


#Plot geochem and isotopes seperately
# Build Plot 1 grid
plot1_geochem <- (
  (salinity2_panel | pH2_panel) /
    (temp2_panel     | TDN2_panel) /
    (NPDOC2_panel    | CH4_plot_panel) /
    site_legend_panel
) +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 25, face = "bold"))

# Print & Save Plot 1
print(plot1_geochem)


ggsave("plot1_geochem.pdf", plot = plot1_geochem, width = 16, height = 22, units = "in", device = "pdf", limitsize = FALSE)


#isotope plot
# Build Plot 2 grid
plot2_isotopes <- (
  # Top row: equal widths for both panels
  (carbon_isotope_panel | tapp_plot2_panel) / 
    # Bottom row: water isotope is 2/3, legends is 1/3
    ((water_isotope_panel | combined_legends_panel) | patchwork::plot_layout(widths = c(2, 1)))
) +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 25, face = "bold"))

# Print & Save Plot 2
print(plot2_isotopes)
ggsave("plot2_isotopes.pdf", plot = plot2_isotopes, width = 12, height = 10, units = "in", device = "pdf", limitsize = FALSE)




#Plot sulfate vs depth


# -------- SO4 Supplementary Plot: SO4 vs Depth (y, reversed) --------
SO_cor <- cor.test(metadata$Depth..m., metadata$Sulfate.mg.L., method = "spearman")
SO_rho <- round(SO_cor$estimate, 3) #0.203
SO_p <- signif(SO_cor$p.value, 3) #0.7

lm_SO<- lm(Depth..m. ~ Sulfate.mg.L., data = metadata)
r2_SO <- summary(lm_SO)$r.squared #0.20

SO4_plot <- ggplot(metadata, aes(y = Depth..m.)) +
  # LM regression line (solid black line)
  geom_smooth(aes(x = Sulfate.mg.L.), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 1) +
  geom_point(
    aes(x = Sulfate.mg.L., shape = Sample, fill = Sample), 
    size = 5, color = "black", stroke = 1
  ) +
  scale_shape_manual(name = "BNG Well", values = custom_shapes) +
  scale_fill_manual(name = "BNG Well", values = custom_fills) +
  scale_x_continuous(name = "Sulfate (mg/L)", position = "top") +
  scale_y_reverse(
    name = "Depth Below Surface (m)",
    limits = depth_limits_rev,
    breaks = desired_depth_breaks,
    labels = desired_depth_breaks
  ) +
  theme_minimal() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    axis.title.x.top  = element_text(size = 18, color = "black", vjust = 0),
    axis.text.x.top   = element_text(size = 16, color = "black"),
    axis.title.y      = element_text(size = 18, color = "black"),
    axis.text.y       = element_text(size = 18, color = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_blank(),
    legend.title      = element_text(size = 18),
    legend.text       = element_text(size = 16),
    legend.key.size   = grid::unit(1.2, "cm"),
    plot.margin       = margin(5, 5, 5, 5)
  )

print(SO4_plot)

# -------- combine supplemental plot --------
# Assemble final combined plot with new name
combined_plot_6panel_supp <- (
  (pH2_panel | temp2_panel) /
    (NPDOC2_panel | TDN2_panel) /
    (CH4_plot_panel | SO4_plot_panel)
) +
  patchwork::plot_annotation(
    tag_levels = list(c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"))
  ) &
  theme(plot.tag = element_text(size = 25, face = "bold"))

print(combined_plot_6panel_supp)

# Save as PDF
ggsave("combined_geochem_supplemental.pdf",
       plot = combined_plot_6panel_supp,
       width = 16, height = 16, units = "in", device = "pdf", limitsize = FALSE)


# -------- maintext geochem plot  --------
# Compose the lower row, setting widths only there:
lower_row <- (carbon_isotope_panel | tapp_plot2_panel | gas_legend_panel) + 
  plot_layout(widths = c(1, 1, 0.3))

# Combine all together:
lower_row <- (carbon_isotope_panel | tapp_plot2_panel | gas_legend_panel) + 
  plot_layout(widths = c(1, 1, 0.3))

# Combine all together:
combined_plot_4panel_maintext <- (
  ((salinity2_panel | water_isotope_panel) / lower_row) + site_legend_panel
) +
  patchwork::plot_annotation(
    tag_levels =  list(c("(A)", "(B)", "(C)", "(D)"))
  ) &
  theme(plot.tag = element_text(size = 25, face = "bold"))


combined_plot_4panel_maintext

# Save as PDF

ggsave("combined_geochem_maintext.pdf",
       plot = combined_plot_4panel_maintext,
       width = 18, height = 18, units = "in", device = "pdf", limitsize = FALSE)
