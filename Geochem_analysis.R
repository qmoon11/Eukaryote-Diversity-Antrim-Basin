# Thisn script processes the chemical and geochemical analyses of the Antrim gradient
#### Plot Geochemistry vs Depth ####

# Fit a linear regression model for TDS vs Depth
lm_model <- lm(TDS ~ Depth, data = metadata_ITS)

# View the summary of the linear model
summary(lm_model)

# Run Spearman correlation between Depth and TDS
spearman_test <- cor.test(metadata_ITS$Depth, metadata_ITS$TDS, method = "spearman")

# Extract values for annotation
rho_val <- round(spearman_test$estimate, 2)
p_val <- signif(spearman_test$p.value, 3)
r_squared <- summary(lm_model)$r.squared

# Plot TDS vs Depth with regression line and annotation
TDS_plot <- ggplot(metadata_ITS, aes(x = Depth, y = TDS)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  annotate("text", 
           x = 500,    # position text slightly left of max depth
           y = 173000, # position text near top of TDS range
           label = paste0("Spearman's \u03C1 = ", round(rho_val, 2),
                          "\n", "p = ", round(p_val, 3),
                          "\n", "R² = ", round(r_squared, 2)),
           size = 4, hjust = 0) +
  labs(y = "TDS (mg/L)") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750), limits = c(800, 1850)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Plot Temperature vs Depth
temp_plot <- ggplot(metadata_ITS, aes(x = Depth, y = temperature)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Temperature (°C)") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Plot pH vs Depth
pH_plot <- ggplot(metadata_ITS, aes(x = Depth, y = pH)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "pH") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Plot Ammonium vs Depth
Ammonium_plot <- ggplot(metadata_ITS, aes(x = Depth, y = Ammonium)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Ammonium (mg/L)") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank()
  )

# Plot NPDOC vs Depth
NPDOC_plot <- ggplot(metadata_ITS, aes(x = Depth, y = NPDOC)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "NPDOC (mg/L)", x = "Sample Depth (m)") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Plot Alkalinity vs Depth
Alk_plot <- ggplot(metadata_ITS, aes(x = Depth, y = alkinity)) +
  geom_point(color = "steelblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(y = "Alkalinity (mg/L)", x = "Sample Depth (m)") +
  scale_x_continuous(breaks = c(1000, 1250, 1500, 1750)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Combine all plots into a 2-column, 3-row layout with tags (A), (B), ...
combined_chem_plot <- (TDS_plot | temp_plot) /
  (pH_plot | Ammonium_plot) /
  (NPDOC_plot | Alk_plot) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

# Save combined plot as SVG
ggsave("combined_chem_plot.svg",
       plot = combined_chem_plot,
       width = 15, height = 15, units = "in", device = "svg", limitsize = FALSE)
