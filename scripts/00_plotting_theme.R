# 00_plotting_theme.R
# Author: Pramod BR
# Date:   2026-07-24
# Shared plotting style so every figure in the project looks consistent.
# Sourced by 02, 06 and 07.

library(ggplot2)

# --- Palette ----------------------------------------------------------------
pk_col   <- "#2C6E9B"   # concentration (blue)
pd_col   <- "#B5453B"   # effect / PCA (red)
accent   <- "#E4A700"   # highlights
ink      <- "#1A1A1A"   # text / lines

# --- Theme ------------------------------------------------------------------
theme_pk <- function(base_size = 13) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 3,
                                     hjust = 0, margin = margin(b = 4)),
      plot.subtitle   = element_text(colour = "grey35", size = base_size - 1,
                                     hjust = 0, margin = margin(b = 10)),
      plot.caption    = element_text(colour = "grey55", size = base_size - 3,
                                     hjust = 1, margin = margin(t = 8)),
      axis.title      = element_text(colour = "grey20", face = "bold"),
      axis.text       = element_text(colour = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
      strip.text      = element_text(face = "bold", size = base_size,
                                     margin = margin(4, 4, 4, 4)),
      strip.background = element_rect(fill = "grey95", colour = NA),
      plot.margin     = margin(14, 16, 12, 14),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "top",
      legend.title    = element_text(face = "bold", size = base_size - 2),
      legend.text     = element_text(size = base_size - 2)
    )
}

theme_set(theme_pk())

# --- Helper: consistent, high-resolution export -----------------------------
save_fig <- function(plot, file, width = 7.5, height = 5) {
  ggsave(file, plot, width = width, height = height, dpi = 300,
         bg = "white")
}
