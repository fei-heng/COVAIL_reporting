library(dplyr)
library(R.matlab)
library(R.utils)
library(ggplot2)
library(ggpubr)
library(GGally)
library(gridExtra)
library(glue)
library(cowplot)
library(grid)

#### data preparation ####

df_marks_all <- read.csv("covail_rbd_dms_escape_marks_v6.csv")

df <- read.csv("covail_data_processed_20250524.csv")

non_case_ids <- df %>%
  filter(EventIndPrimaryD15 == 0) %>%
  pull(Ptid)

# Remove non-cases from df_marks
df_marks <- df_marks_all %>%
  filter(swab.num == "SWAB1") %>%
  filter(!ptid %in% non_case_ids)

df <- df %>%
  left_join(df_marks, by = c("Ptid" = "ptid")) %>%
  filter(!arm %in% c(3, 10, 11, 13, 14, 15),
         TwophasesampIndD15 == 1,
         Immunemarkerset == 1)

# define immune marker
# Day 15 nAb ID50 levels against (D614G, BA.1, BA.4/5, XBB.1.5) 
df$IM <- df$Day15pseudoneutid50_D614G # Prototype
df$IM[df$TrtB==0] <- df$Day15pseudoneutid50_BA.1[df$TrtB==0] # Omicron
df$IM[df$arm==17] <- df$Day15pseudoneutid50_BA.4.BA.5[df$arm==17] # Omicron - arm 17

### rename distance marks for Prototype and Omicron groups

df <- df %>% 
  mutate(mark.rbd.rescaled = ifelse(arm %in% c(1, 7, 10, 11),
                                    wuhan.rbd.rescaled, 
                                    hybrid.vx.rbd.rescaled),
         mark.rbd.denoised.rescaled = ifelse(arm %in% c(1, 7, 10, 11),
                                             wuhan.rbd.denoised.rescaled,
                                             hybrid.vx.rbd.denoised.rescaled),
         mark.dadonaite.rescaled = hybrid.vx.dadonaite.rescaled)

df <- df %>% 
  mutate(mark.rbd= ifelse(arm %in% c(1, 7, 10, 11),
                          wuhan.rbd, 
                          hybrid.vx.rbd),
         mark.rbd.denoised = ifelse(arm %in% c(1, 7, 10, 11),
                                    wuhan.rbd.denoised,
                                    hybrid.vx.rbd.denoised),
         mark.dadonaite = hybrid.vx.rbd.dadonaite)

# rescale rule verification
# wuhan.rbd/max(wuhan.rbd, na.rm = T)
# wuhan.rbd.denoised/max(wuhan.rbd.denoised, na.rm = T)
# hybrid.vx.rbd/max(hybrid.vx.rbd, na.rm = T)
# hybrid.vx.rbd.denoised/max(hybrid.vx.rbd.denoised, na.rm = T)
# hybrid.vx.rbd.dadonaite/max(hybrid.vx.rbd.dadonaite, na.rm = T)

rescale_prototype <- c(1,1,max(df_marks_all$wuhan.rbd, na.rm = T),
                       max(df_marks_all$wuhan.rbd.denoised, na.rm = T),
                       max(df_marks_all$hybrid.vx.rbd.dadonaite, na.rm = T))

rescale_omicron <- c(1,1,
                     max(df_marks_all$hybrid.vx.rbd, na.rm = T),
                     max(df_marks_all$hybrid.vx.rbd.denoised, na.rm = T),
                     max(df_marks_all$hybrid.vx.rbd.dadonaite, na.rm = T))

# create vaccine type indicator
df$vac_type <- 'Pfizer'
ind_Moderna <- df$arm %in% c(1, 2, 4, 5, 6)
df$vac_type[ind_Moderna] <- 'Moderna'

#### Table 1 ####

# GM (95% CI)
# logx1, ..., logxn
# GM = exp(mean(logxi))
# LL = exp(mean(logxi) - t0.975,n-1 * sd(logxi) / sqrt(n))
# UL = exp(mean(logxi) + t0.975,n-1 * sd(logxi) / sqrt(n))

## group by naive (get numbers for naive and non-naive separately)
df_summary <- df %>%
  filter(arm==17) %>% # TrtB==0 (Omicron), 1 (Prototype); arm==1,7,2,4,...,17
  group_by(naive, EventIndPrimaryD15) %>%
  reframe(
    n = sum(!is.na(IM)),
    n_viral = sum(!is.na(seqid)),
    mean_logIM = mean(IM, na.rm = TRUE),
    sd_logIM = sd(IM, na.rm = TRUE),
    t_value = qt(0.975, df = n - 1),
    GM = 10^mean_logIM,
    LL = 10^(mean_logIM - t_value * sd_logIM / sqrt(n)),
    UL = 10^(mean_logIM + t_value * sd_logIM / sqrt(n))
  ) %>%
  select(naive, EventIndPrimaryD15, n, GM, LL, UL, n_viral)
df_summary

df_summary %>%
  mutate(
    summary_text = glue("{n}, {round(GM, 0)} ({round(LL, 0)}, {round(UL, 0)})")
  ) %>%
  select(naive, EventIndPrimaryD15, summary_text) %>%
  print(n = Inf) 


## naive and non-naive combined
df_summary <- df %>%
  filter(arm==17) %>% # TrtB==0, 1; arm==1,7,2,4,...,17
  group_by(EventIndPrimaryD15) %>%
  reframe(
    n = sum(!is.na(IM)),
    n_viral = sum(!is.na(seqid)),
    mean_logIM = mean(IM, na.rm = TRUE),
    sd_logIM = sd(IM, na.rm = TRUE),
    t_value = qt(0.975, df = n - 1),
    GM = 10^mean_logIM,
    LL = 10^(mean_logIM - t_value * sd_logIM / sqrt(n)),
    UL = 10^(mean_logIM + t_value * sd_logIM / sqrt(n))
  ) %>%
  select(EventIndPrimaryD15, n, GM, LL, UL, n_viral)
df_summary

df_summary %>%
  mutate(
    summary_text = glue("{n}, {round(GM, 0)} ({round(LL, 0)}, {round(UL, 0)})")
  ) %>%
  select(EventIndPrimaryD15, summary_text) %>%
  print(n = Inf) 


#### Figure 1: Viral Sequence Distances by Vaccine Type and Prior Infection ####
df_plot <- df %>%
  mutate(
    group = paste(vac_type, ifelse(naive == 1, "Naïve", "Non-naïve")),
    group = factor(group, levels = c("Moderna Naïve", "Moderna Non-naïve", "Pfizer Naïve", "Pfizer Non-naïve"))
  )

# add "All"
df_plot_all <- df_plot %>%
  mutate(group = "All")

df_plot <- bind_rows(df_plot, df_plot_all)
df_plot$group <- factor(df_plot$group, 
                        levels = c("Moderna Naïve", "Moderna Non-naïve", 
                                   "Pfizer Naïve", "Pfizer Non-naïve", "All"))


# lineage
df_plot <- df_plot %>%
  mutate(lineage.label = case_when(
    grepl("^BA\\.2", lineage.coarse) ~ "2",
    grepl("^BA\\.4", lineage.coarse) ~ "4",
    grepl("^BA\\.5", lineage.coarse) ~ "5",
    grepl("^XZ", lineage.coarse) ~ "XZ",
    grepl("^XBB", lineage.coarse) ~ "XB", 
    is.na(lineage.coarse) | lineage.coarse == "" ~ "M",
    TRUE ~ "M"  # fallback
  ))

# with(subset(df, EventIndPrimaryD15 == 1), table(lineage.coarse, TrtB, useNA = "ifany"))
#            TrtB
# lineage    0  1
# BA.2      28  7
# BA.4      11  3
# BA.5      67 22
# XBB.1.5    2  0
# XBB??0.10  1  0
# XZ         1  0
# <NA>      42 11

V_list <- c("zhd.spike", 
            "zhd.spike", 
            "zhd.rbd", 
            "zhd.rbd", 
            "mark.rbd", 
            "mark.rbd.denoised", 
            "mark.dadonaite")

TrtB_list <- c(1, 0, 1, 0, 0, 0, 0)

ylab_list <- c("Spike Hamming Distance",
               "Spike Hamming Distance",
               "RBD Hamming Distance",
               "RBD Hamming Distance",
               "DMS-escape RBD-1",
               "DMS-escape RBD-2",
               "DMS-escape RBD-3")

ylim_list <- list(c(0,40),c(0,40),c(0,20),c(0,20),c(0,4),c(0,3),c(0,3))

color_palette <- c("Moderna Naïve" = "#E69F00", 
                   "Moderna Non-naïve" = "#D55E00",
                   "Pfizer Naïve" = "#56B4E9", 
                   "Pfizer Non-naïve" = "#0072B2",
                   "All" = "gray50")


# not used
get_percent_text <- function(df, var) {
  df %>%
    mutate(positive = .data[[var]] > 0) %>%
    group_by(group) %>%
    summarise(
      n_pos = sum(positive, na.rm = TRUE),
      n_total = sum(!is.na(.data[[var]])),
      pct = round(100 * n_pos / n_total, 1),
      .groups = "drop"
    ) %>%
    mutate(pct_label = paste0(group, ": ", pct, "% (", n_pos, "/", n_total, ")")) %>%
    pull(pct_label) %>%
    paste(collapse = "\n") %>%
    paste("% of distances > 0:", ., sep = "\n")
}


# plot 
create_panel_plot <- function(df, var, title, add_pct = FALSE, ylim_range) {
  p <- ggplot(df, aes(x = group, y = .data[[var]], fill = group, color = group)) +
    geom_violin(alpha = 0.5, width = 0.9, trim = TRUE) +
    #geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.6, size = 1.2, shape = 16) +
    geom_text(
      aes(label = lineage.label),
      position = position_jitter(width = 0.25, height = 0),
      size = 3.0,
      alpha = 0.8
    ) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    # scale_y_continuous(labels = function(y) {
    #   if (ylim_range[2] < 5) {
    #     sprintf("%.2f", y)
    #   } else {
    #     as.character(y)
    #   }
    # }) +
    labs(x = NULL, y = title) +
    coord_cartesian(ylim = ylim_range) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  if (add_pct == 1) {
    percent_text <- get_percent_text(df, var)
    y_max <- max(df[[var]], na.rm = TRUE)
    p <- p + annotate("text", x = 5, y = 3, label = percent_text,
                      hjust = 1, vjust = 1, size = 3.5, fontface = "italic")
  } else if (add_pct == 2){
    percent_text <- get_percent_text(df, var)
    y_max <- max(df[[var]], na.rm = TRUE)
    p <- p + annotate("text", x = 5, y = 3, label = percent_text,
                      hjust = 1, vjust = 1, size = 3.5, fontface = "italic")
    
  }
  
  return(p)
}

plots <- lapply(seq_along(V_list), function(i) {
  create_panel_plot(df_plot %>% filter(TrtB == TrtB_list[i]), 
                    V_list[i], 
                    ylab_list[i], 
                    F, 
                    ylim_list[[i]])
})

# arrange plots with legend as the 8th panel
legend_plot <- create_panel_plot(df_plot %>% filter(TrtB == TrtB_list[2]),
                                 V_list[1], ylab_list[1], F, ylim_list[[1]]) +
  theme(legend.position = "right")
library(cowplot)
legend <- get_legend(legend_plot)
legend_grob <- as_ggplot(legend)

all_plots <- c(plots[1],
               plots[2],
               plots[3],
               plots[4],
               plots[5],
               plots[6],
               plots[7],
               list(legend_grob))

panel_labels <- c("A - Prototype",
                  "B - Omicron",
                  "C - Prototype",
                  "D - Omicron",
                  "E - Omicron",
                  "F - Omicron",
                  "G - Omicron",
                  "")

combined_plot <- ggarrange(plotlist = all_plots,
                           ncol = 2, nrow = 4,
                           labels = panel_labels,
                           font.label = list(size = 12, face = "bold"),
                           common.legend = F, 
                           legend = NULL)


final_plot <- annotate_figure(
  combined_plot,
  top = text_grob("Figure 1: Viral Sequence Distances by Vaccine Type and Prior Infection \n",
                  face = "bold", size = 12)
)

ggsave("Figure1_violin_7panels.pdf", final_plot, width = 12, height = 16)



#### Figure 2: Pairwise Correlations of Viral Sequence Distances ####

df_plot <- df %>%
  mutate(
    group = paste(vac_type, ifelse(naive == 1, "Naïve", "Non-naïve")),
    group = factor(group, levels = c("Moderna Naïve", "Moderna Non-naïve", 
                                     "Pfizer Naïve", "Pfizer Non-naïve"))
  )

V_list1 <- c("zhd.spike", 
             "zhd.rbd")
V_list2 <- c("zhd.spike", 
            "zhd.rbd", 
            "mark.rbd", 
            "mark.rbd.denoised", 
            "mark.dadonaite")

ylab_list1 <- c("Spike Hamming Distance",
                "RBD Hamming Distance")
ylab_list2 <- c("Spike Hamming Distance",
               "RBD Hamming Distance",
               "DMS-escape RBD-1",
               "DMS-escape RBD-2",
               "DMS-escape RBD-3")

color_palette <- c("Moderna Naïve" = "#E69F00", 
                   "Moderna Non-naïve" = "#D55E00",
                   "Pfizer Naïve" = "#56B4E9", 
                   "Pfizer Non-naïve" = "#0072B2")

shape_palette <- c("Moderna Naïve" = 15, 
                   "Moderna Non-naïve" = 16,
                   "Pfizer Naïve" = 15, 
                   "Pfizer Non-naïve" = 16)

# Spearman rho + 95% CI
cor_ci_text <- function(data, mapping, method = "spearman", R = 1000, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  rho_obs <- cor(x, y, method = method, use = "complete.obs")
  n <- length(x)
  
  boot_rho <- replicate(R, {
    idx <- sample(1:n, replace = TRUE)
    cor(x[idx], y[idx], method = method, use = "complete.obs")
  })
  
  ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)
  label <- sprintf("rho = %.2f\n[%.2f, %.2f]", rho_obs, ci[1], ci[2])
  
  ggally_text(label, ...)
}

# scatter + linear fit
lower_with_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.6, size = 1.5, ...) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", size = 0.3, ...)
}

# plot
create_pairs_plot <- function(df_sub, title, V_list, ylab_list) {
  ggpairs(
    df_sub,
    columns = V_list,
    columnLabels = ylab_list,
    mapping = aes(color = group, shape = group),
    upper = list(continuous = GGally::wrap(cor_ci_text, size = 3)),
    lower = list(continuous = lower_with_lm),
    diag = list(continuous = "blank"),
    title = title
  ) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
}

p1 <- create_pairs_plot(df_plot %>% filter(TrtB == 1), 
                        "Panel A: Prototype Vaccine Group", V_list1, ylab_list1)
p2 <- create_pairs_plot(df_plot %>% filter(TrtB == 0), 
                        "Panel B: Omicron Vaccine Group", V_list2, ylab_list2)

legend_df <- data.frame(
  group = factor(c("Moderna Naïve", "Moderna Non-naïve", 
                   "Pfizer Naïve", "Pfizer Non-naïve"),
                 levels = c("Moderna Naïve", "Moderna Non-naïve", 
                            "Pfizer Naïve", "Pfizer Non-naïve")),
  x = 1, y = 1
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = group, shape = group)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette, drop = FALSE) +
  scale_shape_manual(values = shape_palette, drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

g <- ggplotGrob(legend_plot)
legend_grob <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]


pdf("Figure2_pairs_rho.pdf", width = 12, height = 16)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 5, heights = unit(c(1, 4, 10), "null"))))

# Row 1: title
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:5))
grid.text("Figure 2: Pairwise Correlations of Viral Sequence Distances",
          gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# Row 2: Panel A and legend
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(p1 + theme(legend.position = "none"), newpage = FALSE)  # legend removed
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3:5))
grid.draw(legend_grob)
popViewport()

# Row 3: Panel B
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:5))
print(p2 + theme(legend.position = "none"), newpage = FALSE)
popViewport()

dev.off()

#### Figure 3: Viral Distances vs Day 15 nAb ID50 Titers ####

V_list <- c("zhd.spike", 
            "zhd.spike", 
            "zhd.rbd", 
            "zhd.rbd", 
            "mark.rbd", 
            "mark.rbd.denoised", 
            "mark.dadonaite")

TrtB_list <- c(1, 0, 1, 0, 0, 0, 0)

ylab_list <- c("Spike Hamming Distance",
               "Spike Hamming Distance",
               "RBD Hamming Distance",
               "RBD Hamming Distance",
               "DMS-escape RBD-1",
               "DMS-escape RBD-2",
               "DMS-escape RBD-3")

xlim_list <- list(c(3.5, 5),c(2,5),c(3.5, 5),c(2,5),c(2,5),c(2,5),c(2,5))
ylim_list <- list(c(0,40),c(0,40),c(0,20),c(0,20),c(0,4),c(0,3),c(0,3))


color_palette <- c("Moderna Naïve" = "#E69F00", 
                   "Moderna Non-naïve" = "#D55E00",
                   "Pfizer Naïve" = "#56B4E9", 
                   "Pfizer Non-naïve" = "#0072B2")

shape_palette <- c("Moderna Naïve" = 15, 
                   "Moderna Non-naïve" = 16,
                   "Pfizer Naïve" = 15, 
                   "Pfizer Non-naïve" = 16)

# number of unique values
get_unique_text <- function(df, var) {
  n_unique <- length(unique(na.omit(df[[var]])))
  paste("Unique values:", n_unique)
}

# add rho
get_rho_ci_text <- function(x, y, R = 1000) {
  rho_obs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- length(x)
  
  boot_rho <- replicate(R, {
    idx <- sample(1:n, replace = TRUE)
    cor(x[idx], y[idx], method = "spearman", use = "complete.obs")
  })
  
  ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)
  sprintf("rho = %.2f\n[%.2f, %.2f]", rho_obs, ci[1], ci[2])
}

# plot
create_scatter_plot <- function(df, var, title, annotate_unique = FALSE,
                                xlim_range, ylim_range, annotate_spearman = TRUE) {
  p <- ggplot(df, aes(x = IM, y = .data[[var]], color = group, shape = group)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(aes(group = 1), method = "loess", se = FALSE, size = 0.8, color = "black") +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    labs(x = "Day 15 nAb ID50 titer", y = title, title = "") +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
          plot.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  if (annotate_unique) {
    txt <- get_unique_text(df, var)
    p <- p + annotate("text", x = 5, y = max(ylim_range), label = txt,
                      hjust = 1, vjust = 1, size = 3.5, fontface = "italic")
  }
  
  if (annotate_spearman) {
    rho_text <- get_rho_ci_text(df$IM, df[[var]])
    p <- p + annotate("text", x = min(xlim_range), y = max(ylim_range),
                      label = rho_text, hjust = 0, vjust = 1, size = 3.5, fontface = "bold")
  }
  
  return(p)
}

df_plot <- df %>%
  filter(EventIndPrimaryD15 == 1) %>%
  mutate(
    group = paste(vac_type, ifelse(naive == 1, "Naïve", "Non-naïve")),
    group = factor(group, levels = names(color_palette))
  )

plots <- lapply(seq_along(V_list), function(i) {
  create_scatter_plot(df_plot %>% filter(TrtB == TrtB_list[i]), 
                      V_list[i], 
                      ylab_list[i], 
                      annotate_unique = TRUE, 
                      xlim_range = xlim_list[[i]],
                      ylim_range = ylim_list[[i]],
                      annotate_spearman = TRUE)
})

# arrange plots with legend as the 8th panel
legend_plot <- create_scatter_plot(df_plot %>% filter(TrtB == TrtB_list[2]),
                                   V_list[2], ylab_list[2], 1, xlim_list[[2]],ylim_list[[2]]) +
  theme(legend.position = "right")
legend <- get_legend(legend_plot)
legend_grob <- as_ggplot(legend)

all_plots <- c(plots[1],
               plots[2],
               plots[3],
               plots[4],
               plots[5],
               plots[6],
               plots[7],
               list(legend_grob))

panel_labels <- c("A - Prototype",
                  "B - Omicron",
                  "C - Prototype",
                  "D - Omicron",
                  "E - Omicron",
                  "F - Omicron",
                  "G - Omicron",
                  "")

combined_plot <- ggarrange(plotlist = all_plots,
                           ncol = 2, nrow = 4,
                           labels = panel_labels,
                           font.label = list(size = 12, face = "bold"),
                           common.legend = F, 
                           legend = NULL)

final_plot <- annotate_figure(combined_plot,
                              top = text_grob("Figure 3: Viral Distances vs Day 15 nAb ID50 Titers\n",
                                              face = "bold", size = 16))

ggsave("Figure3_mark_vs_IM_with_rho.pdf", final_plot, width = 12, height = 16)

#### Figure S1 ####
color_palette <- c("Moderna" = "#E69F00", "Pfizer" = "#56B4E9")
shape_palette <- c("Moderna" = 15, "Pfizer" = 15)

df_plot <- df %>%
  filter(EventIndPrimaryD15 == 1,
         naive == 1) %>%  
  mutate(
    group = vac_type, 
    group = factor(group, levels = c("Moderna", "Pfizer"))
  )


plots <- lapply(seq_along(V_list), function(i) {
  create_scatter_plot(df_plot %>% filter(TrtB == TrtB_list[i]), 
                      V_list[i], 
                      ylab_list[i], 
                      1, 
                      xlim_list[[i]],
                      ylim_list[[i]])
})

# arrange plots with legend as the 8th panel
legend_plot <- create_scatter_plot(df_plot %>% filter(TrtB == TrtB_list[2]),
                                   V_list[2], ylab_list[2], 1, xlim_list[[2]],ylim_list[[2]]) +
  theme(legend.position = "right")
legend <- get_legend(legend_plot)
legend_grob <- as_ggplot(legend)

all_plots <- c(plots[1],
               plots[2],
               plots[3],
               plots[4],
               plots[5],
               plots[6],
               plots[7],
               list(legend_grob))

panel_labels <- c("A - Prototype",
                  "B - Omicron",
                  "C - Prototype",
                  "D - Omicron",
                  "E - Omicron",
                  "F - Omicron",
                  "G - Omicron",
                  "")

combined_plot <- ggarrange(plotlist = all_plots,
                           ncol = 2, nrow = 4,
                           labels = panel_labels,
                           font.label = list(size = 12, face = "bold"),
                           common.legend = F, 
                           legend = NULL)

final_plot <- annotate_figure(combined_plot,
                              top = text_grob("Supplementary Figure 1: Viral Distances vs Day 15 nAb ID50 Titers in Naïve COVID-19 Cases\n",
                                              face = "bold", size = 16))

ggsave("FigureS1_mark_vs_IM_naive.pdf", final_plot, width = 12, height = 16)


#### Figure 4 ####

pdf("Figure4_HR_panels.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 4, 0), mar = c(5.1, 5.1, 4, 5))

panel_titles <- c("A - Prototype", "B - Omicron", "C - Omicron", "D - Omicron")

# p-value 
pval_texts <- c(
  "p-values    \nT1m: 0.0586 \nT2m: 0.3796",
  "p-values    \nT1m: 0.0255 \nT2m: 0.5075",
  "p-values    \nT1m: 0.0017 \nT2m: 0.1931",
  "p-values    \nT1m: <0.0001\nT2m: 0.4746"
)

V_list <- list(
  "zhd.spike",                # A analysis 2.1
  "zhd.spike",                # B analysis 1.1
  "zhd.rbd",                  # C analysis 1.1
  "hybrid.vx.rbd.rescaled"   # D analysis 1.1
)

lab_list <- c("Spike Hamming Distance",
              "Spike Hamming Distance",
              "RBD Hamming Distance",
              "DMS-escape RBD-1")

mean_list <- c(30.7615, 24.6243, 8.3480, 0.3141*max(df_marks_all$hybrid.vx.rbd, na.rm = T))
sd_list <- c(7.1080, 4.0990, 2.6380, 0.1960*max(df_marks_all$hybrid.vx.rbd, na.rm = T))

mat_files <- list(
  "res/a2_1_V1_L5M10_naive1.mat",   # A
  "res/a1_1_V1_L5M10_naive1.mat",   # B
  "res/a1_1_V2_L5M10_naive1.mat",   # C
  "res/a1_1_V3_L5M10_naive1.mat"    # D
)

df_list <- list(
  df %>% filter(TrtB == 1),  # A
  df %>% filter(TrtB == 0),  # B
  df %>% filter(TrtB == 0),  # C
  df %>% filter(TrtB == 0)   # D
)

# 
for (iv in 1:4) {
  var <- V_list[[iv]]
  df_sub <- df_list[[iv]]
  matfilename <- mat_files[[iv]]
  pval_text <- pval_texts[iv]
  
  result <- readMat(matfilename)
  
  xlim_u <- 1
  xlim_l <- 0
  xlim_v <- seq(xlim_l,xlim_u,by=0.1)
  ylim_u <- 5
  ylim_l <- -5
  ylim_b <- c(-6,6)
  ax_b <- c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
  
  i <- 1
  
  v=result$v
  lv=length(v)
  eff=result$eff.beta.hat[((i-1)*lv+1):(i*lv)]
  ratio=exp(eff)
  sig.eff=result$sig.eff.beta[((i-1)*lv+1):(i*lv)]
  v=v[sig.eff[1:lv]!=0]
  lv=length(v)
  
  sv<-c(v,rev(v))
  effEst<-eff[1:lv]
  effU<-eff[1:lv]+1.96*sig.eff[1:lv]
  effL<-eff[1:lv]-1.96*sig.eff[1:lv]
  
  ymax <- ceiling(max(c(abs(effU), abs(effL))))
  
  
  histogram_color <- adjustcolor("#bddcbd", alpha.f = 0.9)
  
  temp <- na.omit(df_sub[, var])
  mean <- mean(temp)
  sd <- sd(temp)
  temp <- pnorm(temp, mean=mean, sd=sd)
  h = hist(temp, plot=F, breaks = seq(0,1,0.1))
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,xlab="", ylab="", ylim=c(0,100),
       axes=FALSE, col=histogram_color, main="", border="#bddcbd")
  
  
  ax <- seq(0,100,25)
  axis(4,las=1,at=ax,labels=paste0(ax,"%"), cex.axis=1)
  mtext(side=4, line=3,  "Percentage", cex=0.9)
  
  par(new=TRUE)
  plot(v,effEst, xlim=c(xlim_l,xlim_u),ylim=c(-ymax,ymax), xlab=" ", ylab=" ",
       lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,cex.sub=1.1,cex.main=1.5, axes=F, type="n")
  
  polygon(sv,c(effL,rev(effU)),col=adjustcolor("#2babd5", alpha.f = 0.10),border=FALSE)
  lines(v,effU,col="#2babd5",lty=3,lwd=1)
  lines(v,effL,col="#2babd5",lty=3,lwd=1)
  lines(v,effEst,col="#2babd5",lty=1,lwd=1.5)
  
  if (iv==1){
    # Define desired log(HR) positions
    yticks <- log(c(0.1, 0.5, 1, 2, 3, 5))
    # Corresponding HR labels
    ytick_labels <- c("0.1", "0.5", "1", "2", "3", "5")
  } else {
    # Define desired log(HR) positions
    yticks <- log(c(0.5, 0.7, 1, 1.5, 2))
    # Corresponding HR labels
    ytick_labels <- c("0.5", "0.7", "1", "1.5", "2")
  }
  
  axis(2, at = yticks, labels = ytick_labels)
  
  #mtext(side=2, line=2.5, expression(paste(beta[2],'(v)')),cex=1.0)
  mtext(side=2, line=2.5, "HR per unit change in Day 15 nAb ID50 titer",cex=0.9)
  
  #axis(1, xlim_v, cex.axis=1)
  if (iv>0){
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv],1), 
         cex.axis=0.9)
  } else{
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv]), 
         cex.axis=0.9)
  }
  
  # mtext(side=1, line=2.5, 'AA sequence distance v', cex=1.5)
  mtext(side=1, line=2.2, lab_list[iv], cex=0.9)
  
  
  abline(h=0,lty=1,lwd=1, col = "gray")
  box()

  mtext(panel_titles[iv], side = 3, line = 1.5, adj = 0, font = 2, cex = 1)
  text(x = 0.35, y = ymax*2/3, labels = pval_text, cex = 1, pos = 3)
}

mtext("Figure 4: Hazard Ratios for Day 15 nAb ID50 Titers by Viral Distance",
      outer = TRUE, cex = 1.5, font = 2)

dev.off()

#### Figure S2 ####

pdf("FigureS2_HR_panels_naive.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 4, 0), mar = c(5.1, 5.1, 4, 5))

panel_titles <- c("A - Prototype", "B - Omicron", "C - Omicron", "D - Omicron")

# p-value 
pval_texts <- c(
  "p-values    \nT1m: 0.1542 \nT2m: 0.3041",
  "p-values    \nT1m: 0.0682 \nT2m: 0.5469",
  "p-values    \nT1m: 0.0270 \nT2m: 0.1684",
  "p-values    \nT1m: 0.0064 \nT2m: 0.2820"
)

V_list <- list(
  "zhd.spike",                # A
  "zhd.spike",                # B
  "zhd.rbd",                  # C
  "hybrid.vx.rbd.rescaled"   # D
)

lab_list <- c("Spike Hamming Distance",
              "Spike Hamming Distance",
              "RBD Hamming Distance",
              "DMS-escape RBD-1")

mean_list <- c(30.8502, 24.6977, 8.2338, 0.3110*max(df_marks_all$hybrid.vx.rbd, na.rm = T))
sd_list <- c(7.2396, 3.7868, 2.5885, 0.1884*max(df_marks_all$hybrid.vx.rbd, na.rm = T))

mat_files <- list(
  "res/a2_3_V1_L5M10.mat",   # A
  "res/a1_3_V1_L5M10.mat",   # B
  "res/a1_3_V2_L5M10.mat",   # C
  "res/a1_3_V3_L5M10.mat"    # D
)

df_list <- list(
  df %>% filter(TrtB == 1, naive==1),  # A
  df %>% filter(TrtB == 0, naive==1),  # B
  df %>% filter(TrtB == 0, naive==1),  # C
  df %>% filter(TrtB == 0, naive==1)   # D
)


for (iv in 1:4) {
  var <- V_list[[iv]]
  df_sub <- df_list[[iv]]
  matfilename <- mat_files[[iv]]
  pval_text <- pval_texts[iv]
  
  result <- readMat(matfilename)
  
  xlim_u <- 1
  xlim_l <- 0
  xlim_v <- seq(xlim_l,xlim_u,by=0.1)
  ylim_u <- 5
  ylim_l <- -5
  ylim_b <- c(-6,6)
  ax_b <- c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
  
  i <- 1
  
  v=result$v
  lv=length(v)
  eff=result$eff.beta.hat[((i-1)*lv+1):(i*lv)]
  ratio=exp(eff)
  sig.eff=result$sig.eff.beta[((i-1)*lv+1):(i*lv)]
  v=v[sig.eff[1:lv]!=0]
  lv=length(v)
  
  sv<-c(v,rev(v))
  effEst<-eff[1:lv]
  effU<-eff[1:lv]+1.96*sig.eff[1:lv]
  effL<-eff[1:lv]-1.96*sig.eff[1:lv]
  
  ymax <- ceiling(max(c(abs(effU), abs(effL))))
  
  
  histogram_color <- adjustcolor("#bddcbd", alpha.f = 0.9)
  
  temp <- na.omit(df_sub[, var])
  mean <- mean(temp)
  sd <- sd(temp)
  temp <- pnorm(temp, mean=mean, sd=sd)
  h = hist(temp, plot=F, breaks = seq(0,1,0.1))
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,xlab="", ylab="", ylim=c(0,100),
       axes=FALSE, col=histogram_color, main="", border="#bddcbd")
  
  
  ax <- seq(0,100,25)
  axis(4,las=1,at=ax,labels=paste0(ax,"%"), cex.axis=1)
  mtext(side=4, line=3,  "Percentage", cex=0.9)
  
  par(new=TRUE)
  plot(v,effEst, xlim=c(xlim_l,xlim_u),ylim=c(-ymax,ymax), xlab=" ", ylab=" ",
       lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,cex.sub=1.1,cex.main=1.5, axes=F, type="n")
  
  polygon(sv,c(effL,rev(effU)),col=adjustcolor("#2babd5", alpha.f = 0.10),border=FALSE)
  lines(v,effU,col="#2babd5",lty=3,lwd=1)
  lines(v,effL,col="#2babd5",lty=3,lwd=1)
  lines(v,effEst,col="#2babd5",lty=1,lwd=1.5)
  
  if (iv==1){
    # Define desired log(HR) positions
    yticks <- log(c(0.1, 0.5, 1, 2, 3, 5))
    # Corresponding HR labels
    ytick_labels <- c("0.1", "0.5", "1", "2", "3", "5")
  } else if (iv==2) {
    # Define desired log(HR) positions
    yticks <- log(c(0.3, 0.5, 1, 2, 3))
    # Corresponding HR labels
    ytick_labels <- c("0.3", "0.5", "1", "2", "3")
  } else {
    # Define desired log(HR) positions
    yticks <- log(c(0.5, 0.7, 1, 1.5, 2))
    # Corresponding HR labels
    ytick_labels <- c("0.5", "0.7", "1", "1.5", "2")
  }
  
  axis(2, at = yticks, labels = ytick_labels)
  
  #mtext(side=2, line=2.5, expression(paste(beta[2],'(v)')),cex=1.0)
  mtext(side=2, line=2.5, "HR per unit change in Day 15 nAb ID50 titer",cex=0.9)
  
  #axis(1, xlim_v, cex.axis=1)
  #axis(1, xlim_v, cex.axis=1)
  if (iv>0){
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv],1), 
         cex.axis=0.9)
  } else{
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv]), 
         cex.axis=0.9)
  }
  
  # mtext(side=1, line=2.5, 'AA sequence distance v', cex=1.5)
  mtext(side=1, line=2.2, lab_list[iv], cex=0.9)
  
  
  abline(h=0,lty=1,lwd=1, col = "gray")
  box()
  
  mtext(panel_titles[iv], side = 3, line = 1.5, adj = 0, font = 2, cex = 1)
  text(x = 0.35, y = ymax*2/3, labels = pval_text, cex = 1, pos = 3)
}

mtext("Supplementary Figure 2: Hazard Ratios for Day 15 nAb ID50 Titers by Viral Distance",
      outer = TRUE, cex = 1.5, font = 2)

dev.off()

#### Figure 5 ####

pdf("Figure5_CIF_panels.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 4, 0), mar = c(5.1, 5.1, 4, 5))

panel_titles <- c("A - Prototype", "B - Omicron", "C - Omicron", "D - Omicron")

V_list <- list(
  "zhd.spike",                # A
  "zhd.spike",                # B
  "zhd.rbd",                  # C
  "hybrid.vx.rbd.rescaled"   # D
)

lab_list <- c("Spike Hamming Distance",
              "Spike Hamming Distance",
              "RBD Hamming Distance",
              "DMS-escape RBD-1")

mean_list <- c(30.7615, 24.6243, 8.3480, 0.3141*max(df_marks_all$hybrid.vx.rbd, na.rm = T))
sd_list <- c(7.1080, 4.0990, 2.6380, 0.1960*max(df_marks_all$hybrid.vx.rbd, na.rm = T))

mat_files <- list(
  "res/a2_1_V1_L5M10_naive1.mat",   # A
  "res/a1_1_V1_L5M10_naive1.mat",   # B
  "res/a1_1_V2_L5M10_naive1.mat",   # C
  "res/a1_1_V3_L5M10_naive1.mat"    # D
)

df_list <- list(
  df %>% filter(TrtB == 1),  # A
  df %>% filter(TrtB == 0),  # B
  df %>% filter(TrtB == 0),  # C
  df %>% filter(TrtB == 0)   # D
)


for (iv in 1:4) {
  var <- V_list[[iv]]
  df_sub <- df_list[[iv]]
  matfilename <- mat_files[[iv]]
  
  result <- readMat(matfilename)
  
  xlim_u <- 1
  xlim_l <- 0
  xlim_v <- seq(xlim_l,xlim_u,by=0.1)
  ylim_u <- 0.4
  ylim_l <- 0
  ylim_v <- seq(ylim_l,ylim_u,by=0.1)
  
  v=result$v
  lv=length(result$v)
  F0=result$F.0
  F50=result$F.50
  F90=result$F.90
  
  ymax <- ceiling(max(c(F0, F50, F90)*20*1.3)/4)*4/20
  ylim_v <- seq(0,ymax,by=ymax/4)
  
  
  histogram_color <- adjustcolor("#bddcbd", alpha.f = 0.9)
  
  temp <- na.omit(df_sub[,var])
  mean <- mean(temp)
  sd <- sd(temp)
  temp <- pnorm(temp, mean=mean, sd=sd)
  h = hist(temp, plot=F, breaks = seq(0,1,0.1))
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,xlab="", ylab="", ylim=c(0,100), 
       axes=FALSE, col=histogram_color, main="", border="#bddcbd")
  
  
  ax <- seq(0,100,25)
  axis(4,las=1,at=ax,labels=paste0(ax,"%"), cex.axis=1)
  mtext(side=4, line=2.5,  "Percentage", cex=0.9)
  
  
  par(new=TRUE)
  plot(v,F0, xlab=" ",ylab=" " , xlim=c(xlim_l,xlim_u), ylim=c(0, ymax), 
       lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5, cex.sub=1.5,cex.main=1.5,axes=F, type='l')
  
  lines(v,F0, lty=1,lwd=2, col = "#2babd5" )
  lines(v,F50, lty=2,lwd=2, col = "#2babd5" )
  lines(v,F90, lty=3,lwd=2, col = "#2babd5" )
  
  abline(h=1,lty=1,lwd=1, col = "gray")
  
  #axis(1, xlim_v, cex.axis=1)
  if (iv>0){
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv],1), 
         cex.axis=0.9)
  } else{
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv]), 
         cex.axis=0.9)
  }
  
  # mtext(side=1, line=2.5, 'AA sequence distance v', cex=1.5)
  mtext(side=1, line=2.2, lab_list[iv], cex=0.9)
  
  axis(2,las=1,at=ylim_v,labels=ylim_v,cex.axis=1)
  mtext(side=2, line=2.5, "CIF rate",cex=1.0)
  
  legend(0.032, ymax, c(paste("D15 IM 10th percentile:", round(10^result$qtl[1], 0)),
                        paste("D15 IM 50th percentile:", round(10^result$qtl[2], 0)),
                        paste("D15 IM 90th percentile:", round(10^result$qtl[3], 0))),
         col=c("#2babd5","#2babd5","#2babd5") ,lwd=c(2,2,2),lty=c(1,2,3),cex=0.8, bty="n")
  box()  
  
  mtext(panel_titles[iv], side = 3, line = 1.5, adj = 0, font = 2, cex = 1)
}

mtext("Figure 5: Distance-specific COVID-19 CIF rates for Day 15 nAb ID50 titer",
      outer = TRUE, cex = 1.5, font = 2)

dev.off()

#### Figure S3 ####

pdf("FigureS3_CIF_panels_naive.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), oma = c(2, 2, 4, 0), mar = c(5.1, 5.1, 4, 5))

panel_titles <- c("A - Prototype", "B - Omicron", "C - Omicron", "D - Omicron")


V_list <- list(
  "zhd.spike",                # A
  "zhd.spike",                # B
  "zhd.rbd",                  # C
  "hybrid.vx.rbd.rescaled"   # D
)

lab_list <- c("Spike Hamming Distance",
              "Spike Hamming Distance",
              "RBD Hamming Distance",
              "DMS-escape RBD-1")

mean_list <- c(30.8502, 24.6977, 8.2338, 0.3110*max(df_marks_all$hybrid.vx.rbd, na.rm = T))
sd_list <- c(7.2396, 3.7868, 2.5885, 0.1884*max(df_marks_all$hybrid.vx.rbd, na.rm = T))

mat_files <- list(
  "res/a2_3_V1_L5M10.mat",   # A
  "res/a1_3_V1_L5M10.mat",   # B
  "res/a1_3_V2_L5M10.mat",   # C
  "res/a1_3_V3_L5M10.mat"    # D
)

df_list <- list(
  df %>% filter(TrtB == 1, naive == 1),  # A
  df %>% filter(TrtB == 0, naive == 1),  # B
  df %>% filter(TrtB == 0, naive == 1),  # C
  df %>% filter(TrtB == 0, naive == 1)   # D
)

for (iv in 1:4) {
  var <- V_list[[iv]]
  df_sub <- df_list[[iv]]
  matfilename <- mat_files[[iv]]
  
  result <- readMat(matfilename)
  
  xlim_u <- 1
  xlim_l <- 0
  xlim_v <- seq(xlim_l,xlim_u,by=0.1)
  ylim_u <- 0.4
  ylim_l <- 0
  ylim_v <- seq(ylim_l,ylim_u,by=0.1)
  
  v=result$v
  lv=length(result$v)
  F0=result$F.0
  F50=result$F.50
  F90=result$F.90
  
  ymax <- ceiling(max(c(F0, F50, F90)*20*1.3)/4)*4/20
  ylim_v <- seq(0,ymax,by=ymax/4)
  
  
  histogram_color <- adjustcolor("#bddcbd", alpha.f = 0.9)
  
  temp <- na.omit(df_sub[,var])
  mean <- mean(temp)
  sd <- sd(temp)
  temp <- pnorm(temp, mean=mean, sd=sd)
  h = hist(temp, plot=F, breaks = seq(0,1,0.1))
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,xlab="", ylab="", ylim=c(0,100), 
       axes=FALSE, col=histogram_color, main="", border="#bddcbd")
  
  
  ax <- seq(0,100,25)
  axis(4,las=1,at=ax,labels=paste0(ax,"%"), cex.axis=1)
  mtext(side=4, line=2.5,  "Percentage", cex=0.9)
  
  
  par(new=TRUE)
  plot(v,F0, xlab=" ",ylab=" " , xlim=c(xlim_l,xlim_u), ylim=c(0, ymax), 
       lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5, cex.sub=1.5,cex.main=1.5,axes=F, type='l')
  
  lines(v,F0, lty=1,lwd=2, col = "#2babd5" )
  lines(v,F50, lty=2,lwd=2, col = "#2babd5" )
  lines(v,F90, lty=3,lwd=2, col = "#2babd5" )
  
  abline(h=1,lty=1,lwd=1, col = "gray")
  
  #axis(1, xlim_v, cex.axis=1)
  if (iv>0){
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv],1), 
         cex.axis=0.9)
  } else{
    axis(1, 
         at=seq(0.1,0.9,0.1), 
         labels=round(qnorm(seq(0.1,0.9,0.1))*sd_list[iv]+mean_list[iv]), 
         cex.axis=0.9)
  }
  
  # mtext(side=1, line=2.5, 'AA sequence distance v', cex=1.5)
  mtext(side=1, line=2.2, lab_list[iv], cex=0.9)
  
  axis(2,las=1,at=ylim_v,labels=ylim_v,cex.axis=1)
  mtext(side=2, line=2.5, "CIF rate",cex=1.0)
  
  legend(0.032, ymax, c(paste("D15 IM 10th percentile:", round(10^result$qtl[1], 0)),
                        paste("D15 IM 50th percentile:", round(10^result$qtl[2], 0)),
                        paste("D15 IM 90th percentile:", round(10^result$qtl[3], 0))),
         col=c("#2babd5","#2babd5","#2babd5") ,lwd=c(2,2,2),lty=c(1,2,3),cex=0.8, bty="n")
  box()
  
  mtext(panel_titles[iv], side = 3, line = 1.5, adj = 0, font = 2, cex = 1)
}

mtext("Supplementary Figure 3: Distance-specific COVID-19 CIF rates for Day 15 nAb ID50 titer",
      outer = TRUE, cex = 1.5, font = 2)

dev.off()

#### Figure 6: Pairwise Correlations of Day 15 nAb ID50 Titers ####

# Prepare data
df_plot <- df %>%
  filter(TrtB == 0) %>%
  filter(!is.na(Day15pseudoneutid50_D614G),
         !is.na(Day15pseudoneutid50_BA.1),
         !is.na(Day15pseudoneutid50_BA.4.BA.5)) %>%
  mutate(
    group = ifelse(naive == 1, "Naïve", "Non-naïve"),
    group = factor(group, levels = c("Naïve", "Non-naïve"))
  )

# Color and shape
color_palette <- c("Naïve" = "#E69F00", "Non-naïve" = "#D55E00")
shape_palette <- c("Naïve" = 15, "Non-naïve" = 16)

# Spearman rho + 95% CI annotation
cor_label <- function(x, y) {
  rho <- cor(x, y, method = "spearman", use = "complete.obs")
  boot_rho <- replicate(1000, {
    idx <- sample(seq_along(x), replace = TRUE)
    cor(x[idx], y[idx], method = "spearman", use = "complete.obs")
  })
  ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)
  sprintf("rho = %.2f\n[%.2f, %.2f]", rho, ci[1], ci[2])
}

# Custom scatter plot with rho text
create_plot <- function(df, xvar, yvar, xlabel, ylabel, xlim_range, ylim_range) {
  label <- cor_label(df[[xvar]], df[[yvar]])
  ggplot(df, aes_string(x = xvar, y = yvar, color = "group", shape = "group")) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", size = 0.4) +
    annotate("text", x = xlim_range[1], y = ylim_range[2], hjust = 0, vjust = 1,
             label = label, size = 3.5, fontface = "bold") +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    labs(x = xlabel, y = ylabel) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.margin = margin(5, 5, 5, 5),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# Set axis limits
xlim_range <- c(1, 5.5)
ylim_range <- c(1, 5.5)

# 3 scatter plots
p1 <- create_plot(df_plot, "Day15pseudoneutid50_D614G", "Day15pseudoneutid50_BA.1", "D614G", "BA.1", xlim_range, ylim_range)
p2 <- create_plot(df_plot, "Day15pseudoneutid50_D614G", "Day15pseudoneutid50_BA.4.BA.5", "D614G", "BA.4/5", xlim_range, ylim_range)
p3 <- create_plot(df_plot, "Day15pseudoneutid50_BA.1", "Day15pseudoneutid50_BA.4.BA.5", "BA.1", "BA.4/5", xlim_range, ylim_range)

# Extract legend
legend_plot <- ggplot(df_plot, aes(x = 1, y = 1, color = group, shape = group)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme_void() +
  theme(legend.position = "right")
legend_grob <- get_legend(legend_plot)
legend_gg <- as_ggplot(legend_grob)

# Arrange panels
panel_grid <- grid.arrange(p1, legend_gg, p2, p3,
                           ncol = 2, nrow = 2,
                           widths = c(1, 1), heights = c(1, 1))

# Add title
final_plot <- annotate_figure(panel_grid,
                              top = text_grob("Figure 6: Pairwise Correlations of Day 15 nAb ID50 Titers",
                                              face = "bold", size = 14)
)

# Save
ggsave("Figure6_pairs_IM_rho_2x2.pdf", final_plot, width = 10, height = 8)
