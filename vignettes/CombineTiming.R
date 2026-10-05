require(SLIMpaper)
library(dplyr)
library(forcats)
library(ggplot2)
library(gridExtra)
library(scales)

name.directory <- "Output_timing"
files <- list.files(name.directory, full.names = TRUE,recursive = TRUE)

outputs <- do.call("rbind", lapply(files, function(f) {
  o <- readRDS(f)
  return(o)
}))

timings <- outputs %>%
  group_by(n, p, n.samp, method) %>%
  summarise(mean_time = mean(time), sd_time = sd(time),
            lwr = quantile(time, 0.025),
            upr = quantile(time, 0.975)) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = c("binary program", "L.R. B.P.", "relaxed B.P."), labels = c("B.P.", "Lagrangian\nB.P.", "Relaxed\nB.P.")))

pal <- ggsci::pal_jama()(7)[c(1,6,2)]

#### plots by N ####
plot_n <- timings %>%
  filter(p == 21) %>%
  filter(n.samp == 100) %>%
  ggplot(aes(x = n, y = mean_time, fill = method, color = method, group = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = 0.2) +
  labs( x = "Sample Size (n)", y = "Time (s)") +
  theme_bw(11) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.6)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 2^16)) +
  scale_y_continuous(trans = "sqrt",
                     limits = c(0,25),
                     breaks = c(0, 1, 4, 10, 25),
                     expand = c(0,0)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal)

#### plots by P ####
plot_k_width <- timings %>%
  filter(n == 1024) %>%
  filter(n.samp == 100) %>%
  ggplot(aes(x = p, y = mean_time, fill = method, color = method, group = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = 0.2) +
  labs(title = "a.", x = "", y = "Time (s)") +
  coord_cartesian(ylim = c(0,11)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,504)) +
  theme_bw(11) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,1, 5, 10),
                     expand = c(0,0)
                     ) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  ggplot2::theme(legend.position="none")

plot_k_both <- timings %>%
  filter(n == 1024) %>%
  filter(n.samp == 100) %>%
  ggplot(aes(x = p, y = mean_time, fill = method, color = method, group = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = 0.2) +
  labs(title = "b.", x = "Number of Parameters (k)", y = "") +
  coord_cartesian(xlim = c(20,71), ylim = c(0.01,60)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(trans = "log", breaks = c(0.01, 0.3, 2, 50, 1000, 22000)) +
  theme_bw(11) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  ggplot2::theme(legend.position="bottom")

k_legend <- get_legend(plot_k_both)

plot_k_both <- plot_k_both + theme(legend.position="none")

plot_k_all <- timings %>%
  filter(n == 1024) %>%
  filter(n.samp == 100) %>%
  ggplot(aes(x = p, y = mean_time, fill = method, color = method, group = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = 0.2) +
  labs(title = "c.", x = "", y = "") +
  coord_cartesian(xlim = c(20,51)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(trans = "log", breaks = c(0.01,0.1, 2,50, 1000, 22000), labels = scales::label_comma()) +
  theme_bw(11) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  ggplot2::theme(legend.position="none")

plot_k_p <- arrangeGrob(plot_k_width, plot_k_both, plot_k_all, ncol = 3)
plot_k <- arrangeGrob(plot_k_p, k_legend, ncol = 1, heights = c(1,0.1))

#### plots by T (samples) ####
plot_T <- timings %>%
  filter(n == 1024) %>%
  filter(p == 21) %>%
  ggplot(aes(x = n.samp, y = mean_time, fill = method, color = method, group = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr,color = NULL), alpha = 0.2) +
  labs(x = "Number of Samples (T)", y = "Time (s)") +
  theme_bw(11) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4000)) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,4, 25,100, 250),
                     limits = c(0,300),
                     expand = c(0,0)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal)

#### Save plots ####
pdf(file.path("inst", "figure", "timing", "timing_n.pdf"),
    width = 7, height = 3.5)
print(plot_n)
dev.off()

pdf(file.path("inst", "figure", "timing", "timing_k.pdf"),
    width = 7, height = 3.5)
grid.arrange(plot_k)
dev.off()

pdf(file.path("inst", "figure", "timing", "timing_T.pdf"),
    width = 7, height = 3.5)
print(plot_T)
dev.off()

