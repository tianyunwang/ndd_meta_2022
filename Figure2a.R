
library(dplyr)
library(ggplot2)
library(qqman)
library(ggrepel)
library(tidyr)
library(scatterpie)

df_input = read.csv("46621NDD_min_p_q.csv")
head(df_input)

df_input$MinQraw_log_adj = -log10(df_input$MinQ) 
df_input$MinQraw_log = -log10(df_input$MinQ_raw)
#rawlog_rank = order(df_input$MinQraw_log, decreasing = TRUE)
#indexes = which(df_input$MinQraw_log[rawlog_rank] >50)
#df_input$MinQraw_log_adj = df_input$MinQraw_log
#df_input$MinQraw_log_adj[indexes] = 50 + 2*seq(length(indexes), 1, -1)

################

min_max_norm <- function(x, values) {
  result <- (values-min(x)) / (max(x)-min(x))
  return(result)
}
s_factor = 100

## compute the cumulative position of SNP
df = df_input %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df_input, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP and normalize axis
  # because BPcum (x-axis) and log10(MinQ)(y-axis) are on too different scales 
  # they need to be normalized so the piecharts are draw properly
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot,
         xnorm = s_factor*min_max_norm(BPcum, BPcum),
         ynorm = s_factor*min_max_norm(MinQraw_log_adj, MinQraw_log_adj))

## prepare axis
x_axis = df %>% group_by(CHR) %>% summarize(center = (max(xnorm) + min(xnorm)) / 2 )

# prepare thresholds for hlines:
fdr5 = -log10(5e-2)
fwer5 = -log10(5e-5)
fdr5_norm = s_factor*min_max_norm(df$MinQraw_log_adj, -log10(5e-2)) # FDR 5% threshold
fwer5_norm = s_factor*min_max_norm(df$MinQraw_log_adj, -log10(5e-5)) # FWER 5% threshold



ylabels <- seq(0,50,10)
yticks <- min_max_norm(df_input$MinQraw_log_adj/s_factor, ylabels)
ylabels <- c(ylabels, round(max(df_input$MinQraw_log)))
yticks <- c(yticks, s_factor)

ggplot(df, aes(x = xnorm, y = ynorm, label = SNP)) +
  
  geom_point(aes(size = DNM_All/8), alpha = 0) +
  # Show all points
  geom_point(data = subset(df, MinQraw_log <= fdr5), aes(color = as.factor(CHR), size = DNM_All/8), alpha = 0.8, shape=16) +
  scale_color_manual(values = c(rep(c("grey77", "grey29"), 11), c("grey"))) + 
  #scale_color_manual(values = c(rep(c("deepskyblue", "darkblue"), 11 ), "deepskyblue")) +
  
  # add pie chars
  geom_scatterpie(data = subset(df, MinQraw_log >= fdr5), aes(x = xnorm, y = ynorm, r = DNM_All/32),
                  cols = c("ASD_LGD", "ASD_MIS", "DD_LGD", "DD_MIS"),
                  color = NA) +
  scale_fill_manual(values = c("mediumpurple4", "mediumpurple", "darkgreen", "mediumseagreen")) +
  geom_scatterpie_legend(radius = df$DNM_All/32, x = 96, y = 65, n=5, labeller = function(x) round(x*32*1.1)) +
  
  # gene labels
  geom_text(data = subset(df, FWER138 >0), size = 3, fontface = "italic", 
            nudge_x = 3, nudge_y = 0, check_overlap = TRUE) +
  
  # custom axis:
  scale_x_continuous(label = x_axis$CHR, breaks = x_axis$center) +
  scale_y_continuous(label = ylabels, breaks = yticks, expand = c(0, 0), name = "-log10(MinQ)") +
  
  # threshold for hline:
  geom_hline(yintercept = fdr5_norm, color = "black", linetype = "dashed", size = 0.8) + # FDR 5% threshold
  geom_hline(yintercept = fwer5_norm, color = "red", linetype = "dashed", size = 0.8) + # FWER 5% threshold
  
  # Custom the theme:
  theme_bw() +
  theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  guides(color = FALSE, size = FALSE)

ggsave("manhattan_piechart_overFDR.pdf", width = 10*1.309565, height = 10, device='pdf', dpi=700)

