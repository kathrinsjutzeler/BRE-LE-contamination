# Fst statistics for Contamination project
# Author: Kathrin Jutzeler
# Date: October , 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(ggridges)
library(tidyverse)
library(locfit)
library(ggpubr)
library(patchwork)

#--------------#
# Functions ####
#--------------#

# Calculate genome wide Fst 
f_fst_gw <- function(df){ 
  df %>%
    filter(!CHROM %in% c("SM_V10_WSR", "SM_V10_Z", "SM_V10_MITO")) %>%
    separate(FST, c("pops", 'fst'), sep = "=") %>%
    mutate(fst = as.numeric(fst)) %>%
    summarize(mean = mean(fst))
}

# Calculate Fst across genome in 20000 base windows
f_fst <- function(df){
  df %>%
    separate(FST, c("pops", 'fst'), sep = "=") %>%
    mutate(fst = as.numeric(fst)) %>%
    group_by(CHROM, sex, pos = cut(POSITION, breaks = seq(0, max(POSITION), 20000))) %>%
    summarize(mean = mean(fst))
}

# Calculate Fst across genome per base
f_fst_no_window <- function(df){
  df %>%
    separate(FST, c("pops", 'fst'), sep = "=") %>%
    mutate(fst = as.numeric(fst)) 
}

#---------------------#
# Import the data ####
#---------------------#
setwd('fst_BRE_LE')

files <- list.files()
list <- as.list(files[7:24])
names(list) <- str_sub(files[7:24], end =-5)

ldf <- lapply(list, read.delim,header =F)

ldf <- lapply(ldf, setNames, c('CHROM', 'POSITION', 'v3', 'V4', 'COVERAGE', 'FST'))


ldf[c(1,3,5,7,10,13,15, 17)] <- lapply(ldf[c(1,3,5,7,10,13,15,17)], function(x) x %>% mutate(sex = "f"))
ldf[c(2,4,6,8,9,11,12, 14,16,18)] <- lapply(ldf[c(2,4,6,8,9,11,12,14,16,18)], function(x) x %>% mutate(sex = "m"))

genome_wide_fst <- lapply(ldf, f_fst_gw)

gw_fst_df <- bind_rows(genome_wide_fst, .id = 'comp')
write.csv(gw_fst_df, 'gw_fst_df_auto.csv')

BRE_LE_list <- lapply(ldf, f_fst) 

BRE_LE_df <- bind_rows(BRE_LE_list, .id ='sample')

BRE_LE_df <- BRE_LE_df %>%
  mutate(time = str_sub(sample, start = 3, end =8))

BRE_LE_df <- BRE_LE_df %>%
  filter(CHROM != 'SM_V10_MITO') %>%
  mutate(smooth = runmed(mean, 11))

# Put time in actual order
BRE_LE_levels <- c("110216","062118", "051420", "112320", "070521", "092921", "102621",  "122121", 
                   "070522","021523")
BRE_LE_df$time <- factor(BRE_LE_df$time, levels= BRE_LE_levels)

# Rename time
#LE_time.labs <- c("t0_v_t1", "t0_v_t2", "t0_v_t3", "t0_v_t4", "t0_v_t5", "t0_v_t6", "t0_v_t7")
#names(LE_time.labs) <- LE_levels

# with mean 

my_color <- c("m" = "black", "f" = "red")

chrom.labs <- c('1', '2', '3', '4', '5', '6', '7', 'W', 'Z')
names(chrom.labs) <- levels(factor(BRE_LE_df$CHROM))

# Plot Fig 2 - FST BRE v LE####
p_BRE_LE <- ggplot(subset(BRE_LE_df, CHROM != "SM_V10_MITO"), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.5,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x", space = 'free_x', labeller = labeller(CHROM= chrom.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Chromosome", y = bquote(bold("Mean  " * F[ST] * " (LE vs BRE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16)) +
  ggh4x::force_panelsizes(cols = c(8,4.5,4.5,4.5,3,3,2,0.8,8))

jpeg('Figure2_smooth.jpeg', width =10, height =12, unit ='in', res = 300, bg = 'white')
p_BRE_LE
dev.off()

p_BRE_LE_smooth <- ggplot(subset(BRE_LE_df, CHROM != "SM_V10_MITO"), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.5,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x") + #, labeller = labeller(time = time.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Genomic Position", y = "Mean Fst (LE vs BRE)") +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

ggplot(subset(BRE_LE_df, CHROM== "SM_V10_1"), aes(pos, mean, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.5,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x") + #, labeller = labeller(time = time.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Genomic Position", y = "Mean Fst (LE vs BRE)") +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())


chr1 <- BRE_LE_df %>%
  filter(CHROM == 'SM_V10_1', time == '021523') 
summary(chr1$mean)
filter(chr1, mean > 0.3)


chr7 <- BRE_LE_df %>%
  filter(CHROM == 'SM_V10_7', time == '021523') 
summary(chr7$mean)
filter(chr7, mean > 0.5)

chrz <- BRE_LE_df %>%
  filter(CHROM == 'SM_V10_Z', time == '021523') 
summary(chrz$mean)
filter(chrz, mean > 0.4)


#pdf("Fst_BRE_LE_by_time.pdf", width = 10, height = 12)

#(p_BRE_LE) &  theme(plot.tag = element_text(face = 'bold'))
#+ plot_annotation(tag_levels = c("BRE")) 

#dev.off()

png("Fst_BRE_LE_by_time.png", width = 10, height = 12, units = 'in', res = 600)

(p_BRE_LE) &  theme(plot.tag = element_text(face = 'bold'))
#+ plot_annotation(tag_levels = c("BRE")) 

dev.off()

png("Fst_BRE_LE_by_time_smooth.png", width = 10, height = 12, units = 'in', res = 150)

(p_BRE_LE_smooth) &  theme(plot.tag = element_text(face = 'bold'))
#+ plot_annotation(tag_levels = c("BRE")) 

dev.off()


# Plot Fst by time ####
files <- list.files('fst_time/')
list <- as.list(paste0('fst_time/', files[1:34]))
names(list) <- str_sub(files[1:34], end =-5)

ldf2 <- lapply(list, read.delim,header =F)

ldf2 <- lapply(ldf2, setNames, c('CHROM', 'POSITION', 'v3', 'V4', 'COVERAGE', 'FST'))

ldf2[c(1:8, 18:25)] <- lapply(ldf2[c(1:8, 18:25)], function(x) x %>% mutate(sex = "f"))
ldf2[c(9:17, 26:34)] <- lapply(ldf2[c(9:17, 26:34)], function(x) x %>% mutate(sex = "m"))


# Genome wide across time 
BRE_fst_wg <- lapply(ldf2[1:17], f_fst_gw) 
LE_fst_wg <- lapply(ldf2[18:34], f_fst_gw) 

# Simplified plot
BRE_df_wg <- bind_rows(BRE_fst_wg, .id ='comp')
BRE_df_wg <- BRE_df_wg %>%
  mutate(time = str_sub(comp, start = 16, end =-3), sex = c(rep('f', 8), rep('m', 9)), pop = 'BRE')
BRE_df_wg$time <- factor(BRE_df_wg$time, levels= BRE_levels)

LE_df_wg <- bind_rows(LE_fst_wg, .id ='comp')
LE_df_wg <- LE_df_wg %>%
  mutate(time = str_sub(comp, start = 14, end =-3), sex = c(rep('f', 8), rep('m', 9)), pop = 'LE')
LE_df_wg$time <- factor(LE_df_wg$time, levels= LE_levels)

levels(LE_df_wg$time) <- c("June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                           "Dec 2021","July 2022", "Feb 2023")

levels(BRE_df_wg$time) <-c("June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                           "Dec 2021","July 2022", "Feb 2023")

time_df <- bind_rows(BRE_df_wg, LE_df_wg)

# Figure 3 Plot ####

ggplot(time_df, aes( time, mean,  group = interaction(sex, pop), linetype = sex, color = pop)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  #scale_color_manual(values = c("#b0d4a4", "#F98400" )) +
  scale_color_manual(values = c("black", "red" )) +
  labs(y = bquote(bold("Mean " * F[ST])), x = 'Time comparison', linetype = 'Sex', color = 'Population') + 
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
        axis.line.x = element_line(), axis.line.y = element_line(),
        text = element_text(size =16), axis.text.x = element_text(angle =45, hjust = 1),
        panel.grid = element_blank(), axis.title = element_text(face = 'bold'), axis.ticks = element_line())

ggsave('Figure3v2.jpg', width = 8, height = 6, dpi = 300, units = 'in', bg = 'white')  


# Detailed plot by time ####
BRE_list <- lapply(ldf2[1:17], f_fst) 
LE_list <- lapply(ldf2[18:34], f_fst) 

BRE_df <- bind_rows(BRE_list, .id ='sample')
LE_df <- bind_rows(LE_list, .id ='sample')

BRE_df <- BRE_df %>%
  mutate(time = str_sub(sample, start = 16, end =-3))

LE_df <-  LE_df %>%
  mutate(time = str_sub(sample, start = 14, end =-3))

BRE_df <- BRE_df %>%
  mutate(smooth = runmed(mean, 21))

LE_df <- LE_df %>%
  mutate(smooth = runmed(mean, 21))

# Put time in actual order
BRE_levels <- c("062718", "051420", "112320", "070521", "092921", "102621", "122121", "070522","021523")
BRE_df$time <- factor(BRE_df$time, levels= BRE_levels)

LE_levels <- c("062118", "051420", "112320", "070521", "092921", "102621",  "122121", "070522","021523")
LE_df$time <- factor(LE_df$time, levels= LE_levels)


# Rename time
comp_time_order <- c("June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                     "Dec 2021","July 2022", "Feb 2023")

BRE_time.labs <- comp_time_order
names(BRE_time.labs) <- BRE_levels

LE_time.labs <- comp_time_order
names(LE_time.labs) <- LE_levels

#my_color <- c("m" = "#3030D0", "f" = "#D12959")

## Supplement Fst by time ####
p_BRE_by_time <- ggplot(subset(BRE_df, CHROM != "SM_V10_MITO"), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.25,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x", space = 'free_x',
             labeller = labeller(time= BRE_time.labs, CHROM = chrom.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              linewidth =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Chromosome", y = bquote(bold("Mean " * F[ST] * " (BRE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
             axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16)) +
   ggh4x::force_panelsizes(cols = c(8,4.5,4.5,4.5,3,3,2,0.8,8))
       
ggsave('FigureS1.jpg', width = 10, height =12, units = 'in', dpi = 300, bg = 'white')

#pdf("Fst_BRE_by_time.pdf", width = 10, height = 12)

#(p_BRE_by_time)
#+ plot_annotation(tag_levels = c("BRE")) &  theme(plot.tag = element_text(face = 'bold'))

#dev.off()


p_LE_by_time <- ggplot(subset(LE_df, CHROM != "SM_V10_MITO"), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.25,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x", space = 'free_x',
             labeller = labeller(time= LE_time.labs, CHROM = chrom.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Chromosome", y = bquote(bold("Mean " * F[ST] * " (LE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16)) +
  ggh4x::force_panelsizes(cols = c(8,4.5,4.5,4.5,3,3,2,0.8,8))

ggsave('FigureS2.jpg', width = 10, height =12, units = 'in', dpi = 300, bg = 'white')


jpeg("Fst_by_time.jpg", width = 20, height = 12, units = 'in', res = 600)

((p_LE_by_time | p_BRE_by_time) + plot_layout(guides = "collect")) + plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'), legend.position = 'bottom')

dev.off()

#pdf("Fst_LE_by_time.pdf", width = 10, height = 12)

#(p_LE_by_time) & #plot_annotation(tag_levels = c("LE")) 
 #   theme(plot.tag = element_text(face = 'bold'))

#dev.off()


# LE specific only ####
LE_list_SNP <- lapply(ldf2[18:34], f_fst_no_window) 

LE_df_SNP <- bind_rows(LE_list_SNP, .id ='sample')  

LE_df_SNP <-  LE_df_SNP %>%
  mutate(time = str_sub(sample, start = 14, end =-3))

LE_df_SNP %>%
  dplyr::select(CHROM, POSITION) %>%
  unique()

##### Get a list of SNPs with high FST ####
SNP_df <- LE_df_SNP %>%
  arrange(CHROM, POSITION) %>%
  group_by(CHROM, POSITION) %>%
  count() %>%
  filter(n == 17)

dim(SNP_df)

high_fst_SNPs <- LE_df_SNP %>%
  right_join(SNP_df, by = c("CHROM", "POSITION")) %>%
  arrange(desc(fst))

head(high_fst_SNPs)


#select_SNPs %>%
#  filter(interval == "SM_V10_2:43301696")

t <- high_fst_SNPs %>%
  filter(fst >= 0.2)

dim(t)

select_SNPs <- high_fst_SNPs[1:749,] %>%
  arrange(CHROM, POSITION)
  
select_SNPs %>%
  group_by(CHROM, POSITION) %>%
  count() %>%
  arrange(desc(n))


select_SNPs$interval <- paste0(select_SNPs$CHROM, ":", 
                                    select_SNPs$POSITION, '-', select_SNPs$POSITION)


SNP_list <- unique(select_SNPs$interval)

LE_df_SNP %>%
  filter(CHROM == 'SM_V10_2', POSITION == 43301696)

ggplot(select_SNPs, aes(POSITION, fst, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 1,alpha = .1) +
  facet_grid(time~CHROM, switch = "x", scales = "free_x") + #, labeller = labeller(time = time.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Genomic Position", y = "Mean Fst (LE)", color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())


write_delim(data.frame(SNP_list), 'select_LE_SNPs.list')

##### Histograms to show metric for selection of SNPs ####
p_full_histo_fst <- ggplot(LE_df_SNP, aes(fst, fill = sex)) +
  geom_histogram(position = position_dodge(0.01)) + 
  #facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  labs(y = 'Count', x = bquote(bold(F[ST])), fill = "Sex") +
  theme_minimal() +
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title = element_text(face = 'bold'), panel.grid = element_blank())

p_select_histo_fst <- 
ggplot(select_SNPs, aes(fst, fill = sex)) +
  geom_histogram(position = position_dodge(0.01)) + 
  #facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  labs(y = 'Count', x = bquote(bold(F[ST])), fill = "Sex") +
  theme_minimal()+
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title = element_text(face = 'bold'), panel.grid = element_blank())


jpeg('Figure5b.jpg', width = 10, height =3, units = 'in', res =300)
(p_full_histo_fst | p_select_histo_fst)  + plot_layout(guides = "collect") 
  #plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))
dev.off()

pdf('Figure_histograms.pdf', width = 8, height = 10)

(p_full_histo_fst + theme(legend.position = 'none')) / p_select_histo_fst + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A")) &  
  theme(plot.tag = element_text(face = 'bold', size = 18), 
        legend.position = 'bottom', panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'))


dev.off()

#Not used
# Lollipot plot for distribution of SNPs 
unique_SNPs <- select_SNPs %>%
  group_by(CHROM, POSITION) %>%
  reframe(mean_fst = mean(fst))

unique_SNPs$POSITION / 0.5
  
ggplot(unique_SNPs, aes(POSITION, mean_fst)) +
  geom_segment( aes(x=POSITION, xend=POSITION, y=0, yend=mean_fst)) +
  geom_point(color = LE_col, alpha = .8) +
  facet_grid(~CHROM, switch = "x", scales = 'free') +
  theme(axis.text.x = element_text(angle =90)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Genomic Position", y = "Mean Fst", color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(),
        axis.title = element_text(size = 14, face = 'bold'))

ggsave('SNP_dist.png', width = 8, height = 6)

  

