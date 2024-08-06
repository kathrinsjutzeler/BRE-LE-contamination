# Fst statistics for Contamination project
# Author: Kathrin Jutzeler
# Date: October , 2023
# Last updated: August 6, 2024
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

#----------------------#
# FST by population  ####
#----------------------#
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

# with mean 

my_color <- c("m" = "black", "f" = "red")

# Labels
chrom.labs <- c('1', '2', '3', '4', '5', '6', '7', 'W', 'Z')
names(chrom.labs) <- levels(factor(BRE_LE_df$CHROM))

comp_time_order <- c("Nov 2016", "June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                     "Dec 2021","July 2022", "Feb 2023")

time.labs <- comp_time_order
names(time.labs) <- BRE_LE_levels

## Plot Fig 2 - FST BRE v LE####
p_BRE_LE <- ggplot(subset(BRE_LE_df, !CHROM %in% c("SM_V10_WSR", "SM_V10_MITO")), 
                   aes(pos, smooth, color=sex, group=sex)) +
  geom_point(size = 0.5,alpha = .1) +
  ggh4x::facet_grid2(time~CHROM, switch = "x", scales = "free_x", space = 'free_x', 
             labeller = labeller(CHROM= chrom.labs, time = time.labs), axes = 'x', remove_labels = "all") +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x= "Chromosome", y = bquote(bold("Mean  " * F[ST] * " (LE vs BRE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16),
        axis.line = element_line()) 

jpeg('Figure2_smooth.jpeg', width =10, height =12, unit ='in', res = 300, bg = 'white')
p_BRE_LE
dev.off()


#------------------#
# FST by time   ####
#------------------#

files <- list.files('fst_time/')
list <- as.list(paste0('fst_time/', files[1:34]))
names(list) <- str_sub(files[1:34], end =-5)

ldf2 <- lapply(list, read.delim,header =F)

ldf2 <- lapply(ldf2, setNames, c('CHROM', 'POSITION', 'v3', 'V4', 'COVERAGE', 'FST'))

ldf2[c(1:8, 18:25)] <- lapply(ldf2[c(1:8, 18:25)], function(x) x %>% mutate(sex = "f"))
ldf2[c(9:17, 26:34)] <- lapply(ldf2[c(9:17, 26:34)], function(x) x %>% mutate(sex = "m"))

ldf2 <- lapply(ldf2, function(x) x %>% filter(CHROM != 'SM_V10_WSR'))

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

## Figure 3 Plot ####

jpeg('Figure3.jpg', width = 6, height =4, units = 'in', res =300)

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
dev.off()


## Plot Fig S1 and S2 Fst by time ####
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

## Plot BRE####
p_BRE_by_time <- ggplot(subset(BRE_df, !CHROM %in% c("SM_V10_WSR", "SM_V10_MITO")), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.25,alpha = .1) +
  ggh4x::facet_grid2(time~CHROM, switch = "x", scales = "free_x", space = 'free_x', 
                     labeller = labeller(CHROM= chrom.labs, time = BRE_time.labs), axes = 'x', remove_labels = "all") +
  #facet_grid(time~CHROM, switch = "x", scales = "free_x", space = 'free_x',
   #          labeller = labeller(time= BRE_time.labs, CHROM = chrom.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              linewidth =1, se = F, aes(color = sex)) +
  theme_minimal() +
  #scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = my_color) +
  labs(x= "Chromosome", y = bquote(bold("Mean " * F[ST] * " (BRE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
             axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16),
        axis.line = element_line())
       
ggsave('FigureS1.jpg', width = 10, height =12, units = 'in', dpi = 300, bg = 'white')

## Plot LE####
p_LE_by_time <- ggplot(subset(LE_df, !CHROM %in% c("SM_V10_WSR", "SM_V10_MITO")), aes(pos, smooth, color=sex, group=sex)) +
  #geom_line() + 
  geom_point(size = 0.25,alpha = .1) +
  ggh4x::facet_grid2(time~CHROM, switch = "x", scales = "free_x", space = 'free_x', 
                     labeller = labeller(CHROM= chrom.labs, time = LE_time.labs), axes = 'x', remove_labels = "all") +
  #facet_grid(time~CHROM, switch = "x", scales = "free_x", space = 'free_x',
   #          labeller = labeller(time= LE_time.labs, CHROM = chrom.labs)) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
              size =1, se = F, aes(color = sex)) +
  theme_minimal() +
  scale_color_manual(values = my_color) +
  labs(x= "Chromosome", y = bquote(bold("Mean " * F[ST] * " (LE)")), color = 'Sex') +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 16),
        axis.line = element_line())

ggsave('FigureS2.jpg', width = 10, height =12, units = 'in', dpi = 300, bg = 'white')

#---------------------#
# LE specific only ####
#---------------------#

files <- list.files('fst_time_LE_all/')
list <- as.list(paste0('fst_time_LE_all/', files[1:17]))
names(list) <- str_sub(files[1:17], end =-5)

ldf3 <- lapply(list, read.delim,header =F)

ldf3 <- lapply(ldf3, setNames, c('CHROM', 'POSITION', 'v3', 'V4', 'COVERAGE', 'FST'))

ldf3[c(1:8)] <- lapply(ldf3[c(1:8)], function(x) x %>% mutate(sex = "f"))
ldf3[c(9:17)] <- lapply(ldf3[c(9:17)], function(x) x %>% mutate(sex = "m"))

ldf3 <- lapply(ldf3, function(x) x %>% filter(CHROM != 'SM_V10_WSR'))

LE_list_SNP <- lapply(ldf3, f_fst_no_window) 

LE_df_SNP <- bind_rows(LE_list_SNP, .id ='sample')  

LE_df_SNP <-  LE_df_SNP %>%
  mutate(time = str_sub(sample, start = 14, end =-3))

LE_df_SNP %>%
  dplyr::select(CHROM, POSITION) %>%
  distinct() %>%
  nrow()

##### Get a list of SNPs with high FST ####
SNP_df <- LE_df_SNP %>%
  arrange(CHROM, POSITION) %>%
  group_by(CHROM, POSITION ) %>%
  count() %>%
  filter(n == 17)

dim(SNP_df)

high_fst_SNPs <- LE_df_SNP %>%
  right_join(SNP_df, by = c("CHROM", "POSITION")) %>%
  arrange(desc(fst))

head(high_fst_SNPs)

t <- high_fst_SNPs %>%
  filter(fst >= 0.2)

dim(t)

select_SNPs <- high_fst_SNPs[1:352883,] %>%
  arrange(CHROM, POSITION)
  
select_SNPs %>%
  group_by(CHROM, POSITION) %>%
  filter(sex == 'f') %>%
  count() %>%
  arrange(desc(n)) 

# This generates the list of SNPs used for SmLE allele frequencies
select_SNPs$interval <- paste0(select_SNPs$CHROM, ":", 
                                    select_SNPs$POSITION, '-', select_SNPs$POSITION)

SNP_list <- unique(select_SNPs$interval)

write_delim(data.frame(SNP_list), 'select_LE_SNPs.list')

#---------------------------------#
## Figure 5b separated by sex ####
#---------------------------------#

# Import all sex specific SmLE variants and merge with FST table
LE_var_m <- read.table("./LE_m_var.table", sep = '\t',header = TRUE)
LE_var_f <- read.table("./LE_f_var.table", sep = '\t',header = TRUE)

# This merges the list to retain NA values and assing sex 
LE_m_merged <- lapply(LE_list_SNP[9], function(x) LE_var_m %>% left_join(x, by = c('CHROM', "POS" = "POSITION"))) %>% na.omit()
LE_f_merged <- lapply(LE_list_SNP[1], function(x) LE_var_f %>% left_join(x, by = c('CHROM', "POS" = 'POSITION')))

LE_m_merged <- lapply(LE_m_merged, function(x) x %>%
                           mutate(fst = if_else(is.na(fst), 0, fst))) 

LE_f_merged <- lapply(LE_f_merged, function(x) x %>%
                        mutate(fst = if_else(is.na(fst), 0, fst))) 

LE_m_merged <- lapply(LE_m_merged, function(x) x %>% mutate(sex = "m"))
LE_f_merged <- lapply(LE_f_merged, function(x) x %>% mutate(sex = "f"))

LE_df_m <- bind_rows(LE_m_merged, .id ='sample')  
LE_df_f <- bind_rows(LE_f_merged, .id = 'sample')

LE_df_merged <- bind_rows(LE_df_m, LE_df_f)
LE_df_merged %>% group_by(sex) %>% count()

# A tibble: 2 × 2
# Groups:   sex [2]
#sex        n
#<chr>  <int>
#1 f     701935
#2 m     700375

## Fst histograms  ####
p_full_histo_fst <- ggplot(LE_df_merged, aes(fst, fill = sex)) +
  geom_histogram(position = position_dodge(0.01)) + 
  #facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  labs(y = 'Count', x = bquote(bold(F[ST])), fill = "Sex") +
  theme_minimal() +
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title= element_text(size = 18, face = 'bold'), panel.grid = element_blank())

# Export with LE_af.R ####
