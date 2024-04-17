# Allele frequencis for Contamination project
# Author: Kathrin Jutzeler
# Date: February 19 , 2023
# Last updated: April 17, 2024
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(magrittr)
library(tidyverse)
library(ggpubr)
library(fields)
library(locfit)

#LE_col <- "#00A08A"

# Import AF of SNPs with high Fst ####
high_fst <- read.table("./LE_AF/LE_high_Fst.table", sep = '\t',header = TRUE)

dim(high_fst)
#[1] 246   43

# Extract AD and DP
AD <- high_fst[,seq(6,43,2)]
DP <- high_fst[,seq(7,43,2)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)

#[1] 246  19  

# Get allele frequencies 
refFre.AD <- matrix(ncol=19,nrow=246)
for(j in 1:19){
  refFre.AD[,j]<- sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
}

# Drop all rows with DP < 20
refFre.AD[DP<20]<-NA

# Make data frame
colnames(refFre.AD) <- colnames(AD)	
refFre.AD <- cbind(high_fst[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

dim(refFre.AD)
#[1] 246  23

# Reorganize data frame
high_fst_df <- gather(refFre.AD, 'sample', 'AF', 5:23)

high_fst_df <- high_fst_df %>%
  mutate(time = str_sub(sample, start = 3, end =8), sex =str_sub(high_fst_df$sample, start = 10, end =10))

LE_levels <- c("110216","062118", "051420", "112320", "070521", "092921", "102621",  "122121", 
                   "070522","021523")


# AF change across time ####
high_fst_df$time <- factor(high_fst_df$time, levels = LE_levels)

# Isolate SNPs with the highest change
change_by_time <- high_fst_df %>%
  #filter( AF == '0') %>%
  arrange(CHROM, POS, sex, time) %>%
  group_by(CHROM, POS,sex) %>%
  reframe(change = abs(diff(AF)))

## Plot histogram of change over time ####
p_change_over_time <-
ggplot(change_by_time, aes(change, fill = sex)) +
  geom_histogram(position = position_dodge(0.02)) + 
#  facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  theme_minimal() + 
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title= element_text(size = 18, face = 'bold'), panel.grid = element_blank()) +
  labs(y = 'Count',  x = 'Consecutive allele frequency change', fill = "Sex")

ggsave('change_by_time.pdf', width = 8, height = 6)

# AF change overall ####
# Isolate SNPs with the highest change
change_overall <- high_fst_df %>%
  #filter( AF == '0') %>%
  arrange(CHROM, POS, sex, time) %>%
  group_by(CHROM, POS, sex) %>%
  reframe(min = min(AF), max = max(AF)) %>%
  mutate(change = max - min)

## Plot histogram of change overall ####
p_change_overall <- 
ggplot(subset(change_overall, change != 0), aes(change, fill = sex)) +
  geom_histogram(position = position_dodge(0.02)) + 
  #  facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  theme_minimal() + 
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title= element_text(size = 18, face = 'bold'), panel.grid = element_blank()) +
  labs(y = 'Count',  x = 'Allele frequency change overall', fill = "Sex")

#ggsave('change_overall.pdf', width = 8, height = 6)

jpeg('Figure4c.jpg', width = 10, height =3, units = 'in', res =300)
(p_change_over_time | p_change_overall) + plot_layout(guides = "collect") 
dev.off()

# Look for SNPs with maximum allele frequency change
#max_change <- change_overall %>%
#  arrange(desc(change), CHROM, POS)

#max_change_sel <- max_change[c(2,7,22,27,35,38),]

# AF for each SNP across time ####

# Check
by_time <- high_fst_df %>%
  arrange(CHROM, POS, time) 

by_time %>%
  group_by(CHROM, POS) %>%
  count() %>%
  filter(n != 19)


# 
# v <- t[seq(1, length(t), 19)]
# 
# ggplot(by_time[c(3288:3306),], aes(factor(time,levels = LE_levels), AF, group = interaction(CHROM, POS, sex),
#                                   color =  interaction(CHROM, POS))) +
#   geom_line(aes(linetype = sex)) + 
#   annotate("rect", xmin='092921', xmax='102621', ymin=0, ymax=1, alpha=0.2, fill="grey") +
#   ylim(0,1) +
#   labs(x= "Time", y = "AF", color = 'SNP', linetype = 'Sex') +
#   theme_minimal() 

# 19 entries for each SNP, 5 SNPs per plot = 95
v <- seq(1,4674,19)

v2 <- v[c(3,12,17,66,75,108,124,176,195,197,213,236)]

t <- list()

## Plot AF change for select SNPs ####
for (i in v2){
  t[[i]] <- by_time[c(i:(i+18)),]
  
}

t <- t %>% discard(is.null)

t_df <- bind_rows(t)

int_levels <- t_df %>% arrange(CHROM, POS) %>% distinct(CHROM, POS) %>% mutate(int = factor(paste0(CHROM,'.',POS)))
  
#t <- by_time %>%
 # filter(CHROM == 'SM_V10_5', POS == 11368169)

# Plot AF across time for select SNPs 
ggplot(t_df, aes(factor(time,levels = LE_levels), AF, group = interaction(CHROM, POS, sex))) +
  geom_line(aes(linetype = sex), linewidth = 1) + 
  annotate("rect", xmin='092921', xmax='102621', ymin=0, ymax=1, alpha=0.2, fill="grey") +
  facet_wrap(~ factor(interaction(CHROM, POS), levels = int_levels$int)) +
  ggh4x::facet_wrap2(~ factor(interaction(CHROM, POS), levels = int_levels$int), axes = 'all',
                     remove_labels = 'all') +
  #scale_color_manual(values = "#b0d4a4" ) +
  ylim(0,1) +
  labs(x= "Time", y = "AF", color = 'SNP', linetype = 'Sex') +
  theme_minimal() +
  theme(axis.line = element_line(), axis.text.x = element_text(size = 12, angle =45, hjust =1), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 14),
        panel.grid = element_blank(), axis.ticks = element_line())

ggsave('Figure5a.jpg', width = 10, height =6, units = 'in', dpi =600, bg="white")

# Plot AF change for each SNP
pdf('temp.pdf', width = 8, height =6)
for (i in v) {
  p <- ggplot(by_time[c(i:(i+18)),], aes(factor(time,levels = LE_levels), AF, group = interaction(CHROM, POS, sex),
                                         color =  interaction(CHROM, POS))) +
    geom_line(aes(linetype = sex)) + 
    annotate("rect", xmin='092921', xmax='102621', ymin=0, ymax=1, alpha=0.2, fill="grey") +
    scale_color_manual(values = "#b0d4a4" ) +
    ylim(0,1) +
    labs(x= "Time", y = "AF", color = 'SNP', linetype = 'Sex') +
    theme_minimal() +
    theme(axis.line = element_line())
   
   print(p)
}

dev.off()


# Export plots ####

jpeg('Figure5d.jpg', width = 10, height =3, units = 'in', res =300)
(p_change_over_time | p_change_overall)  + plot_layout(guides = "collect") 
#plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))
dev.off()

pdf('Figure_LE_AF_change.pdf', width = 8, height = 10)

p_change_over_time / p_change_overall + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A")) &  
  theme(plot.tag = element_text(face = 'bold', size = 18), 
        legend.position = 'bottom', panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'))

dev.off()

