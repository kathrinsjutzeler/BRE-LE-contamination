# Allele frequencis for Contamination project
# Author: Kathrin Jutzeler
# Date: February 19 , 2023
# Last updated: August 6, 2024
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(magrittr)
library(tidyverse)
library(ggpubr)
library(fields)
library(locfit)
library(patchwork)

#---------------#
# Figure 5C ####
#---------------#

# Function to calculate allele frequencies from AD and DP
f_AF <- function(df, var){
  # Extract AD and DP
  AD <- df[,seq(7,10,2)]
  DP <- df[,seq(8,10,2)]
  ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
  
  # Get allele frequencies 
  refFre.AD <- matrix(ncol=2,nrow=var)
    for(j in 1:2){
  refFre.AD[,j]<- 1 - (sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j]))
  }

  # Drop all rows with DP < 10
  #refFre.AD[DP<10]<-NA

  # Make data frame
  colnames(refFre.AD) <- colnames(AD)	
  refFre.AD <- cbind(df[,1:4],refFre.AD)
  colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

  return(refFre.AD)
}

# Run function on data frames generated in conta_fst.R
dim(LE_df_m)
AF_m <- f_AF(LE_df_m,661996)

dim(LE_df_f)
AF_f <- f_AF(LE_df_f, 657592)

# Reorganize and merge data frames
AF_m <- gather(AF_m, 'sample', 'AF', 5:6)
AF_f <- gather(AF_f, 'sample', 'AF', 5:6)

AF <- bind_rows(AF_m, AF_f)

fst_df <- AF %>%
  mutate(time = str_sub(sample, start = 3, end =8), sex =str_sub(sample, start = 10, end =10))

## AF change overall ####
# Isolate SNPs with the highest change
change_overall <- fst_df %>%
  #filter( AF == '0') %>%
  arrange(CHROM, POS, time) %>%
  group_by(CHROM, POS, sex) %>%
  reframe(change = abs(diff(AF)))
  
### Plot histogram of change overall ####
p_change_overall <- 
  ggplot(change_overall, aes(change, fill = sex)) +
  geom_histogram(position = position_dodge(0.02)) + 
  #  facet_wrap(~ sex) +
  scale_fill_manual(values = my_color) +
  theme_minimal() + 
  theme(axis.line = element_line(), text = element_text(size =16), 
        axis.title= element_text(size = 18, face = 'bold'), panel.grid = element_blank()) +
  labs(y = 'Count',  x = 'Allele frequency change', fill = "Sex")

change_overall %>%
  group_by(sex) %>%
  summarize(mean = mean(change, na.rm = T), min = min(change, na.rm = T), max = max(change, na.rm = T))


change_overall %>% group_by(sex) %>% count()

# Number of variants with AF change > 0.8
change_overall %>% 
  group_by(sex) %>%
  filter(change > 0.8) %>%
  count()
#2216
2216/706496 *100

# Number of variants with AF change = 1
change_overall %>% 
  group_by(sex) %>%
  filter(change == 1) %>%
  count()
#583
583/706496 *100

#--------------#
# Figure 5A ####
#--------------#

## AF change across time ####
# Import AF table with all samples
high_fst <- read.table("./LE_AF/LE_high_Fst.table", sep = '\t',header = TRUE)

## Run AF function ####
dim(high_fst)
#246 43

AD <- high_fst[,seq(6,43,2)]
DP <- high_fst[,seq(7,43,2)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}

dim(AD)
#246  19
# Get allele frequencies 
refFre.AD <- matrix(ncol=19,nrow=246)
for(j in 1:19){
  refFre.AD[,j]<- 1 - (sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j]))
}

# Make data frame
colnames(refFre.AD) <- colnames(AD)	
refFre.AD <- cbind(high_fst[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

# Reorganize data frame
high_fst_df <- gather(refFre.AD, 'sample', 'AF', 5:23)

high_fst_df <- high_fst_df %>%
  mutate(time = str_sub(sample, start = 3, end =8), sex =str_sub(high_fst_df$sample, start = 10, end =10))

high_fst_df <- high_fst_df %>% filter(CHROM != "SM_V10_WSR")

# Check
by_time <- high_fst_df %>%
  arrange(CHROM, POS, time) 

by_time %>%
  group_by(CHROM, POS) %>%
  count() %>%
  filter(n != 19)

# 19 entries for each SNP, 5 SNPs per plot = 95
v <- seq(1,4674,19)

# Select SNPs for plot
v2 <- v[c(3,12,17,66,75,108,124,176,195,197,213,236)]

t <- list()

## AF change for select SNPs ####
for (i in v2){
  t[[i]] <- by_time[c(i:(i+18)),]
  
}

t <- t %>% discard(is.null)

t_df <- bind_rows(t) 
t_df <- t_df %>% mutate(CHROM = str_remove(CHROM, "SM_V10_"))

int_levels <- t_df %>% arrange(CHROM, POS) %>% distinct(CHROM, POS) %>% 
  mutate(int = factor(paste0(CHROM,'.',POS)))

comp_time_order <- c("Nov 2016", "June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                     "Dec 2021","July 2022", "Feb 2023")

names(comp_time_order) <- levels(factor(t_df$time, levels = LE_levels))

## Plot AF across time for select SNPs ####
p_AF_time <- ggplot(t_df, aes(factor(time,levels = LE_levels), AF, group = interaction(CHROM, POS, sex))) +
  geom_line(aes(linetype = sex), linewidth = 1) + 
  ggh4x::facet_wrap2(~ factor(interaction(CHROM, POS), levels = int_levels$int), axes = 'all',
                     remove_labels = 'all') +
  #scale_color_manual(values = "#b0d4a4" ) +
  ylim(0,1) +
  scale_x_discrete(labels = comp_time_order) +
  labs(x= "Time", y = "Allele frequency", color = 'SNP', linetype = 'Sex') +
  theme_minimal() +
  theme(axis.line = element_line(), axis.text.x = element_text(size = 12, angle =45, hjust =1), 
        axis.title = element_text(size = 18, face = 'bold'), text = element_text(size = 14),
        panel.grid = element_blank(), axis.ticks = element_line())

ggsave('Figure5a.jpg', width = 10, height =6, units = 'in', dpi =600, bg="white")

# # Plot AF change for each SNP
# pdf('temp.pdf', width = 8, height =6)
# for (i in v) {
#   p <- ggplot(by_time[c(i:(i+18)),], aes(factor(time,levels = LE_levels), AF, group = interaction(CHROM, POS, sex),
#                                          color =  interaction(CHROM, POS))) +
#     geom_line(aes(linetype = sex)) + 
#     annotate("rect", xmin='092921', xmax='102621', ymin=0, ymax=1, alpha=0.2, fill="grey") +
#     scale_color_manual(values = "#b0d4a4" ) +
#     ylim(0,1) +
#     labs(x= "Time", y = "AF", color = 'SNP', linetype = 'Sex') +
#     theme_minimal() +
#     theme(axis.line = element_line())
#    
#    print(p)
# }
# 
# dev.off()

#-----------------#
# Export plots ####
#-----------------#

jpeg('Figure5.jpg', width = 10, height =9, units = 'in', res =300)
  p_AF_time /
  ((p_full_histo_fst | p_change_overall) + plot_layout(guides = "collect")) +
    plot_layout(height = c(6, 3)) +
    plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))
dev.off()

