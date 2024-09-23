# Selection coefficient 
# Author: Kathrin Jutzeler
# Date: February 19 , 2023
# Last updated: August 6, 2024
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(GenomicRanges)
library(magrittr)
library(tidyverse)
library("ggpubr")
library(fields)
library(locfit)
library(patchwork)

#--------------------#
# Import the data ####
#--------------------#

LE_AF <- read.table("LE_specific.filter.table", sep = '\t',header = TRUE)

dim(LE_AF)
#[1] 96485   80


#LE_AF <- filter(LE_AF, CHROM == 'SM_V10_1')

# Calculate (pooled) allele frequencies
AD <- LE_AF[,seq(5,80,4)]
DP <- LE_AF[,seq(6,80,4)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)

#[1] 96485  19  

# Make and organize data frame 

# refFre.AD <- matrix(ncol=19,nrow=96485)
# for(j in 1:19){
# 	refFre.AD[,j]<- sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
# 	}
# 
# refFre.AD[DP<20]<-NA
# 
# colnames(refFre.AD)<-colnames(AD)	
# refFre.AD <- cbind(LE_AF[,1:4],refFre.AD)
# colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

SNPFre.AD <- matrix(ncol=19,nrow=96485)
for(j in 1:19){
  SNPFre.AD[,j]<- 1 - (sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j]))
}

SNPFre.AD[DP<20]<-NA

colnames(SNPFre.AD)<-colnames(AD)	
SNPFre.AD <- cbind(LE_AF[,1:4],SNPFre.AD)
colnames(SNPFre.AD) <- gsub(".AD", "", colnames(SNPFre.AD))

df <- gather(SNPFre.AD, 'sample', 'AF', 5:23)

df$pop <- str_sub(df$sample, start = 1L, end =2L)
df$time <- df$time <- gsub("\\D", "", df$sample)
df$sex <- str_sub(df$sample, start = -1)

df <- df %>%
  mutate(smooth = runmed(AF,51))
  
# Put time in actual order
df_levels <- c("110916","062718", "051420", "112320", "070521", "092921", "102621",  "122121", 
                   "070522","021523")

df$time <- factor(df$time, levels= df_levels)


#----------------------------------------------------------#
# Plot AF change in BRE across time for each chromosome ####
#----------------------------------------------------------#

mean_df <- df %>% 
  filter(!CHROM %in% c('SM_V10_WSR', 'SM_V10_MITO')) %>%
  group_by(time, CHROM) %>%
  summarize(mean = mean(AF, na.rm =T))

#Rename chroms

chrom.labs <- c('Chr. 1', 'Chr. 2', 'Chr. 3', 'Chr. 4', 'Chr. 5', 'Chr. 6', 'Chr. 7', 
                'Chr. Z')
names(chrom.labs) <- levels(factor(mean_df$CHROM))

comp_time_order <- c("Nov 2016", "June 2018", "May 2020", "Nov 2020", "July 2021", "Sep 2021", "Oct 2021", 
                     "Dec 2021","July 2022", "Feb 2023")

names(comp_time_order) <- levels(factor(mean_df$time))
  
mean_AF <- ggplot(mean_df, aes(time, mean, group = CHROM)) +
  geom_point() +
  geom_line() +
  geom_text(data =subset(mean_df, time == '021523'), aes(label = round(mean,2)),
            nudge_x = 0, nudge_y = -0.2) +
  #geom_smooth(size =1, se = F) +
  ggh4x::facet_wrap2(CHROM~., labeller = labeller(CHROM = chrom.labs), axes = 'all',
                     remove_labels = 'all', ncol = 2) +
  theme_minimal() +
  scale_x_discrete(labels = comp_time_order) +
  labs(x= "Time", y = "Mean allele frequency")  +
  theme(axis.text.x = element_text(angle =90, hjust = 1), axis.line = element_line(), 
        legend.position = 'none',
        text = element_text(size =16), panel.grid = element_blank(), 
        axis.title = element_text(face =  'bold'), axis.ticks = element_line()) 

ggsave('Figure4a.jpg', width = 8, height =6, units = 'in', dpi =600, bg="white")


WSR <- SNPFre.AD %>%
  filter(CHROM == 'SM_V10_WSR')

dim(WSR)

WSR <- WSR %>%
  gather('sample', 'AF', 5:23) 

WSR$time <- WSR$time <- gsub("\\D", "", WSR$sample)

WSR$time <- factor(WSR$time, levels = df_levels)

WSR %>%
  arrange(POS, time)

#---------------------------#
# Selection coefficient  ####
#---------------------------#

scInfo <- read.csv('conta_metadata.csv') 
scInfo <- scInfo %>% filter(population == 'BRE')

dim(scInfo)
#[1] 7  5

AF.BRE <- subset(refFre.AD, select= as.character(scInfo$sampleID))
rownames(AF.BRE) <- paste(SNPFre.AD$CHROM, ".", SNPFre.AD$POS,sep="")

dim(AF.BRE)
#[1] 96485   7

Ln.ratio <- log(AF.BRE/(1-AF.BRE))
Ln.ratio <- cbind(scInfo,t(Ln.ratio))
is.na(Ln.ratio) <- sapply(Ln.ratio, is.infinite)
dim(Ln.ratio)
#[1]   9  96490


#write.csv(cbind(scInfo,t(AF.BRE)), file = "AF.BRE.csv",row.names=FALSE)
#write.csv(Ln.ratio, file = "Ln.ratio.csv",row.names=FALSE)

#Ln.ratio <- read.csv('Ln.ratio.csv')

Ln.ratio[1:2,1:10]

r <- Ln.ratio %>%
  gather('position', 's', 6:96490)

r2 <- r %>% separate(position, into = c('CHROM', 'POS'), sep = '\\.')

head(r2)

#-------------------#
# Plot ln ratio ####
#-------------------#

r_mean <- r2 %>%
  filter(!sampleID %in% c('BRE092921_f', 'BRE092921_m'), !CHROM %in% c('SM_V10_WSR', 'SM_V10_MITO')) %>%
  group_by(LCS) %>%
  reframe(mean_s = mean(s, na.rm = T))

#jpeg("Figure3b.jpg", width = 3, height = 3, units = 'in', res =300)  
  
p_ratio <- ggplot(r_mean, aes(LCS, -mean_s)) +
  #geom_point(color = "#00A08A", size =2) +
  geom_point(size =2) +
  geom_smooth(method = 'lm', se =F, color ='red') +
  stat_regline_equation(label.y = 1.8, label.x = 0) +
  stat_cor(method = 'pearson',  digits =3, #label.y = 2e+5, label.x = 0.01,
           aes(label = paste(after_stat(rr.label), sep = "~`,`~"))) +
  #facet_wrap(~ CHROM, scales ='free') +
  #scale_x_continuous(labels=c("092921", "102621","122121",  "070522",  "021523")) +
  #scale_color_manual(values= conta_colors) +
  labs(x = "Sexual life cycles", y = "Ln(genotype ratio)") +
  theme_bw() +
  theme(text = element_text(size =16), axis.title = element_text(face = 'bold'),
        panel.grid = element_blank()) 

ggsave('Figure3b.jpg', width = 4, height =4, units = 'in', dpi =600, bg="white")

#dev.off()

#ggsave('ln_ratio_slope.png', width =8, heigh =6)

# Calculate correlation coefficient ####

BRE <- Ln.ratio[Ln.ratio$population == "BRE",]
#BRE <- BRE[-1:-2,] 

#BRE$LCS <- as.integer(c('0', '0', '1', '1', '5', '10', '10'))

dim(Ln.ratio)

dim(BRE)

# Make data frame
sc <- matrix(ncol=1,nrow=96485)
rownames(sc) <- rownames(AF.BRE)
colnames(sc) <- c('BRE')

# BRE
for(j in 6:96490){
  ifelse(
    all(is.na(BRE[,j])), 
    sc[j-4,1] <- NA, 
    sc[j-4,1] <- -coef(lm( as.numeric(BRE[,j]) ~ BRE$LCS))[2]
  )	  
}

sc2 <- cbind(refFre.AD[,1:4],sc)

sc3 <- sc2 %>%
  group_by(CHROM) %>%
  mutate(smooth_BRE = runmed(BRE, 51)) %>%
  arrange(CHROM, POS)

# Plot selection coefficient across genome ####
chrom.labs <- c('1', '2', '3', '4', '5', '6', '7', 'Z')
names(chrom.labs) <- levels(factor(mean_df$CHROM))

p_selection <- ggplot(subset(sc3, !CHROM %in% c("SM_V10_WSR", 'SM_V10_MITO')), aes(POS, smooth_BRE)) +
  #geom_line() + 
  geom_point(color = 'black', size = 0.5,alpha = .1) +
  geom_smooth(method='locfit', method.args = list(deg=0, alpha=0.3),
            linewidth =1, se = F, color = "red") +
  
  facet_grid(~CHROM, switch = "x", scales = "free_x", space = 'free_x', labeller = labeller(CHROM= chrom.labs)) +
  theme_minimal() +
  ylim(0,0.45) +
  geom_hline(yintercept = 0.229, linetype = 'dashed', linewidth =1) + # shows mean, but median = 0.228
  labs(x= "Chromosome", y = expression(italic("S"))) +
  theme(axis.text.x = element_blank(), panel.grid.minor = element_blank(), legend.position = 'bottom',
        axis.line.x = element_line(), axis.line.y = element_line(),
        text = element_text(size =16), strip.text.x = element_text( size = 12), strip.placement = "outside",
        panel.grid = element_blank(), axis.title = element_text(face = 'bold'))

ggsave('Figure3c.jpg', width = 10, height =6, units = 'in', dpi =600, bg="white")

#-----------------#
# Export plots ####
#-----------------#

jpeg('Figure4.jpg', width = 10, height =12, units = 'in', res = 300)

#(mean_AF | p_ratio) / p_selection
((mean_AF | p_ratio) + plot_layout(widths = c(6, 4))) / p_selection +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()


# Find genes of interest (test section) ####
t <- sc3 %>%
  filter(CHROM == 'SM_V10_2', smooth_BRE > 0.3) #%>%
#  summarize(max = max(smooth_BRE, na.rm = T))

t$chr <- t$CHROM
t$start <- t$POS
t$end <- t$POS

t <- t[,7:9]

sc3 %>%
  filter(CHROM == "SM_V10_Z") %>%
  summary()

lines <- which(sc3$smooth_BRE >= 0.365, sc3$CHROM)

sc3[lines,] %>%
  filter(CHROM == 'SM_V10_Z')

library(GenomicRanges)

#gff <- rtracklayer::import.gff('schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3')

gt_t <- makeGRangesListFromDataFrame(t)

overlaps <- findOverlaps(gt_t, gff)
overlapping_genes <- subjectHits(overlaps)

overlapping_genes <- unique(subjectHits(overlaps))

# Retrieve gene information
overlapping_genes_info <- as.data.frame(mcols(gff)$ID[overlapping_genes])

fil <- grep(pattern = '^gene', overlapping_genes_info$`mcols(gff)$ID[overlapping_genes]`)

gen <- data.frame(gsub('gene:', '', overlapping_genes_info[fil,]))
