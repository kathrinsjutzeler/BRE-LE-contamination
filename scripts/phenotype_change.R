# Cercarial shedding 
# Author: Kathrin Jutzeler
# Date: May 15, 2023
# Last updated: April 17, 2024
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(patchwork)



############ Cercarial shedding ######################
# Import data
BRE15 <- read_excel('BRE_LE_pheno.xlsx', sheet = 'BRE_2015')
BRE23 <- read_excel('BRE_LE_pheno.xlsx', sheet = 'BRE_2023')

BRE15 <- BRE15 %>%
  mutate(year = "2015") %>%
  mutate_at(c('shed1', 'shed2', 'shed3', 'shed4'), as.numeric)

BRE23 <- BRE23 %>%
  mutate(year = "2023") %>%
  mutate_at(c('shed1', 'shed2', 'shed3', 'shed4'), as.numeric)

BRE_df <- bind_rows(BRE15, BRE23)

BRE_df <- BRE_df %>%
  gather('shed', 'n', 3:6)


# Calculate stats
signif <- BRE_df %>%
  group_by(shed) %>%
  wilcox_test(n ~ year) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(p.col = "p.adj")

t_fold <- BRE_df %>% 
  group_by(year, shed) %>%
  summarize(mean = mean(n, na.rm =T)) 
  
t_fold %>%
  spread(year, mean) %>%
  mutate(fold = `2023` / `2015`)

# Plot the data
p_shed <- ggplot(BRE_df, aes(shed, n, fill = year)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Week post infection", y = "Cercarial production", fill = "Year") + 
  theme( text = element_text(size =18), legend.position = "none",
         axis.line = element_line(), panel.grid = element_blank(), 
         axis.title = element_text(face = 'bold')) + 
    scale_x_discrete(labels=c("shed1" = "4", "shed2" = "5", 'shed3' = '6', 'shed4' = '7')) +
  scale_fill_manual(values = c('white', "gray" )) +
  #scale_fill_manual(values = c('#b0d4a4', "white" )) +
  #scale_color_manual(values = c('black', 'white')) +
  guides(color = FALSE) +
  stat_compare_means(method = "wilcox", label = 'p.signif', label.y = 8000, size =8)

#ggsave('Plots/pheno_BRE.png', width = 8, height = 6)


############ Worm burden ######################

# Import data
worms <- read_excel('phenotype.xlsx')

worms_df <- worms %>%
  filter(`Sm strain` == 'SmBRE', `Rodent type` == "Hamster") %>%
  mutate(year = ifelse(`Dissection date` < '2016-01-01' & `Dissection date` > '2014-12-31', '2015', 
                       ifelse(`Dissection date` > '2022-12-31',"2023", "other")),
         total = `Number of females worms` + `Number of male worms`) %>%
  filter(year %in% c("2015", "2023"))

# Calculate stats
worms_df %>%
  group_by(year) %>%
  shapiro_test(total)

worms_df$`Number of cercariae`<- as.numeric(worms_df$`Number of cercariae`)

temp <- worms_df %>%
  group_by(year) %>%
  summarize(mean = mean(`Number of cercariae`), norm = total / `Number of cercariae`) %>%
  ungroup()

signif2 <- 
  temp %>%
  t_test(norm ~ year) %>%
  add_significance(p.col = "p")

temp %>%
  filter(norm != 'NA') %>%
  count()

temp %>%
  group_by(year) %>%
  summarize(mean = mean(norm, na.rm =T))

# Plot worms BRE ####
p_worms <- ggplot(temp, aes(year, norm, fill = year, group = year)) + # only use color to generate transition plot
  geom_boxplot() +
  scale_fill_manual(values = c('white', "gray" )) +
  #scale_fill_manual(values = c('#b0d4a4', 'white' )) +
  #scale_color_manual(values = c('black', 'white')) +
  theme_minimal() +
  labs(x = "Year", y = "Normalized worm burden") + 
  theme( text = element_text(size =18), legend.position = "none",
         axis.line = element_line(), panel.grid = element_blank(), 
         axis.title = element_text(face = 'bold')) + 
  stat_compare_means(method = "t.test", label = 'p.signif', size =8, label.x = 1.5)

ggsave('Plots/pheno2_BRE.png', width = 8, height = 6)


# Worm burden LE ####

LE_df <- worms %>%
  filter(`Sm strain` == 'SmLE', `Rodent type` == "Hamster") %>%
  mutate(year = ifelse(`Dissection date` < '2016-01-01' & `Dissection date` > '2014-12-31', '2015', 
                       ifelse(`Dissection date` > '2022-12-31',"2023", "other")),
         total = `Number of females worms` + `Number of male worms`) %>%
  filter(year %in% c("2015", "2023"))

# Calculate stats
LE_df$`Number of cercariae`<- as.numeric(LE_df$`Number of cercariae`)

LE_df <- LE_df %>%
  group_by(year) %>%
  summarize(mean = mean(`Number of cercariae`), norm = total / `Number of cercariae`)

LE_df%>%
  filter(norm != 'NA') %>%
  count()

LE_df %>%
  group_by(year) %>%
  summarize(mean = mean(norm, na.rm =T)) 

LE_df %>% 
  ungroup() %>%
  t_test(norm ~ year)

# Plot 
p_worms_LE <- ggplot(LE_df, aes(year, norm, fill = year, group = year)) + # only use color to generate transition plot
  geom_boxplot() +
  scale_fill_manual(values = c('#FFB94D', "#F98400" )) +
  #scale_fill_manual(values = c('#b0d4a4', 'white' )) +
  #scale_color_manual(values = c('black', 'white')) +
  theme_minimal() +
  labs(x = "Year", y = "Normalized worm burden") + 
  theme( text = element_text(size =18), legend.position = "none",
         axis.line = element_line()) + 
  stat_compare_means(method = "t.test", label = 'p.signif', size =6, label.x = 1.5)


# Export plots ####
jpeg('Figure1.jpg', width = 8, height = 12, unit = 'in', res = 300)

((p_shed & theme(legend.position = 'bottom'))/ p_worms) + 
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()

 
