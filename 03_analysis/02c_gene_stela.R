library(tidyverse)
library(reshape2)
detach("package:dplyr", character.only = TRUE)
library("dplyr", character.only = TRUE)

sample_pairs <- read_csv("/path/to/project/sample_pairs.csv")
sample_pairs <- sample_pairs %>% add_column(sample = paste(sample_pairs$normal_db, sample_pairs$tumor_db, sep = "_"))
sample_pairs <- sample_pairs %>% mutate(diff = (normal_stela - tumor_stela))
sample_pairs

sample_genes <- read_csv("/path/to/project/data/sample_gene_table.csv")
sample_genes

hist(sample_pairs$tumor_stela)

sample_pairs$bin <- as.integer(cut(sample_pairs$tumor_stela, breaks = seq(2, 8, by = 1), labels = 2:7))+1
#sample_pairs

gene_bin_table <- tibble(gene=sample_genes$gene)
for (i in seq(2,7,by=1)){
  gene_bin_table[[as.character(i)]] <- integer(length(gene_bin_table$gene))
}

gene_bin_table

for (i in colnames(sample_genes)[c(-1,-2,-length(colnames(sample_genes)))]) {
  bin <- sample_pairs %>% filter(sample == i) %>% select(bin) %>% pull()
  gene_bin_table[as.character(bin)] <- gene_bin_table[as.character(bin)] + sample_genes[i]
}

gene_filter_thresh <-40

gene_bin_table <- gene_bin_table[rowSums(gene_bin_table[,-1]) > gene_filter_thresh,]

gene_bin_table

hist(lapply(gene_bin_table, is.numeric))

melted_bin <- melt(gene_bin_table, id="gene")

melted_bin %>%
  #mutate(sample=factor(sample, levels=unique(sample))) %>%
  ggplot(aes(x=variable, y=value, fill=gene)) +
  geom_bar(stat = 'identity') +
  xlab("Tumour telomere length (floor)") +
  ylab("Gene count") #+
  #theme(legend.position = "none")

#as.integer(sample_pairs$bin)+1


(sample_genes %>% filter(gene == "BRCA1") %>% select("count"))
(sample_genes %>% filter(gene == "BRCA2") %>% select("count"))
(sample_genes %>% filter(gene == "TP53") %>% select("count"))

#sample_genes %>% filter(gene == "NCHL1") # not present


cancer_vs_total <- read_csv("/path/to/project/data/cancer_total_10.csv")
cancer_non <- cancer_vs_total %>% mutate(non = (total-cancer)) #%>% filter(!sample %in% c("DB170_DB171", "DB143_DB144"))

gcount_tlen <- left_join(cancer_non, sample_pairs)

cor.test(gcount_tlen$total, gcount_tlen$tumor_stela)
cor.test(gcount_tlen$total, gcount_tlen$diff)
cor.test(gcount_tlen$total, abs(gcount_tlen$diff))

gcount_tlen <- gcount_tlen %>% dplyr::rename("num_cancer_mut" = "cancer", "num_total_mut" = "total", "num_non_cancer_mut" = "non")

cancer_non <- melt(left_join(cancer_non, sample_pairs) %>% select(tumor_db, cancer, non, tumor_stela, diff), id=c('tumor_db', 'tumor_stela', 'diff'))



zoom_len <- cancer_non %>%
  arrange((tumor_stela)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("Sample ordered by tumour telomere length") + 
  ylab("Number of mutations") +
  labs(fill = "") +
  scale_fill_discrete(labels = c("tumour", "normal")) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0,10000, by =2500), seq(10000,50000, by =5000)) ) +
  theme_foundation(base_size=14, base_family="helvetica") +
  coord_cartesian(ylim=c(0,10000)) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
        axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        
        plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        #axis.title.y = element_text(angle=90,vjust =2),
        #axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major =element_blank(),  #element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(), #element_line(colour="#f0f0f0", size = 1),#
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        #legend.margin = unit(0, "cm"),
        #legend.title = element_text(face="italic"),,
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))


out_len <- cancer_non %>%
  arrange((tumor_stela)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("") + 
  ylab("") +
  labs(fill = "") +
  scale_fill_discrete(labels = c("cancer", "non-cancer")) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0,10000, by =2500), seq(10000,50000, by =5000)) ) +
  theme_foundation(base_size=14, base_family="helvetica") +
  guides(color=guide_legend(override.aes=list(fill=NA, colour = NA))) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "white"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "white"),
        legend.text = element_text(size = 15, colour = "white"),
        
        plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),

        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"))


grid.arrange(zoom_len, out_len, nrow=1, widths=c(0.7,0.3))




zoom_diff <- cancer_non %>%
  arrange((diff)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("Sample ordered by differnce in telomere length (normal-tumour)") + 
  ylab("Number of mutations") +
  labs(fill = "") +
  scale_fill_discrete(labels = c("tumour", "normal")) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0,10000, by =2500), seq(10000,50000, by =5000)) ) +
  theme_foundation(base_size=14, base_family="helvetica") +
  coord_cartesian(ylim=c(0,10000)) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
        axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        
        plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        #axis.title.y = element_text(angle=90,vjust =2),
        #axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major =element_blank(),  #element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(), #element_line(colour="#f0f0f0", size = 1),#
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        #legend.margin = unit(0, "cm"),
        #legend.title = element_text(face="italic"),,
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))


out_diff <- cancer_non %>%
  arrange((diff)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("") + 
  ylab("") +
  labs(fill = "") +
  scale_fill_discrete(labels = c("cancer", "non-cancer")) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0,10000, by =2500), seq(10000,50000, by =5000)) ) +
  theme_foundation(base_size=14, base_family="helvetica") +
  guides(color=guide_legend(override.aes=list(fill=NA, colour = NA))) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "white"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "white"),
        legend.text = element_text(size = 15, colour = "white"),
        
        plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),
        
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"))


grid.arrange(zoom_diff, out_diff, nrow=1, widths=c(0.7,0.3))




cancer_non %>%
  arrange((diff)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("Sample (low to high difference in telomere length (normal - cancer))") + 
  ylab("Number of mutations") +
  labs(fill = "Gene type") +
  scale_fill_discrete(labels = c("cancer", "non-cancer")) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(0,10000, by =2500), seq(10000,50000, by =5000)) ) +
  theme_foundation(base_size=14, base_family="helvetica") +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        
        plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        #axis.title.y = element_text(angle=90,vjust =2),
        #axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(), #element_line(colour="#f0f0f0", size = 1),#
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        #legend.margin = unit(0, "cm"),
        #legend.title = element_text(face="italic"),,
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))


mean((cancer_vs_total$cancer/cancer_vs_total$total)*100)


lower <- min(cancer_non$value-100) 
upper <- max(cancer_non$value)


cancer_non %>%
  arrange((tumor_stela)) %>% 
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  mutate(value = log10(value)) %>%
  ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab("Sample (low to high tumour telomere length)") + 
  ylab("Number of mutations") +
  labs(fill = "Gene type") +
  scale_fill_discrete(labels = c("cancer", "non-cancer")) +
  #log scale not clear
  scale_y_continuous(expand = c(0,0),
                     breaks = log10(c(seq(0,2500, by=500), seq(2500,10000, by =2500), seq(10000,30000, by =5000) )), #trans_breaks('log10', function(x) 10^x, n =5),#c( seq(2500,10000, by =2500), seq(10000,50000, by =5000)),
                     labels = c(seq(0,2500, by=500), seq(2500,10000, by =2500), seq(10000,30000, by =5000)) ) + #trans_format('log10', math_format(10^.x)))+
  #log10(c( seq(2500,10000, by =2500), seq(10000,50000, by =5000))) ) +#, labels = c( seq(2500,10000, by =2500), seq(10000,50000, by =5000)) )  +
  #annotation_logticks(base = 10, sides = "l", scaled = FALSE) +
  #scale_y_log10( breaks = seq(0,45000, by=5000) ) +
  #ylim(log10(c(lower, upper))) +
  coord_cartesian(ylim=log10(c(lower, upper))) +
  theme_foundation(base_size=14, base_family="helvetica") +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        #axis.text.y = element_text(size = c(15, 15, 15, 15, 15, 0, 0, 15,0,15,0,15,0)),
        #axis.ticks.y = element_line(size = c(1,1,1,1,1,0,0,1,0,1,0,1,0,0)),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1)),
        #axis.title.y = element_text(angle=90,vjust =2),
        #axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(), #element_line(colour="#f0f0f0", size = 1),#
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        #legend.margin = unit(0, "cm"),
        #legend.title = element_text(face="italic"),,
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))


