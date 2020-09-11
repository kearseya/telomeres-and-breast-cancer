library(ltm)
library(tidyverse)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(MutationalPatterns)

library(reshape)

library(RColorBrewer)
library(colorRamps)
library(ggsci)
library(ggthemes)


#nature 2013: Signatures of mutational processes in human cancer
#looking for signatures: 1B, 2, 3, 8, 13
#https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
#looking for signatures: 2, 3, 8, 13, 17, 18, 20, 26, 30.

sample_pairs <- read_csv("/path/to/project/sample_pairs.csv")
sample_pairs

samples <- list.dirs("/path/to/project/neu-out")
samples <- samples[seq(2,length(samples))]

sig_data_list <- c()

for (i in seq_along(samples)) {
  file <- paste(samples[i], "/10_filtered.vcf", sep = "")
  sample <- basename(samples[i])
  print(sample)
  
  data <- read_tsv(file, col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",	"SAMPLE"), skip = 10)
  #data <- rename(data, "#CHROM"="CHROM")
  data <- add_column(data, sample)
  #data <- rename(data, sample="sample")
  if (nrow(data) > 0) {
    sigs.input <- mut.to.sigs.input(mut.ref = as.data.frame(data),
                                    sample = "sample",
                                    chr = "CHROM", 
                                    pos = "POS", 
                                    ref = "REF", 
                                    alt = "ALT",
                                    bsg = BSgenome.Hsapiens.UCSC.hg38)
    
    plot_data <- whichSignatures(tumor.ref = sigs.input, 
                                 signatures.ref = signatures.cosmic[-c(4,9,11,12,14,15,16,19,21,22,23,24,25,27,28,29),],
                                 contexts.needed = TRUE,
                                 sample.id = sample)
    
    sig_data_list[[i]] <- as.data.frame(plot_data$weights)
  }
  #plot_data$weights[which(plot_data$weights!=0)]
  #plotSignatures(plot_data, sub = "COSMIC")
  #makePie(plot_data, sub = "COSMIC")
}


sig_data_list
sig_data = do.call(rbind, sig_data_list)
sig_data <- rownames_to_column(sig_data, var = "sample")

list(sig_data$sample)

sample_pairs <- sample_pairs %>% add_column(sample = paste(sample_pairs$normal_db, sample_pairs$tumor_db, sep = "_"))
sig_data <- sig_data %>% left_join(sample_pairs, by="sample")
sig_data <- sig_data %>% mutate(diff = (normal_stela - tumor_stela))

barplot(sort(sig_data$diff))
ggplot(sig_data, aes(diff)) + geom_histogram()
plot(sig_data$normal_stela, sig_data$tumor_stela)


sig_data %>% arrange(diff)

sig_data <- sig_data %>% mutate(unknown = 1-rowSums(.[2:31]))
?transmute

melted <- melt.data.frame(sig_data, id.vars = c("sample", "normal_db", "tumor_db", "normal_original", "tumor_original", "normal_stela", "tumor_stela", "diff"))
melted <- melted %>% filter(value > 0)
melted <- melted %>% mutate(perc=(value*diff))
melt_sig_plot_data <- melted %>% select(c(sample, diff, variable, perc))

#nature 2013: Signatures of mutational processes in human cancer
#looking for signatures: 1B, 2, 3, 8, 13
#https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
#looking for signatures: 2, 3, (5), 8, 13, 17, 18, 20, 26, 30.

melt_sig_plot_data %>%
  arrange((diff)) %>%
  mutate(sample=factor(sample, levels=unique(sample))) %>% 
  ggplot(aes(x=sample, y=perc, fill=variable)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  xlab("") +
  ylab("Normal-Tumour telomere length (kb)")




melted <- melt.data.frame(sig_data, id.vars = c("sample", "normal_db", "tumor_db", "normal_original", "tumor_original", "normal_stela", "tumor_stela", "diff"))
melted <- melted %>% filter(value > 0)
melted <- melted %>% mutate(perc=(value*tumor_stela))
melt_sig_plot_data <- melted %>% dplyr::select(c(sample, normal_stela, tumor_stela, variable, perc))

#nature 2013: Signatures of mutational processes in human cancer
#looking for signatures: 1B, 2, 3, 8, 13
#https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
#looking for signatures: 2, 3, (5), 8, 13, 17, 18, 20, 26, 30.

melt_sig_plot_data %>%
  arrange((tumor_stela)) %>%
  mutate(sample=factor(sample, levels=unique(sample))) %>% 
  ggplot(aes(x=sample, y=perc, fill=variable)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  xlab("") +
  ylab("Tumour telomere length (kb)") +
  geom_hline(yintercept=3.81, linetype='dashed', size=0.2) +
  geom_point(aes(y=normal_stela), shape='x', alpha=0.5)

# stnadardised
melted %>%
  arrange((tumor_stela)) %>%
  mutate(sample=factor(sample, levels=unique(sample))) %>% 
  ggplot(aes(x=sample, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  xlab("Sample (low to high telomere lentgh)") +
  ylab("Signature percentage") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(13, "Accent"))(13)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())




#Map signatures to their aetiologies

aetiology_mapping <- tibble(sig = character(), aet = character()) %>%
  add_row(          sig = "Signature.1", aet = "deamination of 5mC") %>%
  add_row(          sig = "Signature.2", aet = "APOBEC") %>%
  add_row(          sig = "Signature.3", aet = "homologous recombination") %>%
  add_row(          sig = "Signature.4", aet = "tobacco") %>%
  add_row(          sig = "Signature.5", aet = "unknown") %>%
  add_row(          sig = "Signature.6", aet = "DNA mismatch repair") %>%
  add_row(          sig = "Signature.7", aet = "ultraviolet light") %>%
  add_row(          sig = "Signature.8", aet = "unknown") %>%
  add_row(          sig = "Signature.9", aet = "polymerase η") %>% 
  add_row(          sig = "Signature.10", aet = "POLE") %>%
  add_row(          sig = "Signature.11", aet = "alkylating agents") %>% 
  add_row(          sig = "Signature.12", aet = "unknown") %>%
  add_row(          sig = "Signature.13", aet = "APOBEC") %>%
  add_row(          sig = "Signature.14", aet = "unknown") %>%
  add_row(          sig = "Signature.15", aet = "DNA mismatch repair") %>% 
  add_row(          sig = "Signature.16", aet = "unknown") %>%
  add_row(          sig = "Signature.17", aet = "unknown") %>%
  add_row(          sig = "Signature.18", aet = "unknown") %>%
  add_row(          sig = "Signature.19", aet = "unknown") %>%
  add_row(          sig = "Signature.20", aet = "DNA mismatch repair") %>%
  add_row(          sig = "Signature.21", aet = "unknown") %>%
  add_row(          sig = "Signature.22", aet = "aristolochic acid") %>%
  add_row(          sig = "Signature.23", aet = "unknown") %>%
  add_row(          sig = "Signature.24", aet = "aflatoxin") %>%
  add_row(          sig = "Signature.25", aet = "unknown") %>%
  add_row(          sig = "Signature.26", aet = "DNA mismatch repair") %>%
  add_row(          sig = "Signature.27", aet = "unknown") %>%
  add_row(          sig = "Signature.28", aet = "unknown") %>%
  add_row(          sig = "Signature.29", aet = "tobacco") %>%
  add_row(          sig = "Signature.30", aet = "unknown")

aetiologies_list <- aetiology_mapping$aet %>% as.list() %>% unique()

aetiologies <- sig_data
for (a in aetiologies_list) {
  aetiologies[a] = 0
  print(a)
  list <- ( aetiology_mapping %>% filter(aet == a) %>% select(sig) %>% as.list() )
  print(list)
  if (length(list$sig) > 1) {
    aetiologies[[a]] <- rowSums(aetiologies[,list$sig])
    #print(aetiologies %>% mutate(aetiologies[[a]] = rowSums(aetiologies[,list$sig])))
  }
  if (length(list$sig) == 1) {
    aetiologies[[a]] <- aetiologies[,list$sig]
    #aetiologies <- aetiologies %>% mutate(aetiologies[[a]] = aetiologies[,list$sig])
  }
}



melted_aet <- melt.data.frame(aetiologies %>% dplyr::select(c("sample", "normal_db", "tumor_db", "normal_original", "tumor_original", "normal_stela", "tumor_stela", "diff", contains(unlist(aetiologies_list)))), 
                          id.vars = c("sample", "normal_db", "tumor_db", "normal_original", "tumor_original", "normal_stela", "tumor_stela", "diff"))

#aetiologies
#melted_aet

melted_aet %>%
  filter(variable != "unknown") %>%
  filter(value > 0) %>%
  arrange((tumor_stela)) %>% 
  mutate(sample=factor(sample, levels=unique(sample))) %>% 
  ggplot(aes(x=sample, y=value, fill=variable)) +
  geom_bar(stat = 'identity') +
  xlab("Sample (low to high tumour telomere length)") + 
  ylab("Known signature aetiology") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#cor.test(aetiologies$tumor_stela, aetiologies$`homologous recombination`)
#cor.test(aetiologies$`homologous recombination`, aetiologies$tumor_stela)


#Function for creating list of colours in plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


