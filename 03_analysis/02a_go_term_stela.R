library(tidyverse)
library(reshape2)

sample_pairs <- read_csv("/path/to/project/sample_pairs.csv")
sample_pairs <- sample_pairs %>% add_column(sample = paste(sample_pairs$normal_db, sample_pairs$tumor_db, sep = "_"))
sample_pairs <- sample_pairs %>% mutate(diff = (normal_stela - tumor_stela))
sample_pairs

sample_GO_terms <- read_csv("/path/to/project/GORILLA-out/sample_GO_terms_10.csv")
sample_GO_terms

hist(sample_pairs$tumor_stela)

sample_pairs$bin <- as.integer(cut(sample_pairs$tumor_stela, breaks = seq(2, 8, by = 1), labels = 2:7))+1
#sample_pairs

go_bin_table <- tibble(GO_term=sample_GO_terms$GO_term)
for (i in seq(2,7,by=1)){
  go_bin_table[[as.character(i)]] <- integer(length(go_bin_table$GO_term))
}

go_bin_table

for (i in colnames(sample_GO_terms)[-1]) {
  bin <- sample_pairs %>% filter(sample == i) %>% select(bin) %>% pull()
  go_bin_table[as.character(bin)] <- go_bin_table[as.character(bin)] + sample_GO_terms[i]
}

go_filter_thresh <- 0

go_bin_table <- go_bin_table[rowSums(go_bin_table[,-1]) > go_filter_thresh,]

go_bin_table

hist(lapply(go_bin_table, is.numeric))

melted_go_bin <- melt(go_bin_table, id="GO_term")

melted_go_bin %>%
  #mutate(sample=factor(sample, levels=unique(sample))) %>%
  ggplot(aes(x=variable, y=value, fill=GO_term)) +
  geom_bar(stat = 'identity') +
  xlab("Tumour telomere length (floor)") +
  ylab("GO terms count")

#as.integer(sample_pairs$bin)+1
