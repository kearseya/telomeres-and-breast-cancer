library(ltm)
library(tidyverse)
detach("package:dplyr", character.only = TRUE)
library("dplyr", character.only = TRUE)



clin <- read_csv("/path/to/project/data/BC_clinical_and_TL_data.csv")
clin

colnames(clin)



sample_pairs <- read_csv("/path/to/project/sample_pairs.csv")
sample_pairs <- sample_pairs %>% add_column(sample = paste(sample_pairs$normal_db, sample_pairs$tumor_db, sep = "_"))
sample_pairs <- sample_pairs %>% mutate(diff = (normal_stela - tumor_stela))
sample_pairs

sample_genes <- read_csv("/path/to/project/data/sample_gene_table_10.csv")
sample_genes %>% filter(gene == "TP53")

hist(sample_pairs$tumor_stela)
clin$'chromothripsis_p0.04_n8'

length(colnames(sample_genes)[c(-1,-2,-length(colnames(sample_genes)))])

make_con_table <- function(variable, short_thresh) {
  con_table <- tibble(tumour_tel_class=c("short", "long"), gene=0, no_gene=0)
  for (i in colnames(sample_genes)[c(-1,-2,-length(colnames(sample_genes)))]) {
    if (sample_pairs %>% filter(sample == i) %>% select(tumor_stela) %>% pull() <= short_thresh) {
      if (sample_genes %>% filter(gene == variable) %>% select(i) %>% pull == 1) {
        con_table[1,2] <- con_table[1,2] + sample_genes %>% filter(gene == variable) %>% select(i) %>% pull()
      }
      else {
        con_table[1,3] <- con_table[1,3] + 1
      }
    }
    if (sample_pairs %>% filter(sample == i) %>% select(tumor_stela) %>% pull() > short_thresh) {
      if (sample_genes %>% filter(gene == variable) %>% select(i) %>% pull == 1) {
        con_table[2,2] <- con_table[2,2] + sample_genes %>% filter(gene == variable) %>% select(i) %>% pull()
      }
      else {
        con_table[2,3] <- con_table[2,3] + 1
      }
    }
  }
  print(con_table)
  return(fisher.test(con_table[c(2,3)]))
}

#make_con_table("TP53", 3.81)
#make_con_table("BRCA2", 3.81)

fisher_table <- tibble(gene=character(), pval=numeric(), confint_low=numeric(), confint_high=numeric(), estimate=numeric(), nullval=numeric())

#small list from papers
drivers <- c("BARD1", "BRCA1", "BRCA2", "CASP8", "CHEK2", "CTLA4",
             "CYP19A1", "FGFR2", "H19", "LSP1", "MAP3K1", "MRE11A",
             "RAD51C", "STK11", "TERT", "TOX3", "XRCC2", "GATA3",
             "PIK3CA", "AKT1", "CDH1", "RB1", "TP53", "PTEN", "hTERT")

all_genes <- sample_genes %>% filter(count >= 40) %>% select(gene) %>% pull()

for (g in drivers) {
  #if (g %in% sample_genes$gene) {
  tmp <- make_con_table(g, 3.81)
  fisher_table <- add_row(fisher_table, gene=g, pval=tmp$p.value, confint_low=tmp$conf.int[1], confint_high=tmp$conf.int[2], estimate=tmp$estimate, nullval=tmp$null.value)
  #}
}

#write_csv(fisher_table, "/path/to/project/data/10_fisher_table.csv")

fisher_table
fisher_table %>% filter(pval < 0.05)
fisher_table %>% filter(pval < 0.05) %>% select(gene) %>% pull()




bi_cor <- function(gene_in) {
  samples_names <- colnames(sample_genes)[c(-1,-2,-length(colnames(sample_genes)))]
  tumour_tel_len <- sample_pairs %>% filter(sample %in% samples_names) %>% dplyr::select(tumor_stela) %>% sapply(as.numeric)
  gene_present <- sample_genes %>% filter(gene == gene_in) %>% dplyr::select(samples_names)
  return(biserial.cor(tumour_tel_len, gene_present, use = c("all.obs", "complete.obs"), level = 1))
}

bi_cor("PTEN")


bi_cor_table <- tibble(gene=character(), bi_cor_val=numeric())

for (g in drivers) {
  if (g %in% sample_genes$gene) {
    bi_cor_table <- add_row(bi_cor_table, gene=g, bi_cor_val=bi_cor(g))
  }
}



clin <- dplyr::rename(clin, normal_db = Normal_DB_id, tumor_db = Tumor_DB_id)

all <- left_join(aetiologies, clin)
all <- left_join(all, gcount_tlen)

all <- all %>% mutate(GRADE = case_when(all$GRADE == "I" ~ 1, 
                                        all$GRADE == "II" ~ 2,
                                        all$GRADE == "III"  ~ 3))

all$GRADE <- as.factor(all$GRADE)


all



biserial.cor(all$tumor_stela %>% sapply(as.numeric), all$Below_3.81)


biserial.cor(all$tumor_stela, all$Below_3.81, use = c("all.obs", "complete.obs"), level = 1)



all_bi <- function(convar, bivar) {
  both <- all %>% dplyr::select(convar, bivar) %>% drop_na(convar, bivar)
  biserial.cor(both[[convar]], both[[bivar]], use = c("all.obs", "complete.obs"), level = 1)
}

all_bi("tumor_stela", "HER2 STATUS")


colnames(all)


cor(all %>% select(tumor_stela, `AGE AT OPERATION`))



trip_cols <- c('ER STATUS', 'PGR STATUS', 'HER2 STATUS')
trip_col_val_tmp <- all %>% select(sample, trip_cols) %>% rename_all(function(x) gsub(" ", "_", x))

trip_col_val_tmp <- trip_col_val_tmp %>%
  mutate(triple_negative = case_when(ER_STATUS == 'Positive' | PGR_STATUS == 'Positive' | HER2_STATUS == 'Positive' ~ as.integer(0),
                                    (ER_STATUS == 'Negitive' & PGR_STATUS == 'Negitive' & HER2_STATUS == 'Negitive') ~ as.integer(1),
                                     ER_STATUS == NA | PGR_STATUS == NA | HER2_STATUS == NA ~ NA_integer_,
                                      TRUE ~ NA_integer_))

all <- left_join(all, trip_col_val_tmp %>% select("sample", "triple_negative"))

# all %>% select(trip_cols)
# tmp$triple_negative


all_bi("tumor_stela", "triple_negative")



