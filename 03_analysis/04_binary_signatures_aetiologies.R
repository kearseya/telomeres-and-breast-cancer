detach("package:dplyr", character.only = TRUE)
library("dplyr", character.only = TRUE)

aetiologies_list
bin_sig_aet <- all %>% select("sample", signature_names_list, contains(as.character(aetiologies_list)))
bin_sig_aet <- bin_sig_aet[, !apply(bin_sig_aet == 0, 2, all)]

convert_binary <- function(x){
  case_when(x < 0.05 ~ 0, x >= 0.05 ~ 1)
}

bin_sig_aet <- mutate_if(bin_sig_aet, is.numeric, convert_binary)
bin_sig_aet <- bin_sig_aet %>% mutate(sumsig = rowSums(select(., contains(signature_names_list))))

bin_sig_aet <- mutate_if(bin_sig_aet, is.numeric, as.factor)

bin_sig_aet <- left_join(select(all, c("sample", "tumor_stela", "chromothripsis_p0.04_n8")), bin_sig_aet)
bin_sig_aet$sumsig <- as.numeric(bin_sig_aet$sumsig)


bin_sig_aet %>%
  arrange((tumor_stela)) %>%
  select(tumor_stela, sumsig)

cor.test(bin_sig_aet$tumor_stela, bin_sig_aet$sumsig)


simple_bisereal <- function(con_var, bin_var){
  biserial.cor(bin_sig_aet[[con_var]], bin_sig_aet[[bin_var]])
}

bi_cor_sig_aet_table <- tibble(variable = character(), biserial = double())
for (b in colnames(bin_sig_aet)[-c(1,2)]){
  bi_cor_sig_aet_table <- add_row(bi_cor_sig_aet_table, variable = b, biserial = simple_bisereal("tumor_stela", b))
}


colnames(bin_sig_aet)[-c(1:3)] <- paste(colnames(bin_sig_aet)[-c(1:3)], "bin", sep = ".")
bin_sig_aet
all <- left_join(all, bin_sig_aet)


