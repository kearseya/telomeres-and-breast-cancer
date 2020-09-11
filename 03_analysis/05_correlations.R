library(tidyverse) 
library(cluster)
library(factoextra)
library(dendextend)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ggthemes)

unloadNamespace("tidyverse")
library(tidyverse)
detach("package:dplyr", character.only = TRUE)
library("dplyr", character.only = TRUE)

signature_number_list <- seq(1:30)[-c(4,9,11,12,14,15,16,19,21,22,23,24,25,27,28,29)]

make_sig_cor_table <- function() {
  sig_cor_table <- tibble(signature=character(), pval=numeric(), statistic=numeric(), confint_low=numeric(), confint_high=numeric(), estimate=numeric(), nullval=numeric())
  for (i in signature_number_list) {
    sig <- paste("Signature", i, sep=".")
    #print(sig)
    sig_cor <- cor.test(all$tumor_stela, all[[sig]])
    sig_cor_table <- add_row(sig_cor_table, signature=sig, pval=sig_cor$p.value, statistic=sig_cor$statistic, confint_low=sig_cor$conf.int[1], confint_high=sig_cor$conf.int[2], estimate=sig_cor$estimate, nullval=sig_cor$null.value)
  }
  return(sig_cor_table)
}

sig_cor_table <- make_sig_cor_table()

print(sig_cor_table, n=30)

make_aet_cor_table <- function() {
  aet_cor_table <- tibble(aetiology=character(), pval=numeric(), statistic=numeric(), confint_low=numeric(), confint_high=numeric(), estimate=numeric(), nullval=numeric())
  for (i in c("deamination of 5mC", "APOBEC", "homologous recombination")) {
    aet <- paste(i, "", sep = "")
    aet_cor <- cor.test(all$tumor_stela, all[[aet]])
    aet_cor_table <- add_row(aet_cor_table, aetiology=aet, pval=aet_cor$p.value, statistic=aet_cor$statistic, confint_low=aet_cor$conf.int[1], confint_high=aet_cor$conf.int[2], estimate=aet_cor$estimate, nullval=aet_cor$null.value)
  }
  return(aet_cor_table)
}

aet_cor_table <- make_aet_cor_table()


make_sig_cor_table_npi <- function() {
  sig_cor_table <- tibble(signature=character(), pval=numeric(), statistic=numeric(), confint_low=numeric(), confint_high=numeric(), estimate=numeric(), nullval=numeric())
  for (i in signature_number_list) {
    sig <- paste("Signature", i, sep=".")
    #print(sig)
    sig_cor <- cor.test(all$`NPI SCORE`, all[[sig]])
    sig_cor_table <- add_row(sig_cor_table, signature=sig, pval=sig_cor$p.value, statistic=sig_cor$statistic, confint_low=sig_cor$conf.int[1], confint_high=sig_cor$conf.int[2], estimate=sig_cor$estimate, nullval=sig_cor$null.value)
  }
  return(sig_cor_table)
}

sig_cor_table_npi <- make_sig_cor_table_npi()

print(sig_cor_table_npi, n=30)

make_aet_cor_table_npi <- function() {
  aet_cor_table <- tibble(aetiology=character(), pval=numeric(), statistic=numeric(), confint_low=numeric(), confint_high=numeric(), estimate=numeric(), nullval=numeric())
  for (i in c("deamination of 5mC", "APOBEC", "homologous recombination")) {
    aet <- paste(i, "", sep = "")
    aet_cor <- cor.test(all$`NPI SCORE`, all[[aet]])
    aet_cor_table <- add_row(aet_cor_table, aetiology=aet, pval=aet_cor$p.value, statistic=aet_cor$statistic, confint_low=aet_cor$conf.int[1], confint_high=aet_cor$conf.int[2], estimate=aet_cor$estimate, nullval=aet_cor$null.value)
  }
  return(aet_cor_table)
}

aet_cor_table_npi <- make_aet_cor_table_npi()

print(aet_cor_table_npi, n=30)

n_sigx <- tibble(signature=character(), n = numeric())
for (s in signature_names_list) {
  n_sigx <- add_row(n_sigx, signature = s, n = length(all %>% filter(all[[s]] > 0) %>% select(s) %>% pull()))
  #print(cor.test(all$tumor_stela, all[[s]]) )
}
n_sigx

n_aetx <- tibble(aetiology=character(), n = numeric())
for (a in aetiologies_list) {
  print(a)
  n_aetx <- add_row(n_aetx, aetiology = a, n = length(all %>% filter(all[[a]] > 0) %>% select(a) %>% pull()))
  #print(cor.test(all$tumor_stela, all[[s]]) )
}
n_aetx


cor.test(all$`NPI SCORE`, all$sumsig.bin)


signature_names_list <- c()
pos <- 1
for (i in signature_number_list) {
  signature_names_list[[pos]] <- paste("Signature", i, sep = ".")
  pos <- pos+1
}

telomere_length_name_list <- c("normal_stela", "tumor_stela", "diff")

chromothripsis_cols <- c("chromothripsis_p0.04_n8" ,           "chromothripsis_p0.04_n12",           
                         "chromothripsis_p0.04_n20",           "chromothripsis_p0.01_n10",           
                         "chromothripsis_p0.01_n12",           "chromothripsis_from_cn_profile_only")

menopaus_cols <- c("SEX", "Menopausal", "Menopause age")

cancer_cols <- c("TISSUE", "DIAGNOSIS", "GRADE")

status_cols <- c("ER STATUS", "PGR STATUS", "HER2 STATUS")

therapy_cols <- c("NEO-ADJUVANT CHEMOTHERPY", "ADJUVANT CHEMOTHERAPY", "ADJUVANT HORMONE THERAPY")

prognosis_cols <- c("Time to Death (days)"               , "CANCER DEATH RELATED"  ,             
                    "Diagnosis to follow up/death (days)", "Operation to follow up (days)"   ,   
                    "Operation to death (days)"          , "Operation to follow/death (days)"   ,
                    "Censoring Status", "AGE AT OPERATION", "NPI score")

other_cols <- c("Number of bands",                    
                "Mean"        ,                    "St. Deviation",                      
                "Below_3.81"    ,                   "triple_negative"  )



add_gene_status <- function(gene_list = c()) {
  for (g in gene_list) {
    if (g %in% all_genes) {
      tmp_gene_pres <- sample_genes %>% filter(gene == g) %>% t() %>% as.data.frame()
      tmp_gene_pres$sample <- rownames(tmp_gene_pres)
      colnames(tmp_gene_pres) <- c(g, "sample")
      print(tmp_gene_pres)
      tmp_gene_pres <- tmp_gene_pres %>% filter(sample %in% all$sample)
      print(tmp_gene_pres)
      tmp_gene_pres[[g]] <- factor(tmp_gene_pres[[g]])
      #tmp_gene_pres[[g]] %>% lapply(as.factor)
      
      #tmp_gene_pres
      assign("all", left_join(all, tmp_gene_pres), envir = globalenv())
    }
  }
}

opt_k_value_no_bin <- function(var_list=c()) {
  df <- all %>% select("sample", var_list)
  df <- df[, !apply(df == 0, 2, all)]
  df[,-1] <- scale(df[,-1])
  fviz_nbclust(df, FUN = hcut, method = "wss")
}

opt_k_value_with_bin <- function(var_list=c(), bin_vars=c()) {
  df <- all %>% select("sample", var_list)
  df <- df[, !apply(df == 0, 2, all)]
  bin_levels <- factor(c(0,1), levels = 0:1)
  df[bin_vars] <- lapply(df[bin_vars], function(x) factor(x, bin_levels))
  dist_obj <- daisy(df[,-1], metric = "gower", stand = TRUE)
  #not working
  #nbtest <- NbClust(diss = as.matrix(dist_obj), distance = NULL, method = "ward.D", index = "silhouette")
  #fviz_nbclust(nbtest)
  #semi-working
  sil_plot <- (silhouette(all$cl_members %>% sapply(as.integer), dist_obj))
  #sil_plot
  return(sil_plot)
}

hclust_apply <- function(k_val, cl_metric, var_list=c(), bin_vars=c(), add_cl_all=FALSE) {
  df <- all %>% dplyr::select("tumor_db", var_list) 
  df <- df[, !apply(df == 0, 2, all)]
  
  ## for no binary variables
  # df <- mutate_if(df, is.numeric, scale)
  # df[,-1] <- scale(df[,-1])
  # dist_obj <- dist(df)
  # dist_obj <- as.matrix(dist_obj, labels=TRUE)
  # colnames(dist_obj) <- rownames(dist_obj) <- df[['sample']]
  
  ## for when binary variables
  if (cl_metric == "gower") {
    bin_levels <- factor(c(0,1), levels = 0:1)
    df[bin_vars] <- lapply(df[bin_vars], function(x) factor(x, bin_levels))
  }
  
  dist_obj <- daisy(df[,-1], metric = cl_metric, stand = TRUE)
  
  tree <- hclust(dist_obj, method = "ward.D")
  tree$labels <- df$tumor_db
  cl_members <- cutree(tree = tree, k = k_val)

  clust_mem <- as.data.frame(cl_members)
  clust_mem$tumor_db <- row.names(clust_mem)
  clust_all <- left_join(dplyr::select(all, -dplyr::contains("cl_members")), clust_mem)
  clust_all$cl_members <- factor(clust_all$cl_members)

  assign("all", clust_all, envir = globalenv())
  
  tree_cols <- c("#F8766D", "#619CFF", "#00BA38")
  #dend <- set(as.dendrogram(rev(tree)), "branches_lwd", 3)
  temp_col <- tree_cols[as.numeric(all$cl_members)]
  temp_col <- temp_col[order.dendrogram(as.dendrogram(rev(tree)))]
  temp_col <- factor(temp_col, unique(temp_col))
  rev(tree) %>%
    as.dendrogram() %>%
    color_branches(clusters = as.numeric(temp_col), col = levels(temp_col)) %>%
    set("branches_lwd", 2.25) %>%
    plot()
  
  rectangles <- rect.hclust(tree = rev(tree), k = k_val, which = 1:k_val, border = "transparent", cluster = cl_members)
  first_in_clust <- head(cumsum(c(1, lengths(rectangles))), -1)
  mid_y_clust <- weighted.mean(rev(tree$height)[2:3], c(4, 1))
  text(x=first_in_clust, y=mid_y_clust, col=c("#00BA38", "#F8766D", "#619CFF"), labels=c("1b", "1a", "2"), font=2, cex = 2)
  
  # par(mfrow = c(1,1))
  # tree_cols <- gg_color_hue(k_val)
  # plot(x = rev(tree), labels = row.names(tree), cex = 1, xlab = "", sub = "", main = "")
  # 
  # order_appear <- c()
  # or_pos <- 1
  # for (s in rev(tree$labels[tree$order])) {
  #   tmp_cl <- all %>% filter(tumor_db == s) %>% select(cl_members) %>% pull() %>% as.numeric()
  #   if (!tmp_cl %in% order_appear) {
  #     order_appear[[or_pos]] <- tmp_cl
  #     or_pos <- or_pos + 1
  #   }
  # }
  # # rectangles <- rect.hclust(tree = rev(tree), k = k_val, which = 1:k_val, border = tree_cols[order_appear], cluster = cl_members)
  # # first_in_clust <- head(cumsum(c(1, lengths(rectangles))), -1)
  # # mid_y_clust <- weighted.mean(rev(tree$height)[2:3], c(4, 1))
  # # text(x=first_in_clust, y=mid_y_clust, col=tree_cols[order_appear], labels=order_appear, font=2)
  # # #fviz_cluster(list(data = dist_obj, cluster = cl_members))
  # 
  # rectangles_2 <- rect.hclust(tree = rev(tree), k = 2, which = 1:2, border = tree_cols[c(1,2)], cluster = cutree(tree = tree, k = 2))
  # rectangles_3 <- rect.hclust(tree = rev(tree), h=1, which = 1:2, border = tree_cols[order_appear], cluster = cl_members)
  # rectangles <- rect.hclust(tree = rev(tree), k = k_val, which = 1:k_val, border = "transparent", cluster = cl_members)
  # first_in_clust <- head(cumsum(c(1, lengths(rectangles))), -1)
  # mid_y_clust <- weighted.mean(rev(tree$height)[2:3], c(4, 1))
  # text(x=first_in_clust, y=mid_y_clust, col=tree_cols[order_appear], labels=c("1b", "1a", "2"), font=2)
  # #fviz_cluster(list(data = dist_obj, cluster = cl_members))
 
  return(tree)
}

plot_tree_sigs <- function(k_val, tree) {
  cl_members <- cutree(tree = tree, k = k_val)
  #sig_clust_cols <- rainbow(length(signature_names_list))
  sig_clust_cols <- mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(signature_names_list))
  clust_cols <- gg_color_hue(k_val)[c(1,3,2)]
  names(sig_clust_cols) <- signature_names_list
  clust_sig_graph_list <- c()
  
  pos <- 1
  for (i in unique(all$cl_members)) {
    tmp_dat <- all %>% filter(cl_members == i)
    tmp_melted <- melt.data.frame(tmp_dat %>% select("tumor_db", signature_names_list, "tumor_stela"), id.vars = c("tumor_db", "tumor_stela")) %>% filter(value > 0)
    tmp_plot <-tmp_melted %>%
      arrange((tumor_stela)) %>%
      mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
      ggplot(aes(x=tumor_db, y=value, fill=variable)) +
      geom_bar(stat = 'identity') +
      coord_flip() +
      xlab("") +
      ylab(paste("Signature percentage cluster", i, sep = " ")) +
      scale_fill_manual(name = "Signatures", values = sig_clust_cols) +
      #scale_fill_brewer(name = "Signatures", palette = "Set1") +
      scale_y_continuous(limits = c(0,1.01), expand = c(0,0)) +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(color = clust_cols[pos], size = 15, face = "bold"),
            axis.title.y = element_text(color = clust_cols[pos], size = 10),
            panel.grid.major = element_blank(),  #element_line(size = 0.4, color = "grey"),
            panel.grid.minor = element_blank(),  #element_line(size = 0.2, color = "grey"),
            panel.border = element_blank(),
            legend.position = "none")
    
    #assign(paste("cluster_sig_plot", i, sep="_"), tmp_plot)
    clust_sig_graph_list[[pos]] <- tmp_plot
    pos <- pos+1
  }
  
  heights <- c()
  pos <- 1
  total_rows <- nrow(all)
  for (i in unique(all$cl_members)){
    heights[pos] <- (all %>% filter(cl_members == i) %>% nrow() / total_rows)
    pos <- pos+1
  }
  
  # all_plot <- melt.data.frame(all %>% select("tumor_db", signature_names_list, "tumor_stela"), id.vars = c("tumor_db", "tumor_stela")) %>% filter(value > 0) %>%
  #   ggplot(aes(x=tumor_db, y=value, fill=variable)) +
  #   geom_bar(stat = 'identity') +
  #   scale_fill_manual(name = "Signatures", values = sig_clust_cols) +
  #   theme(legend.title = element_text(size = 50),
  #         legend.text = element_text(size = 30),
  #         legend.background = element_blank())

  #tmp <- ggplot_gtable(ggplot_build(all_plot))
  #leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #shared_legend <- tmp$grobs[[leg]]
  #shared_legend <- get_legend(all_plot) # when using add heights = c(0, heights)
  #grid.arrange(grobs = clust_sig_graph_list, heights = heights)
  ggarrange(plotlist = clust_sig_graph_list[c(1,3,2)], heights = heights[c(1,3,2)], ncol = 1, legend = "right", common.legend = TRUE)
}

plot_tree_aetiologies <- function(k_val, tree) {
  cl_members <- cutree(tree = tree, k = k_val)
  present_aet <- all %>% select(!c(contains(".bin"))) %>% select(contains(as.character(aetiologies_list)))
  present_aet <- colnames(present_aet[, !apply(present_aet == 0, 2, all)])
  #aet_clust_cols <- rainbow(length(present_aet)*2)
  clust_cols <- gg_color_hue(k_val)[c(1,3,2)]
  #names(aet_clust_cols) <- present_aet
  clust_aet_graph_list <- c()
  pos <- 1
  for (i in unique(all$cl_members)) {
    #tmp_dat <- all %>% filter(sample %in% names(cl_members[cl_members == i]))
    tmp_dat <- all %>% filter(cl_members == i)
    tmp_melted <- melt.data.frame(tmp_dat %>% select("tumor_db", present_aet, "tumor_stela"), id.vars = c("tumor_db", "tumor_stela")) %>% filter(value > 0)
    tmp_plot <-tmp_melted %>%
      arrange((tumor_stela)) %>%
      mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
      ggplot(aes(x=tumor_db, y=value, fill=variable)) +
      geom_bar(stat = 'identity') +
      coord_flip() +
      xlab("") +
      ylab(paste("Aetiology percentage cluster", i, sep = " ")) +
      #scale_fill_manual(name = "Aetiologies", values = aet_clust_cols) +
      scale_fill_brewer(name = "Aetiologies", palette = "Set1") +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(color = clust_cols[pos], size = 15, face = "bold"),
            axis.title.y = element_text(color = clust_cols[pos], size = 10),
            panel.grid.major = element_blank(),  #element_line(size = 0.4, color = "grey"),
            panel.grid.minor = element_blank(),  #element_line(size = 0.4, color = "grey"),)
            legend.position = "none")
            
    #assign(paste("cluster_sig_plot", i, sep="_"), tmp_plot)
    clust_aet_graph_list[[pos]] <- tmp_plot
    pos <- pos+1
  }
  heights <- c()
  total_rows <- nrow(all)
  pos <- 1
  for (i in unique(all$cl_members)){
    heights[pos] <- (all %>% filter(cl_members == i) %>% nrow() / total_rows)
    pos <- pos+1
  }
  #grid.arrange(grobs = clust_aet_graph_list, heights = heights)
  ggarrange(plotlist = clust_aet_graph_list[c(1,3,2)], heights = heights[c(1,3,2)], ncol = 1, legend = "right", common.legend = TRUE)
}

iter_opt_k <- function(min, max, cl_metric) {
  sil_plot_list <- c()
  n_cols_par <- ceiling((max-min)/2)
  for (k in seq(min, max)) {
    hclust_apply(k, cl_metric, variables_for_hclust, binary_variables)
    if (cl_metric == "gower") {
      sil_plot_list[[k-min+1]] <- opt_k_value_with_bin(variables_for_hclust, binary_variables)
    }
    if (cl_metric == "euclidean") {
      sil_plot_list[[k-min+1]] <- opt_k_value_no_bin(variables_for_hclust)
    }
  }
  par(mfrow=c(2, n_cols_par))
  num_k <- min
  sil_plots <- c()
  for (p in sil_plot_list) {
    hues <- seq(15, 375, length = num_k + 1)
    sil_col <- hcl(h = hues, l = 65, c = 100)[1:num_k]
    sil_plots[[num_k-min+1]] <- fviz_silhouette(p, fill = sil_col) +
      #coord_flip() +
      scale_y_continuous(expand = c(0,0)) +
      labs(title = paste("Average silhouette width:", round(summary(p)$avg.width, 2), sep = " " )) +
      theme_clean() +
      theme(axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(size = 14))
    #plot(p, main = "", col=gg_color_hue(num_k)) #, main = paste("Silhouette plot for k=", (p-1+min), sep = ""))
    num_k <- num_k + 1
  }
  grid.arrange(grobs = sil_plots)
}



gene_status_list <- c("TP53", "BRCA1", "BRCA2", "PTEN", "APOBEC1", "APOBEC3A", "APOBEC3B", "REV1")

add_gene_status(gene_status_list)

bin_sigs <- intersect(colnames(all), paste(signature_names_list, "bin", sep = "."))
bin_aets <- intersect(colnames(all), paste(as.character(aetiologies_list), "bin", sep = "."))



variables_for_hclust <- c("tumor_stela", "chromothripsis_p0.04_n8", "diff", "num_total_mut")
#gene_status_list, "APOBEC", "homologous recombination", "deamination  of 5mC", "AGE AT OPERATION", 
binary_variables <- c("chromothripsis_p0.04_n8")
# gene_status_list, bin_sigs, bin_aets, "chromothripsis_p0.04_n8"
cl_metric <- "gower" #"euclidean"

iter_opt_k(2, 5, cl_metric)

#opt_k_value_with_bin(variables_for_hclust, binary_variables)
#opt_clust

k <- 3

test <- hclust_apply(k, cl_metric, variables_for_hclust, binary_variables, add_cl_all = TRUE)


if (k == 3) {
  all <- all %>% mutate(cl_members = case_when(all$cl_members == 1 ~ "1a",
                                               all$cl_members == 2  ~ "2",
                                               all$cl_members == 3 ~ "1b"))
}


plot_tree_sigs(k, test)
plot_tree_aetiologies(k, test)

#three_boxplot_show()
#three_histogram_show()


make_grid_plot(c("AGE AT OPERATION", "tumor_stela", "NPI SCORE", "diff", "num_total_mut", "chromothripsis_p0.04_n8"),
               c("Age at operation (years)", "Tumour TL (kb)", "NPI Score", "Normal-Tumour TL (kb)", "Number of total mutations", "Chromothripsis status"),
               c("b", "b", "b", "b", "b", "br"), 
               c(NA,NA,NA,NA,NA,NA))

make_gene_hists(c("TP53", "BRCA1", "BRCA2", "APOBEC1", "APOBEC3A", "APOBEC3B", "REV1", "PTEN"))

# make_histogram("tumor_stela", 0.2)
make_boxplot("tumor_stela", "Tumour TL")
# #show_all_sig_box()
make_barplot("GRADE", "Grade")

make_histogram("AGE AT OPERATION", "Age at operation", binw = 5)
# > all %>% filter(cl_members == 1) %>% select(tumor_stela) %>% max()
# [1] 3.735082
# > all %>% filter(cl_members == 3) %>% select(tumor_stela) %>% min()
# [1] 4.640192
# > all %>% filter(cl_members == 2) %>% select(tumor_stela) %>% min()
# [1] 3.74


#Testing difference between variables in clusters

test_diff <- function(group_a, group_b, variable) {
  #t.test(all %>% filter(cl_members == group_a) %>% select(variable), all %>% filter(cl_members == group_b) %>% select(variable))
  wilcox.test(all %>% filter(cl_members == group_a) %>% select(variable) %>% sapply(as.numeric), all %>% filter(cl_members == group_b) %>% select(variable) %>% sapply(as.numeric))
}

test_diff("1a", "2", "tumor_stela")
test_diff("1a", "2", "diff")
test_diff("1a", "2", "AGE AT OPERATION")
test_diff("1a", "2", "NPI SCORE")
test_diff("1a", "2", "PTEN")

test_diff("1a", "2", "Signature.13") # signatures: 2, (and 8?) shows some difference, 5 shows signficant difference between chromothripsis seperated samples, due to other signatures slightly increasing making it less represented
test_diff("1b", "2", "APOBEC")

all_bi("Signature.5", "chromothripsis_p0.04_n8")
#bi_cor_sig_aet_table (var + tumor_stela)
all_bi("tumor_stela", "PTEN")


biserial.cor(all %>% select(num_total_mut, chromothripsis_p0.04_n8) %>% drop_na(num_total_mut, chromothripsis_p0.04_n8) %>% filter(num_total_mut < 10000)%>%select(num_total_mut)%>%pull(), all %>% select(num_total_mut, chromothripsis_p0.04_n8) %>% drop_na(num_total_mut, chromothripsis_p0.04_n8) %>%filter(num_total_mut < 10000) %>% select(chromothripsis_p0.04_n8) %>% pull(), use = c("all.obs", "complete.obs"), level = 1)



t.test(all %>% filter(cl_members == "1a") %>% select("tumor_stela") , all %>% filter(cl_members == "2") %>% select("tumor_stela") %>% sapply(as.numeric))


cor.test(all$tumor_stela, all$Signature.13) ## 13 is close, APOBEC
cor.test(all$tumor_stela, all$APOBEC)

cor.test(all$Signature.5, all$chromothripsis_p0.04_n8) # 8

# # not enough samples
# t.test(all %>% filter(cl_members == "1a") %>% select(num_total_mut) %>% pull(), 
#        all %>% filter(cl_members == "2") %>% select(num_total_mut) %>% pull())
# 
# t.test(all %>% filter(tumor_stela < 3.81) %>% select(num_total_mut), 
#        all %>% filter(tumor_stela >= 3.81) %>% select(num_total_mut))
# 
# t.test(all %>% filter(cl_members == "1a") %>% filter(num_total_mut < 10000) %>% select(num_total_mut), 
#        all %>% filter(cl_members == "2") %>% filter(num_total_mut < 10000) %>% select(num_total_mut))
# 
# t.test(all %>% filter(tumor_stela < 3.81) %>% filter(num_total_mut < 10000) %>% select(num_total_mut), 
#        all %>% filter(tumor_stela >= 3.81) %>% filter(num_total_mut < 10000) %>% select(num_total_mut))

#spearman or pearsons

custom_theme_1 <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                        axis.text.y = element_text(size = 15),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 15),
                        
                        plot.title = element_text(face = "bold",
                                                  size = rel(1.2), hjust = 0.5),
                        text = element_text(),
                        panel.background = element_rect(colour = NA),
                        plot.background = element_rect(colour = NA),
                        panel.border = element_rect(colour = NA),
                        axis.title = element_text(face = "bold",size = rel(1)),
                        
                        axis.text = element_text(), 
                        axis.line = element_line(colour="black"),
                        axis.ticks = element_line(),
                        panel.grid.major = element_line(colour="#f0f0f0"),
                        panel.grid.minor = element_line(colour="#f0f0f0"), #element_blank(),
                        legend.key = element_rect(colour = NA),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.key.size= unit(0.2, "cm"),
                        #legend.margin = unit(0, "cm"),
                        #legend.title = element_text(face="italic"),
                        plot.margin=unit(c(10,5,5,5),"mm"),
                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                        strip.text = element_text(face="bold"),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.minor.y = element_blank())





make_histogram <- function(variable, title, binw) {
  ggplot(all, aes(all[[variable]], fill = cl_members)) +
    geom_histogram(binwidth = binw) +
    labs(title=paste(variable, "by cluster"), sep = " ") +
    xlab("") +
    labs(title=title) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_line(size = 0.1, color = "grey"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(.~cl_members)
}

make_boxplot <- function(variable, title) {
  if (variable != "num_total_mut") {
    return(ggplot(all, aes(all[[variable]], fill = cl_members)) +
      geom_boxplot() +
      coord_flip() +
      #labs(title=paste(variable, "by cluster"), sep = " ") +
      xlab("") +
      labs(title=title) +
      #scale_y_continuous(expand = c(0,0)) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            panel.background = element_blank(), 
            panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
            panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      facet_grid(.~cl_members))
  }
  if (variable == "num_total_mut") {
      return(ggplot( all %>% filter(num_total_mut < 30000), aes(num_total_mut, fill = cl_members)) +
        geom_boxplot() +
        coord_flip() +
        #labs(title=paste(variable, "by cluster"), sep = " ") +
        xlab("") +
        labs(title=title) +
        #scale_y_continuous(expand = c(0,0)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              panel.background = element_blank(),
              panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
              panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        facet_grid(.~cl_members))
  }
}

make_bin_barplot <- function(variable, title) {
  ggplot(all, aes(all[[variable]], fill = cl_members)) +
    geom_bar() +
    scale_x_continuous(breaks = c(0,1), labels = c("-ve", "+ve")) +
    xlab("") +
    labs(title=title) +
    theme(axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(.~cl_members)
}

make_barplot <- function(variable, title) {
  ggplot(all, aes(all[[variable]], fill = cl_members)) +
    geom_bar() +
    xlab("") +
    labs(title=title) +
    theme(axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    facet_grid(.~cl_members)
  
}

make_grid_plot <- function(var_list, title_list, type_list, binw_list) {
  grid_graphs <- c()
  for (i in seq(1:length(var_list))) {
    if (type_list[i] == "b") {
      grid_graphs[[i]] <- make_boxplot(var_list[i], title_list[i])
    }
    if (type_list[i] == "h") {
      grid_graphs[[i]] <- make_histogram(var_list[i], title_list[i], binw_list[i])
    }
    if (type_list[i] == "r") {
      grid_graphs[[i]] <- make_barplot(var_list[i], title_list[i])
    }
    if (type_list[i] == "br") {
      grid_graphs[[i]] <- make_bin_barplot(var_list[i], title_list[i])
    }
  }
  grid.arrange(grobs = grid_graphs)#, layout_matrix = rbind(c(1,2), c(3,4), c(5)))
}


make_gene_hists <- function(gene_list) {
  grid_graphs <- c()
  pos <- 0
  for (g in gene_list) {
    pos <- pos + 1
    grid_graphs[[pos]] <- ggplot(all, aes_string(g, fill = "cl_members")) +
      geom_bar() +
      scale_x_discrete(breaks = c(0,1), labels = c("no", "yes")) +
      labs(title=paste(g, "mutation presence by cluster", sep = " ")) +
      xlab(g) +
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
            panel.grid.minor = element_blank(),
            legend.position = "none") +
      facet_grid(.~cl_members)
    
  }
  grid.arrange(grobs = grid_graphs)
}






#make_grid_plot(c("AGE AT OPERATION", "tumor_stela", "chromothripsis_p0.04_n8"), c("b", "b", "h"), c(0,0,1))

## HISTOGRAMS

ah <- ggplot(all, aes(all$"AGE AT OPERATION", fill = cl_members)) +
  geom_histogram(binwidth = 5) +
  labs(title="Age at operation by cluster") +
  theme(legend.position = "none") +
  facet_grid(.~cl_members)
ah

th <- ggplot(all, aes(tumor_stela, fill = cl_members)) +
  geom_histogram(binwidth = 0.5) +
  labs(title="Tumour telomere length by cluster") +
  theme(legend.position = "none") +
  facet_grid(.~cl_members)
th

ch <- ggplot(all, aes(chromothripsis_p0.04_n8, fill = cl_members)) +
  geom_bar() +
  scale_x_continuous(breaks = c(0,1), labels = c("no", "yes")) +
  labs(title="Chromothrips by cluster") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
ch

grid.arrange(th, ah, ch)

gth <- ggplot(all, aes(TP53, fill = cl_members)) +
  geom_bar() +
  scale_x_discrete(breaks = c(0,1), labels = c("no", "yes")) +
  labs(title="PT53 mutation presence by cluster") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
gth

b1h<- ggplot(all, aes(BRCA1, fill = cl_members)) +
  geom_bar() +
  scale_x_discrete(breaks = c(0,1), labels = c("no", "yes")) +
  labs(title="BRCA1 mutation presence by cluster") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
b1h

b2h <- ggplot(all, aes(BRCA2, fill = cl_members)) +
  geom_bar() +
  scale_x_discrete(breaks = c(0,1), labels = c("no", "yes")) +
  labs(title="BRCA2 mutation presence by cluster") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
b2h

gph <- ggplot(all, aes(PTEN, fill = cl_members)) +
  geom_bar() +
  scale_x_discrete(breaks = c(0,1), labels = c("no", "yes")) +
  labs(title="PTEN mutation presence by cluster") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
gph

grid.arrange(gth, b1h, b2h, gph)

## BOXPLOTS

ab <- ggplot(all, aes(all$"AGE AT OPERATION", fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="Age at operation by cluster") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_line(size = 0.1, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
ab

tb <- ggplot(all, aes(tumor_stela, fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="Tumour telomere length by cluster") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_line(size = 0.1, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
tb

nb <- ggplot(all, aes(all$"NPI SCORE", fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="NPI Score by cluster") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_line(size = 0.1, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
nb

db <- ggplot(all, aes(all$"diff", fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="Telomere difference by cluster") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_line(size = 0.1, color = "grey"),
        legend.position = "none") +
  facet_grid(.~cl_members)
db

grid.arrange(tb, db, ab, nb, ch)





signature_boxplot <- function(sig_num) {
  ggplot(all, aes(all[[paste("Signature", sig_num, sep = ".")]], fill = cl_members)) +
    geom_boxplot() +
    coord_flip() +
    labs(title=paste("Signautre", sig_num, "by cluster", sep = " ")) +
    xlab("perc") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_line(size = 0.1, color = "grey"),
          legend.position = "none") +
    facet_grid(.~cl_members)
}

show_all_sig_box <- function() {
  sig_boxs <- c()
  pos <- 1
  for (s in signature_number_list) {
    sig_boxs[[pos]] <- signature_boxplot(s)
    pos <- pos + 1
  }
  grid.arrange(grobs = sig_boxs)
}
signature_boxplot(1)

show_all_sig_box()





three_boxplot_show <- function(){
  ab <- ggplot(all, aes(all$"AGE AT OPERATION", fill = cl_members)) +
    geom_boxplot() +
    coord_flip() +
    labs(title="Age at operation by cluster") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_line(size = 0.1, color = "grey"),
          legend.position = "none") +
    facet_grid(.~cl_members)
  ab
  
  tb <- ggplot(all, aes(tumor_stela, fill = cl_members)) +
    geom_boxplot() +
    coord_flip() +
    labs(title="Tumour telomere length by cluster") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_line(size = 0.1, color = "grey"),
          legend.position = "none") +
    facet_grid(.~cl_members)
  tb
  
  nb <- ggplot(all, aes(all$"NPI SCORE", fill = cl_members)) +
    geom_boxplot() +
    coord_flip() +
    labs(title="NPI Score by cluster") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          panel.grid.minor = element_line(size = 0.1, color = "grey"),
          legend.position = "none") +
    facet_grid(.~cl_members)
  nb
  
  grid.arrange(tb, ab, nb)
}

three_histogram_show <- function() {
  ah <- ggplot(all, aes(all$"AGE AT OPERATION", fill = cl_members)) +
    geom_histogram(binwidth = 5) +
    labs(title="Age at operation by cluster") +
    theme(legend.position = "none") +
    facet_grid(.~cl_members)
  ah
  
  th <- ggplot(all, aes(tumor_stela, fill = cl_members)) +
    geom_histogram(binwidth = 0.5) +
    labs(title="Tumour telomere length by cluster") +
    theme(legend.position = "none") +
    facet_grid(.~cl_members)
  th
  
  ch <- ggplot(all, aes(chromothripsis_p0.04_n8, fill = cl_members)) +
    geom_bar() +
    scale_x_continuous(breaks = c(0,1), labels = c("no", "yes")) +
    labs(title="Chromothrips by cluster") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_line(size = 0.4, color = "grey"),
          legend.position = "none") +
    facet_grid(.~cl_members)
  ch
  
  grid.arrange(th, ah, ch)
}




all %>%
  arrange((tumor_stela)) %>%
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=tumor_stela, fill=cl_members)) +
  geom_bar(stat = "identity") +
  xlab("Sample") + 
  ylab("Tumour telomere length (kb)") +
  labs(fill = "Cluster") +
  #scale_fill_discrete(labels = c("Chromothripsis +", "Chromothripsis -")) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept=3.81, linetype='dashed', size=0.2) +
  geom_point(aes(y=normal_stela), alpha=0.5, color = "darkgreen", size = 3) +
  theme_foundation(base_size=14, base_family="helvetica") +
  custom_theme_1


all %>%
  arrange((diff)) %>%
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=diff, fill=cl_members)) +
  geom_bar(stat = "identity") +
  xlab("Sample") + 
  ylab("Normal-Tumour TL (kb)") +
  labs(fill = "Cluster") +
  #scale_fill_discrete(labels = c("Chromothripsis +", "Chromothripsis -")) +
  #scale_y_continuous(expand = c(0,0)) +
  theme_foundation(base_size=14, base_family="helvetica") +
  custom_theme_1



all %>%
  #filter(num_total_mut < 10000) %>%
  arrange((num_total_mut)) %>%
  mutate(tumor_db=factor(tumor_db, levels=unique(tumor_db))) %>% 
  ggplot(aes(x=tumor_db, y=num_total_mut, fill=cl_members)) +
  geom_bar(stat = "identity") +
  xlab("Sample") + 
  ylab("Number of mutations") +
  labs(fill = "Cluster") +
  #scale_fill_discrete(labels = c("Chromothripsis +", "Chromothripsis -")) + # for cl = 2
  scale_y_continuous(expand = c(0,0)) +
  theme_foundation(base_size=14, base_family="helvetica") +
  custom_theme_1


# #Not used as grade part of NPI calculation
# ggplot(all, aes(GRADE, fill = cl_members)) +
#   geom_bar() +
#   labs(title="Grade by cluster") +
#   theme(panel.background = element_blank(), 
#        panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
#        legend.position = "none") +
#   facet_grid(.~cl_members)


