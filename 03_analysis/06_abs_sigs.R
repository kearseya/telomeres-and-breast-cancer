all %>% select(signature_names_list, "num_total_mut")

for (s in signature_names_list) {
  abs_sig_name <- paste(s,"abs", sep=".")
  print(abs_sig_name)
  all[[abs_sig_name]] <- all[[s]]*all$num_total_mut
}

cor.test(all$tumor_stela, all$Signature.5.abs)

all %>%
  #filter(num_total_mut < 10000) %>%
  ggplot(aes(tumor_stela, Signature.5.abs)) +
  geom_point(aes(color = cl_members), size=5) +
  ylab("Absolute signature 5 \n (signature x total mutations)") +
  xlab("Tumour TL (kb)") +
  geom_smooth(method = "lm", se = FALSE) +
  theme(panel.background = element_blank(),
        #axis.ticks = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.5),
        legend.background = element_rect(fill="lightblue",
                                         size=0.5, linetype="solid",
                                         colour ="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  stat_cor(method = "pearson", label.x = 4.5, label.y = 5000)


# for sklearn classification
#write_csv(all, "/path/to/project/data/train_clusters.csv")