library(survival)
library(survminer)

pred <- read_csv("/path/to/project/data/clf_mode_pred_cl_w_mode_n3.csv")
pred <- rename(pred, cl_members=mode_pred)

pred_clin <- left_join(clin %>% filter(!tumor_db %in% all$tumor_db), pred) 

surv_data_clf <- rbind(pred_clin %>% select(intersect(colnames(all), colnames(pred_clin))), all %>% select(intersect(colnames(all), colnames(pred_clin))))


nrow(all %>% select(cl_members) %>% filter(cl_members == '1a')) #16
nrow(all %>% select(cl_members) %>% filter(cl_members == '1b')) #12
nrow(all %>% select(cl_members) %>% filter(cl_members == '2')) #16

nrow(surv_data_clf %>% select(cl_members) %>% filter(cl_members == '1a')) #38
nrow(surv_data_clf %>% select(cl_members) %>% filter(cl_members == '1b')) #82
nrow(surv_data_clf %>% select(cl_members) %>% filter(cl_members == '2')) #95


surv_data_clf %>% filter(is.na(cl_members)) %>% select(`AGE AT OPERATION`, `NPI SCORE`, tumor_stela, cl_members)
print(surv_data_clf %>% select(tumor_db, tumor_stela, `AGE AT OPERATION`, `NPI SCORE`, cl_members), n = 216)
print(pred_clin %>% select(tumor_db, tumor_stela, `AGE AT OPERATION`, `NPI SCORE`, cl_members), n = 216)
print(clin %>% select(tumor_db, `AGE AT OPERATION`, `NPI SCORE`), n = 216)
print(all %>% select(tumor_db, tumor_stela, `AGE AT OPERATION`, `NPI SCORE`, cl_members))


#use diagnosis to followup
surv_data_clf <- surv_data_clf %>% filter(!is.na(cl_members))
surv_obj_clf <- Surv(time = surv_data_clf$`Operation to follow/death (days)`, event = surv_data_clf$`Censoring Status`)
surv_fit_clf <- survfit(surv_obj_clf ~ cl_members, data = surv_data_clf)

ggsurvplot(surv_fit_clf, data = surv_data_clf, pval = TRUE, pval.method = TRUE, conf.int = TRUE)


surv_colours <- gg_color_hue(3)

#1a 2
surv_data_clf_1a_2 <- surv_data_clf %>% filter(cl_members %in% c("1a", "2"))
surv_obj_clf_1a_2 <- Surv(time = surv_data_clf_1a_2$`Operation to follow/death (days)`, event = surv_data_clf_1a_2$`Censoring Status`)
surv_fit_clf_1a_2 <- survfit(surv_obj_clf_1a_2 ~ cl_members, data = surv_data_clf_1a_2)

ggsurvplot(surv_fit_clf_1a_2, data = surv_data_clf_1a_2, pval = TRUE, pval.method = TRUE, conf.int = TRUE, palette = surv_colours[c(1, 3)])

#1b 2
surv_data_clf_1b_2 <- surv_data_clf %>% filter(cl_members %in% c("1b", "2"))
surv_obj_clf_1b_2 <- Surv(time = surv_data_clf_1b_2$`Operation to follow/death (days)`, event = surv_data_clf_1b_2$`Censoring Status`)
surv_fit_clf_1b_2 <- survfit(surv_obj_clf_1b_2 ~ cl_members, data = surv_data_clf_1b_2)

ggsurvplot(surv_fit_clf_1b_2, data = surv_data_clf_1b_2, pval = TRUE, pval.method = TRUE, conf.int = TRUE, palette = surv_colours[c(2, 3)])

#1a 1b
surv_data_clf_1a_1b <- surv_data_clf %>% filter(cl_members %in% c("1a", "1b"))
surv_obj_clf_1a_1b <- Surv(time = surv_data_clf_1a_1b$`Operation to follow/death (days)`, event = surv_data_clf_1a_1b$`Censoring Status`)
surv_fit_clf_1a_1b <- survfit(surv_obj_clf_1a_1b ~ cl_members, data = surv_data_clf_1a_1b)

ggsurvplot(surv_fit_clf_1a_1b, data = surv_data_clf_1a_1b, pval = TRUE, pval.method = TRUE, conf.int = TRUE, palette = surv_colours[c(1, 2)])


# length(intersect(clin$Mean, all$tumor_stela))
# nrow(clin) - nrow(clin %>% filter(is.na(Mean)))
# length(intersect(round(surv_data_clf$tumor_stela, 3), round(clin$Mean, 3)))
# print(surv_data_clf %>% filter(!tumor_stela %in% intersect(surv_data_clf$tumor_stela, clin$Mean) ) %>% select(tumor_db, tumor_stela, `AGE AT OPERATION`, `NPI SCORE`, cl_members), n=100)
# all %>% filter(tumor_db == "DB165") %>% select(tumor_db, tumor_stela, `AGE AT OPERATION`, `NPI SCORE`, cl_members)
# all$`Time to Death (days)`


pc_clf <- surv_data_clf %>%
  ggplot(aes(cl_members, fill=cl_members)) +
  geom_histogram(stat="count") +
  #xlabs("Predicted cluster")
  #scale_y_continuous(expand = c(0,1)) +
  theme(axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("Predicted cluster") +
  ylab("Count")

tl_clf <- ggplot(surv_data_clf, aes(tumor_stela, fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="TTL (kb)") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

np_clf <- ggplot(surv_data_clf, aes(surv_data_clf[['NPI SCORE']], fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="NPI Score") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

ag_clf<- ggplot(surv_data_clf, aes(surv_data_clf[['AGE AT OPERATION']], fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="Age (years)") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

grid.arrange(ag_clf, tl_clf, np_clf, pc_clf)

ggsurvplot(surv_fit_clf, data = surv_data_clf, pval = TRUE, pval.method = TRUE, conf.int = TRUE)


max(pred_clin %>% filter(cl_members == "1a") %>% select(tumor_stela))
max(pred_clin %>% filter(cl_members == "1b") %>% select(tumor_stela))
max(pred_clin %>% filter(cl_members == "2") %>% select(tumor_stela))

min(pred_clin %>% filter(cl_members == "1a") %>% select(tumor_stela))
min(pred_clin %>% filter(cl_members == "1b") %>% select(tumor_stela))
min(pred_clin %>% filter(cl_members == "2") %>% select(tumor_stela))

mean(pred_clin %>% filter(cl_members == "1a") %>% select(tumor_stela) %>% pull())
mean(pred_clin %>% filter(cl_members == "1b") %>% select(tumor_stela) %>% pull())
mean(pred_clin %>% filter(cl_members == "2") %>% select(tumor_stela) %>% pull())

min(pred_clin %>% filter(cl_members == "2") %>% select(`AGE AT OPERATION`) %>% pull())
pred_clin %>% filter(`AGE AT OPERATION` == 41) %>% filter(cl_members == "1a") %>% select(tumor_stela, `AGE AT OPERATION`, `NPI SCORE`)
mean(pred_clin %>% filter(cl_members == "1b") %>% select(`AGE AT OPERATION`) %>% pull())
mean(pred_clin %>% filter(cl_members == "2") %>% select(`AGE AT OPERATION`) %>% pull())

colnames(pred_clin)

for (m in pred_clin %>% select(contains("pred")) %>% colnames()) {
  print(m)
  print(nrow(pred_clin[which(pred_clin[[as.character(m)]] == pred_clin$cl_members), ]))
}

(160 + 153 + 165 + 133 + 133)/5
(161 + 154 + 164 + 139 + 132)/5

test_normal <- function(a, b, var){
  t.test(pred_clin %>% filter(cl_members == a) %>% select(var) %>% pull(),
         pred_clin %>% filter(cl_members == b) %>% select(var) %>% pull())
}

test_non_normal <- function(a, b, var){
  wilcox.test(pred_clin %>% filter(cl_members == a) %>% select(var) %>% pull(),
         pred_clin %>% filter(cl_members == b) %>% select(var) %>% pull())
}

test_normal("1a", "1b", "AGE AT OPERATION")
test_normal("1a", "2", "AGE AT OPERATION")
test_normal("1b", "2", "AGE AT OPERATION")

# Normally

shapiro.test(pred_clin %>% filter(cl_members == "1a") %>% select(`AGE AT OPERATION`) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "1b") %>% select(`AGE AT OPERATION`) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "2") %>% select(`AGE AT OPERATION`) %>% pull())

test_non_normal("1a", "1b", "AGE AT OPERATION")
test_non_normal("1a", "2", "AGE AT OPERATION")
test_non_normal("1b", "2", "AGE AT OPERATION")

test_non_normal("1a", "1b", "tumor_stela")
test_non_normal("1a", "2", "tumor_stela")
test_non_normal("1b", "2", "tumor_stela")

test_non_normal("1a", "1b", "NPI SCORE")
test_non_normal("1a", "2", "NPI SCORE")
test_non_normal("1b", "2", "NPI SCORE")

# Not normally

shapiro.test(pred_clin %>% filter(cl_members == "1a") %>% select(tumor_stela) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "1b") %>% select(tumor_stela) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "2") %>% select(tumor_stela) %>% pull())

shapiro.test(pred_clin %>% filter(cl_members == "1a") %>% select(`NPI SCORE`) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "1b") %>% select(`NPI SCORE`) %>% pull())
shapiro.test(pred_clin %>% filter(cl_members == "2") %>% select(`NPI SCORE`) %>% pull())







# 2 clusters

k <- 2

test <- hclust_apply(k, cl_metric, variables_for_hclust, binary_variables, add_cl_all = TRUE)


pred_2 <- read_csv("/path/to/project/data/clf_mode_pred_cl_n2.csv")
pred_2 <- rename(pred_2, cl_members=mode_pred)

pred_clin_2 <- left_join(clin %>% filter(!tumor_db %in% all$tumor_db), pred_2) 

surv_data_clf_2 <- rbind(pred_clin_2 %>% select(intersect(colnames(all), colnames(pred_clin_2))), all %>% select(intersect(colnames(all), colnames(pred_clin_2))))

surv_data_clf_2 %>% filter(is.na(cl_members)) %>% select(`AGE AT OPERATION`, `NPI SCORE`, tumor_stela, cl_members)

surv_data_clf_2 <- surv_data_clf_2 %>% filter(!is.na(cl_members))
surv_obj_clf_2 <- Surv(time = surv_data_clf_2$`Operation to follow/death (days)`, event = surv_data_clf_2$`Censoring Status`)
surv_fit_clf_2 <- survfit(surv_obj_clf_2 ~ cl_members, data = surv_data_clf_2)

ggsurvplot(surv_fit_clf_2, data = surv_data_clf_2, pval = TRUE, pval.method = TRUE, conf.int = TRUE)


pc_clf_n2 <- surv_data_clf_2 %>%
  ggplot(aes(cl_members, fill=cl_members)) +
  geom_histogram(stat="count") +
  #xlabs("Predicted cluster")
  #scale_y_continuous(expand = c(0,1)) +
  theme(axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), #element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(), #element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("Predicted cluster") +
  ylab("Count")

tl_clf_n2 <- ggplot(surv_data_clf_2, aes(tumor_stela, fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="TTL (kb)") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

np_clf_n2 <- ggplot(surv_data_clf_2, aes(surv_data_clf[['NPI SCORE']], fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="NPI Score") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

ag_clf_n2 <- ggplot(surv_data_clf_2, aes(surv_data_clf[['AGE AT OPERATION']], fill = cl_members)) +
  geom_boxplot() +
  coord_flip() +
  labs(title="Age (years)") +
  xlab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),#element_line(size = 0.4, color = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, color = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~cl_members)

grid.arrange(ag_clf_n2, tl_clf_n2, np_clf_n2, pc_clf_n2)

for (m in pred_clin_2 %>% select(contains("pred")) %>% colnames()) {
  print(m)
  print(nrow(pred_clin_2[which(pred_clin_2[[as.character(m)]] == pred_clin_2$cl_members), ]))
}

(170 + 172 + 159 + 171 + 150)/5
(167 + 165 + 164 + 111 + 154)/5




# stratification by single variables

pred_clin <- mutate(pred_clin, strat_npi =  case_when(pred_clin$`NPI SCORE` > 4.4 ~ 1, 
                                                      pred_clin$`NPI SCORE` <= 4.4 ~ 2))

pred_clin <- mutate(pred_clin, strat_ttl =  case_when(pred_clin$Mean <= 3.81 ~ 1, 
                                                      pred_clin$Mean > 3.81 ~ 2))

# pred_clin$strat_npi
# pred_clin$strat_ttl                    

surv_data_strat_npi <- pred_clin
surv_data_strat_npi <- surv_data_strat_npi %>% filter(strat_npi %in% c("1", "2"))
surv_obj_strat_npi <- Surv(time = surv_data_strat_npi$`Operation to follow/death (days)`, event = surv_data_strat_npi$`Censoring Status`)
surv_fit_strat_npi <- survfit(surv_obj_strat_npi ~ strat_npi, data = surv_data_strat_npi)

ggsurvplot(surv_fit_strat_npi, data = surv_data_strat_npi, pval = TRUE, pval.method = TRUE, conf.int = TRUE)


surv_data_strat_ttl <- pred_clin
surv_data_strat_ttl <- surv_data_strat_ttl %>% filter(strat_ttl %in% c("1", "2"))
surv_obj_strat_ttl <- Surv(time = surv_data_strat_ttl$`Operation to follow/death (days)`, event = surv_data_strat_ttl$`Censoring Status`)
surv_fit_strat_ttl <- survfit(surv_obj_strat_ttl ~ strat_ttl, data = surv_data_strat_ttl)

ggsurvplot(surv_fit_strat_ttl, data = surv_data_strat_ttl, pval = TRUE, pval.method = TRUE, conf.int = TRUE)

nrow(pred_clin %>% filter(strat_npi == 1))
nrow(pred_clin %>% filter(strat_npi == 2))
nrow(pred_clin %>% filter(strat_ttl == 1))
nrow(pred_clin %>% filter(strat_ttl == 2))


