library(ggplot2)
library(reshape2)
library(ggsignif)
library(pwr)
library(lsr)
library(ROCR)

df <- read.delim('Flow_data/SLE_Lachie_Results.txt')
df$Disease <- factor(df$Disease)
df$Run <- factor(df$Run)
df$disease_activity <- ''
df$disease_activity[11:20] <- ifelse(df$SLEDAI.2k.score[11:20] < 4, 'inactive', 'active')
df$disease_activity[1:10] <- 'healthy'
df$disease_activity <- factor(df$disease_activity)

MFI <- df[,c(7,9,10,11,12,13,15)]
colnames(MFI) <- gsub('_', ' ', colnames(MFI))
colnames(MFI) <- gsub('Geomean', '', colnames(MFI))
colnames(MFI) <- gsub('  ', ' ', colnames(MFI))

dir.create('normality_plots')

# Check assumptions of normality and homogeneity of variance for T-test
# Iterate over each column and print histogram and qqplot to check normality
# Save plots named by column name to normality_plots directory
for (i in 1:ncol(MFI)) {
  # Save histogram
  pdf(paste0('normality_plots/', gsub(' ', '_', colnames(MFI)[i]), '.histogram.pdf'))
  hist(MFI[,i], main=gsub('_|Geomean', '', colnames(MFI)[i]), xlab='MFI')
  dev.off()
  
  # Save QQ plot
  pdf(paste0('normality_plots/', gsub(' ', '_', colnames(MFI)[i]), '.qqnorm.pdf'))
  qqnorm(MFI[,i], main=gsub('_|Geomean', '', colnames(MFI)[i]))
  qqline(MFI[,i])
  dev.off()
}

normality.test <- lapply(MFI, function(x) shapiro.test(x)$p.value)
normality.test_log <- lapply(MFI, function(x) shapiro.test(log(x))$p.value)
normality.test <- data.frame(celltype=names(normality.test), p.value=unlist(normality.test), p.value.log=unlist(normality.test_log))
write.csv(normality.test, 'normality_test.csv', row.names=FALSE, quote=FALSE)

MFI_log <- log(MFI)

for (i in 1:ncol(MFI)) {
  # Save histogram
  pdf(paste0('normality_plots/', gsub(' ', '_', colnames(MFI_log)[i]), '.histogram_log.pdf'))
  hist(MFI_log[,i], main=gsub('_|Geomean', '', colnames(MFI_log)[i]), xlab='log(MFI)')
  dev.off()
  
  # Save QQ plot
  pdf(paste0('normality_plots/', gsub(' ', '_', colnames(MFI_log)[i]), '.qqnorm_log.pdf'))
  qqnorm(MFI_log[,i], main=gsub('_|Geomean', '', colnames(MFI_log)[i]))
  qqline(MFI_log[,i])
  dev.off()
}

MFI_log <- cbind(disease = df$Disease, Age = df$Age, 
activity = df$disease_activity, Run = df$Run, SLEDAI = df$SLEDAI.2k.score,
MFI_log)
MFI_log$Run <- factor(MFI_log$Run)
MFI_log$activity <- factor(MFI_log$activity)
MFI_log$disease <- factor(MFI_log$disease)

MFI_log_wide <- melt(MFI_log[,-5], id.vars = c('disease', 'Age', 'activity', 'Run'))

pdf('markers_boxplot_log.pdf')
comparisons <- list(c('SLE', 'Healthy'))
ggplot(MFI_log_wide, aes(x=disease, y=log(value), colour=disease)) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill=NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='t.test', color = 'black') +
  ylab('log(MFI)') + xlab('') +
  theme_minimal() +
  facet_wrap(~variable, scales = 'free_y', strip.position = "bottom") +
  theme(legend.position = "none")
dev.off()

# Test for differences with linear moodel
models <- lapply(split(MFI_log_wide, MFI_log_wide$variable), function(x){
    summary(lm(value ~ disease + Age + Run, data = x))
})
model_results <- lapply(models, function(x){
    data.frame(pvalue = x$coefficients[2,4], estimate = x$coefficients[2,1])
})
model_results <- do.call(rbind, model_results)
model_results <- cbind(marker=names(models), model_results)
model_results$padjust <- p.adjust(model_results$pvalue, method = 'bonferroni')
row.names(model_results) <- NULL
write.table(model_results, 'linear_model_results.txt', row.names=FALSE, quote=FALSE)

power <- lapply(models, function(x){
    r_squared <- x$r.squared
    f2 <- r_squared / (1 - r_squared)
    pwr.f2.test(u=3, models[[1]]$df[2], f2=f2, sig.level=0.05, power=NULL)$power
})
power <- data.frame(marker=names(power), power=unlist(power))

sample.size <- lapply(models, function(x){
    r_squared <- x$r.squared
    f2 <- r_squared / (1 - r_squared)
    pwr.f2.test(u=3, v=NULL, f2=f2, sig.level=0.05, power=0.8)$v + 4
})



power <- data.frame(marker=names(models), power=power$power, sample.size=unlist(sample.size))
write.csv(power, 'power_analysis.csv', row.names=FALSE, quote=FALSE)

library(ggpmisc)
# Plot model with two regression lines
pdf('TIMP1_regression.pdf')
ggplot(MFI_log_wide[MFI_log_wide$variable == "Tregs TIMP1",], aes(x = Age, y = value, color = disease)) +
  geom_point() +  # scatter plot
  geom_smooth(method = "lm", aes(group = disease), se = FALSE) +  # regression lines
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = 0.1) +  # Add equation, R², and p-value
  labs(title = "TIMP1 Expression in Tregs by Disease Status",
       x = "Age", y = "log(MFI)",
       color = "Disease Status") +
  theme_minimal()
dev.off()

pdf('PIM2_regression.pdf')
ggplot(MFI_log_wide[MFI_log_wide$variable == "BMEM PIM2",], aes(x = Age, y = value, color = disease)) +
  geom_point() +  # scatter plot
  geom_smooth(method = "lm", aes(group = disease), se = FALSE) +  # regression lines
  labs(title = "PIM2 Expression in BMEM by Disease Status",
       x = "Age", y = "log(MFI)",
       color = "Disease Status") +
  theme_minimal()
dev.off()

# Calculate cohens for TIMP1 and PIM2
cohens_d_timp1 <- cohensD(MFI_log$`Tregs TIMP1`[1:10], MFI_log$`Tregs TIMP1`[11:20])

ggplot

cohens_d_pim2 <- cohensD(MFI_log$`BMEM PIM2`[1:10], MFI_log$`BMEM PIM2`[11:20])

# Correlation between TIMP1 and PIM2 with SLEDAI
TIMP1_cor <- cor.test(MFI_log$`Tregs TIMP1`[11:20], MFI_log$SLEDAI[11:20])
pdf('TIMP1_SLEDAI_correlation.pdf')
ggplot(MFI_log[11:20,], aes(x = `Tregs TIMP1`, y = SLEDAI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "TIMP1 Expression in Tregs vs SLEDAI-2K Score",
  subtitle = paste("Correlation:", round(TIMP1_cor$estimate, 2), ';', "p-value:", round(TIMP1_cor$p.value, 3)),
       x = "log(MFI)", y = "SLEDAI-2K") +
  theme_minimal()
dev.off()

PIM2_cor <- cor.test(MFI_log$`BMEM PIM2`[11:20], MFI_log$SLEDAI[11:20])
pdf('PIM2_SLEDAI_correlation.pdf')
ggplot(MFI_log[11:20,], aes(x = `BMEM PIM2`, y = SLEDAI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "PIM2 Expression in BMEM vs SLEDAI-2K Score",
  subtitle = paste("Correlation:", round(PIM2_cor$estimate, 2), ';', "p-value:", round(PIM2_cor$p.value, 3)),
       x = "log(MFI)", y = "SLEDAI-2K") +
  theme_minimal()
dev.off()

### Plot ROC curve to show predictive ability of markers ###
MFI_log$disease <- ifelse(MFI_log$disease == 'SLE', 1, 0)
MFI_log$disease <- factor(MFI_log$disease)
MFI_log$BMEM_CellCount <- df$BMEM_CellCount
MFI_log$Tregs_CellCount <- df$Tregs_CellCount

# PIM2 - plot all feature combinations
PIM2_model <- glm(disease ~ `BMEM PIM2` + BMEM_CellCount + Age, data = MFI_log[-c(1,18),], family = "binomial")
probabilities <- predict(PIM2_model, type = "response")
pred <- prediction(probabilities, MFI_log$disease[-c(1,18)])
perf <- performance(pred, "tpr", "fpr")
pdf('PIM2_ROC.pdf')
plot(perf, col = "blue", main = "ROC for PIM2 + BMEM Cell Count + Age", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred, measure = "auc")
auc_value <- auc@y.values[[1]]
# Add AUC to the plot
text(0.6, 0.2, labels = paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)
dev.off()
perf <- performance(pred, "sens", "spec")
opt_cutoff <- which.max(perf@y.values[[1]] + perf@x.values[[1]])
sensitivity <- perf@y.values[[1]][opt_cutoff]
specificity <- perf@x.values[[1]][opt_cutoff]

# TIMP1
TIMP1_model <- glm(disease ~ `Tregs TIMP1` + Tregs_CellCount + Age, data = MFI_log, family = "binomial")
probabilities <- predict(TIMP1_model, type = "response")
pred <- prediction(probabilities, MFI_log$disease)
perf <- performance(pred, "tpr", "fpr")
pdf('TIMP1_ROC.pdf')
plot(perf, col = "blue", main = "ROC for TIMP1, Treg cell count, and Age", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred, measure = "auc")
auc_value <- auc@y.values[[1]]
# Add AUC to the plot
text(0.6, 0.2, labels = paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)
dev.off()

perf <- performance(pred, "sens", "spec")
opt_cutoff <- which.max(perf@y.values[[1]] + perf@x.values[[1]])
sensitivity <- perf@y.values[[1]][opt_cutoff]
specificity <- perf@x.values[[1]][opt_cutoff]

all_model <- glm(Disease ~ Tregs_Geomean_TIMP1 + Tregs_CellCount +  BMEM_Geomean_PIM2 + BMEM_CellCount + Age, data = df, family = "binomial")
probabilities <- predict(all_model, type = "response")
pred <- prediction(probabilities, df$Disease)
perf <- performance(pred, "tpr", "fpr")
plot(perf, col = "blue", main = "ROC for all features", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred, measure = "auc")
auc_value <- auc@y.values[[1]]
# Add AUC to the plot
text(0.6, 0.2, labels = paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)


test1 = unique(MFI_log_wide$variable)[7]

pdf('IL2RG_regression.pdf')
ggplot(MFI_log_wide[MFI_log_wide$variable == test1,], aes(x = Age, y = value, 
                                                          color = disease)) +
  geom_point() +  # scatter plot
  geom_smooth(method = "lm", aes(group = disease), se = FALSE) +  # regression lines
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = 0.1) +  # Add equation, R², and p-value
  labs(title = "IL2rG in CD4+ Tcells by Disease Status",
       x = "Age", y = "log(MFI)",
       color = "Disease Status") +
  theme_minimal()
dev.off()


