

library(caret)
library(leaps)

setwd("/media/benni/Data/Benni_work/projects/scCODA_analysis/benchmark_power")

df = read.csv2("./results_grouped_fdr_01.csv", sep=",", header = TRUE, dec = ".")



df2 = df[c("n_controls", "n_cases", "Increase", "log.fold.increase", "mcc")]
df2["total"] =  df$n_controls+df$n_cases
df2["log.total"] =  log(df$n_controls+df$n_cases)
df2["mcc_scaled"] = (df$mcc + 1.0) /2.0
df2["log.increase"] = log(df2["Increase"])
mscale_min = min(df2["mcc_scaled"])
#empirical logit transform based on: https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.1890/10-0340.1
df2["mcc_trans"] = log((df2["mcc_scaled"]+0.00001)/
                            (1-df2["mcc_scaled"]+0.00001))

df2 = df[c("n_controls", "n_cases", "Increase", "log.fold.increase", "tpr")]
df2["total"] =  df$n_controls+df$n_cases
df2["log.control"] = log(df$n_controls)
df2["log.case"] = log(df$n_case)
df2["log.total"] =  log(df$n_controls+df$n_cases)
df2["log.ratio"] =  log(df$n_cases/df$n_controls)
df2["mcc_scaled"] = (df$mcc + 1.0) /2.0
df2["log.increase"] = log(df2["Increase"])

#empirical logit transform based on: https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.1890/10-0340.1
df2["tpr_trans"] = log((df2["tpr"]+0.00001)/
                            (1-df2["tpr"]+0.00001))


df3 = df2[c( "log.total", "log.increase", "log.fold.increase", "tpr")]
set.seed(666)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "repeatedcv", number = 10, repeats=5)
# Train the model
step.model <- train(tpr ~ log.total  + log.fold.increase+ log.total:log.increase+log.increase:log.fold.increase, 
                    data = df3,
                    family=quasibinomial('logit'),
                    method = "glm", 
                    #tuneGrid = data.frame(nvmax = 1:4),
                    trControl = train.control
)
step.model$results
summary(step.model$finalModel)
coef(step.model$finalModel, step.model$bestTune[[1]])