
library(plyr)
library(dplyr)
library(tidyverse)
library(caret)
library(recipes)
library(devtools)
library(data.table)
load_all("useful.stuff.aa")
library(pROC)
#fold="fold10"
#read in the cohort
train<-fread("Cirrhosis_Train.csv",header=T) %>% select(-V1)
test<-fread("Cirrhosis_Test.csv",header=T) %>% select(-V1)
train_meta<-train
test_meta<-test

glmnetGrid <- expand.grid(
    alpha = 1,
    lambda = 10^seq(-5, -1, length.out = 100))

set.seed(1234)
ctrl_all <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 1,
                     verboseIter = TRUE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(train$type, 5, 10),
                     summaryFunction = twoClassSummary)

#######Train all the SSLs
model_name<-"ARTEMIS_DELFI"
recipe_seq <- recipe(type ~ ., data=train) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(-id,-starts_with("ratio"),-type,-starts_with("zscore"),-starts_with("class_"),-starts_with("H3"),-starts_with("H4"),-starts_with("states")) %>%
    step_pca(starts_with("ratio"), prefix = "ratio_pc_",threshold=.9)     %>%
    step_pca(starts_with("H3K27me3"), prefix = "H3K27me3_pc_",threshold=.9)     %>%
    step_pca(starts_with("H3K36me3"), prefix = "H3K36me3_pc_",threshold=.9)     %>%
    step_pca(starts_with("H3K9me3"), prefix = "H3K9me3_pc_",threshold=.9)     %>%
    step_pca(starts_with("H4K20me1"), prefix = "H4K20me1_pc_",threshold=.9)     %>%
    step_pca(starts_with("states_1_5"), prefix = "states_1_5_pc_",threshold=.9)     %>%
    step_pca(starts_with("states_10_13"), prefix = "states_10_13_pc_",threshold=.9)     %>%
    step_pca(starts_with("states_7_9"), prefix = "states_7_9_pc_",threshold=.9)     %>%
    step_pca(starts_with("class_SINE"), prefix = "class_SINE_pc_",threshold=.9)     %>%
    step_pca(starts_with("class_RNA_DNA"), prefix = "class_RNA_DNA_pc_",threshold=.9)     %>%
    step_pca(starts_with("class_LINE"), prefix = "class_LINE_pc_",threshold=.9)     %>%
    step_pca(starts_with("class_LTR"), prefix = "class_LTR_pc_",threshold=.9)     %>%
    step_pca(starts_with("class_Satellite"), prefix = "class_Satellite_pc_",threshold=.9)     %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())


colnames(bake(prep(recipe_seq),train))
set.seed(1234)
model1 <- caret::train(recipe_seq,
                          data = train,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)


train_preds<-get_cv_preds(train,model1)
test_preds<-get_test_preds(test,model1)
train_preds$model<-model_name
test_preds$model<-model_name
preds_cv<-train_preds
preds_val<-test_preds
saveRDS(model1,"LCr_ARTEMIS_DELFI.rds"))



write.csv(preds_cv,"LCr_Cross_Validated_Scores.csv")
write.csv(preds_val,"LCr_Test_Set_Scores.csv")

