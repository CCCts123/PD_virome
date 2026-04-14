library(randomForest)  
library(caret)
library(pROC)
library(dplyr)

###### Set training control parameters for saliva
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                              classProbs = TRUE, summaryFunction = twoClassSummary, sampling = "smote", savePredictions = FALSE)
feature_count <- 105   
if (feature_count <= 15) {
  mtry_values <- 1:feature_count
} else {
  mtry_max <- floor(feature_count / 3)
  mtry_values <- unique(round(seq(1, mtry_max, length.out = 5)))
}
tune_grid <- expand.grid(mtry = mtry_values)

###### MDA calculation
dataset <- read.csv("saliva_dataset.csv", row.names=1)
dataset$Group <- as.factor(dataset$Group)
set.seed(123)
rf_model <- train(Group ~ ., data = dataset, method = "rf", ntree = 2000, 
                  trControl = train_control, tuneGrid = tune_grid, metric = "ROC", importance = TRUE)
importance_scores <- importance(rf_model$finalModel)[, "MeanDecreaseAccuracy"]
importance_rank <- rank(-importance_scores)
importance_df <- data.frame(
  vOTU = names(importance_scores),
  MeanDecreaseAccuracy = importance_scores,
  Rank = importance_rank)
importance_df <- importance_df[order(importance_df$Rank), ]
write.csv(importance_df, "MDA_saliva.csv", row.names = FALSE)

###### intra-cohort
data1 <- read.csv("Belstrøm_2017.csv", row.names=1)
data2 <- read.csv("Belstrøm_2021.csv", row.names=1)
data3 <- read.csv("Moghadam_2023.csv", row.names=1)
data4 <- read.csv("Saito_2022.csv", row.names=1)

k_folds <- 5  
n_repeats <- 5

stratified_cv_auc <- function(data, group_col, k_folds, n_repeats) {  
  results <- list()  
  for (repeat_i in 1:n_repeats) {  
    set.seed(123 + repeat_i)  
    data <- data[sample(seq_len(nrow(data))), ]  
    folds <- createFolds(data[[group_col]], k = k_folds, list = TRUE, returnTrain = FALSE)
  
    repeat_results <- numeric(k_folds)  
    for (fold_i in 1:k_folds) {  
      test_indices <- folds[[fold_i]]  
      train_indices <- setdiff(seq_len(nrow(data)), test_indices)  
      train_data <- data[train_indices, ]  
      test_data <- data[test_indices, ]  
      
      set.seed(123 + repeat_i * 10 + fold_i) 
      auc_val <- tryCatch({
        rf_model <- train(Group ~ ., data = train_data, method = "rf", ntree = 2000, 
                          trControl = train_control, tuneGrid = tune_grid, metric = "ROC")  
        prob_pred <- predict(rf_model, newdata = test_data, type = "prob")[, "periodontitis"]
        roc_obj <- roc(test_data$Group, prob_pred, quiet = TRUE)  
        auc(roc_obj)
      }, error = function(e) {
        NA  
      })
      repeat_results[fold_i] <- auc_val  
    }  
    results[[repeat_i]] <- repeat_results  
  }  
  return(results)  
}


cv_results <- stratified_cv_auc(data1, "Group", k_folds, n_repeats) 
auc_matrix <- do.call(rbind, cv_results)
print(auc_matrix) 
mean_auc <- colMeans(auc_matrix) 
overall_mean_auc <- mean(mean_auc)
print(paste("Overall mean AUC:", overall_mean_auc))

cv_results <- stratified_cv_auc(data2, "Group", k_folds, n_repeats) 
auc_matrix <- do.call(rbind, cv_results)
print(auc_matrix) 
mean_auc <- colMeans(auc_matrix) 
overall_mean_auc <- mean(mean_auc)
print(paste("Overall mean AUC:", overall_mean_auc))

cv_results <- stratified_cv_auc(data3, "Group", k_folds, n_repeats) 
auc_matrix <- do.call(rbind, cv_results)
print(auc_matrix) 
mean_auc <- colMeans(auc_matrix) 
overall_mean_auc <- mean(mean_auc)
print(paste("Overall mean AUC:", overall_mean_auc))

cv_results <- stratified_cv_auc(data4, "Group", k_folds, n_repeats) 
auc_matrix <- do.call(rbind, cv_results)
print(auc_matrix) 
mean_auc <- colMeans(auc_matrix) 
overall_mean_auc <- mean(mean_auc)
print(paste("Overall mean AUC:", overall_mean_auc))

###### cross-cohort
dataset <- read.csv("saliva_dataset with author.csv", row.names=1)
dataset$Group <- as.factor(dataset$Group)
projects <- unique(dataset$author)

results <- data.frame(train_project = character(),
                      test_project = character(),
                      auc = numeric(),
                      stringsAsFactors = F)

for (i in seq_along(projects)) {
  train_idx <- dataset$author == projects[i]
  train_data <- dataset[train_idx, ]
  train_x <- train_data[, !(names(train_data) %in% c("Group", "author"))]
  train_y <- train_data$Group
  
  set.seed(123)
  rf_model <- train(x = train_x, y = train_y, method = "rf", metric = "ROC", trControl = train_control, tuneGrid = tune_grid, ntree = 2000)
  
  for (j in seq_along(projects)) {
    if (i == j) next
    test_idx <- dataset$author == projects[j]
    test_data <- dataset[test_idx, ]
    test_x <- test_data[, !(names(test_data) %in% c("Group", "author"))]
    test_y <- test_data$Group
    
    prob <- predict(rf_model, test_x, type = "prob")[, "periodontitis"]
    roc_obj <- roc(test_y, prob)
    auc_val <- auc(roc_obj)
    
    results <- rbind(results, data.frame(
      train_project = projects[i],
      test_project = projects[j],
      auc = auc_val
    ))
    
    cat("Train:", projects[i], "Test:", projects[j], "AUC =", auc_val, "\n")
  }
}

write.csv(results, "cross_cohort_auc.csv", row.names = FALSE)

###### LOCO
data1 <- read.csv("Belstrøm_2017.csv", row.names=1)
data2 <- read.csv("Belstrøm_2021.csv", row.names=1)
data3 <- read.csv("Moghadam_2023.csv", row.names=1)
data4 <- read.csv("Saito_2022.csv", row.names=1)

data1$Group <- factor(data1$Group, levels = c("healthy", "periodontitis"))
data2$Group <- factor(data2$Group, levels = c("healthy", "periodontitis"))
data3$Group <- factor(data3$Group, levels = c("healthy", "periodontitis"))
data4$Group <- factor(data4$Group, levels = c("healthy", "periodontitis"))

loco_results <- list()
roc_objects <- list()
confusion_matrices <- list() 

for (i in 1:4) {
  cat("\nProcessing LOCO for dataset", i, "\n")
  train_list <- lapply(1:4, function(x) if (x != i) get(paste0("data", x)) else NULL)
  train_data <- do.call(rbind, train_list)
  test_data <- get(paste0("data", i))
  set.seed(123)
  rf_model <- train(Group ~ ., data = train_data, method = "rf", trControl = train_control, tuneGrid = tune_grid, ntree = 2000, metric = "ROC")
  prob_pred <- predict(rf_model, test_data, type = "prob")[, "periodontitis"]
  class_pred <- predict(rf_model, test_data)
  roc_obj <- roc(test_data$Group, prob_pred, levels = rev(levels(test_data$Group)))
  auc_val <- auc(roc_obj)
  cm <- confusionMatrix(class_pred, test_data$Group, positive = "periodontitis")
  metrics <- c(
    AUC = auc_val,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    F1 = cm$byClass["F1"])
  loco_results[[i]] <- metrics
  roc_objects[[i]] <- roc_obj
  confusion_matrices[[i]] <- cm$table
}

results_df <- do.call(rbind, loco_results)
rownames(results_df) <- paste0("Dataset", 1:4)
print(results_df)
write.csv(results_df, "LOCO_performance_metrics.csv", row.names = TRUE)

for (i in 1:4) {
  cm_file <- paste0("LOCO_confusion_matrix_dataset", i, ".csv")
  write.csv(confusion_matrices[[i]], cm_file, row.names = TRUE)
}
