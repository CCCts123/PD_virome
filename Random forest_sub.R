library(randomForest)  
library(caret)
library(pROC)
library(dplyr)

###### Set training control parameters for subgingival plaque
ctrl_smote <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                           classProbs = TRUE, summaryFunction = twoClassSummary,
                           sampling = "smote", savePredictions = FALSE)

ctrl_up_loocv <- trainControl(method = "LOOCV",
                              classProbs = TRUE, summaryFunction = twoClassSummary,
                              sampling = "up", savePredictions = FALSE)

feature_count <- 66
if (feature_count <= 15) {
  mtry_values <- 1:feature_count
} else {
  mtry_max <- floor(feature_count / 3)
  mtry_values <- unique(round(seq(1, mtry_max, length.out = 5)))
}
tune_grid <- expand.grid(mtry = mtry_values)

###### Select the correct trainControl based on the training set
get_train_control <- function(train_data) {
  n_train <- nrow(train_data)
  min_train <- min(table(train_data$Group))
  
  if (n_train >= 20 && min_train > 5) {
    cat("  -> Using SMOTE + repeated 5-fold CV (n_train=", n_train, ", min_class=", min_train, ")\n")
    return(ctrl_smote)
  } else {
    cat("  -> Using UP-SAMPLING + LOOCV (n_train=", n_train, ", min_class=", min_train, ")\n")
    return(ctrl_up_loocv)
  }
}

###### MDA calculation
dataset <- read.csv("sub_dataset.csv", row.names=1)
dataset$Group <- as.factor(dataset$Group)
set.seed(123)
rf_model <- train(Group ~ ., data = dataset, method = "rf", ntree = 2000, 
                  trControl = ctrl_smote, tuneGrid = tune_grid, metric = "ROC", importance = TRUE)
importance_scores <- importance(rf_model$finalModel)[, "MeanDecreaseAccuracy"]
importance_rank <- rank(-importance_scores)
importance_df <- data.frame(
  vOTU = names(importance_scores),
  MeanDecreaseAccuracy = importance_scores,
  Rank = importance_rank)
importance_df <- importance_df[order(importance_df$Rank), ]
write.csv(importance_df, "MDA_sub.csv", row.names = FALSE)

###### intra-cohort 
data1 <- read.csv("Izawa_2021.csv", row.names=1)
data2 <- read.csv("Belstrøm_2021.csv", row.names=1)
data3 <- read.csv("Moghadam_2023.csv", row.names=1)
data4 <- read.csv("Wang_2015.csv", row.names=1)

stratified_cv_auc <- function(data, group_col, n_repeats, max_k = 5) {  
  min_class <- min(table(data[[group_col]]))
  k_folds <- min(max_k, min_class)
  if (k_folds < max_k) {
    cat("Note: Reducing outer folds to", k_folds, "due to small minority class size.\n")
  }
  
  test_auc_matrix <- matrix(NA, nrow = n_repeats, ncol = k_folds)
  train_auc_matrix <- matrix(NA, nrow = n_repeats, ncol = k_folds)
  
  for (repeat_i in 1:n_repeats) {  
    set.seed(123 + repeat_i)  
    data <- data[sample(seq_len(nrow(data))), ]  
    folds <- createFolds(data[[group_col]], k = k_folds, list = TRUE, returnTrain = FALSE)
    
    for (fold_i in 1:k_folds) {  
      test_indices <- folds[[fold_i]]  
      train_indices <- setdiff(seq_len(nrow(data)), test_indices)  
      train_data <- data[train_indices, ]  
      test_data <- data[test_indices, ]  
      
      current_ctrl <- get_train_control(train_data)
      
      set.seed(123 + repeat_i * 10 + fold_i) 
      tryCatch({
        rf_model <- train(Group ~ ., data = train_data, method = "rf", ntree = 2000, 
                          trControl = current_ctrl, tuneGrid = tune_grid, metric = "ROC")  
        
        prob_test <- predict(rf_model, newdata = test_data, type = "prob")[, "periodontitis"]
        roc_test <- roc(test_data$Group, prob_test, quiet = TRUE)
        test_auc_matrix[repeat_i, fold_i] <- auc(roc_test)
        
        prob_train <- predict(rf_model, newdata = train_data, type = "prob")[, "periodontitis"]
        roc_train <- roc(train_data$Group, prob_train, quiet = TRUE)
        train_auc_matrix[repeat_i, fold_i] <- auc(roc_train)
        
      }, error = function(e) {
        cat("Error in repeat", repeat_i, "fold", fold_i, ":", e$message, "\n")
        test_auc_matrix[repeat_i, fold_i] <- NA
        train_auc_matrix[repeat_i, fold_i] <- NA
      })
    }  
  }  
  return(list(train = train_auc_matrix, test = test_auc_matrix))
}

n_repeats <- 5
for (data_name in c("data1", "data2", "data3", "data4")) {
  cat("\nProcessing", data_name, "\n")
  current_data <- get(data_name)
  cv_results <- stratified_cv_auc(current_data, "Group", n_repeats, max_k = 5)
  test_auc <- cv_results$test
  cat("Test AUC:\n"); print(test_auc)
  cat("Mean test AUC:", mean(test_auc, na.rm = TRUE), "\n")
}

###### cross-cohort
dataset <- read.csv("sub_dataset with author.csv", row.names=1)
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
  
  cat("\nTraining project:", projects[i], "\n")
  current_ctrl <- get_train_control(train_data)
  
  set.seed(123)
  rf_model <- train(x = train_x, y = train_y, method = "rf", metric = "ROC",
                    trControl = current_ctrl, tuneGrid = tune_grid, ntree = 2000)
  
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
data1 <- read.csv("Izawa_2021.csv", row.names=1)
data2 <- read.csv("Belstrøm_2021.csv", row.names=1)
data3 <- read.csv("Moghadam_2023.csv", row.names=1)
data4 <- read.csv("Wang_2015.csv", row.names=1)

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
  
  cat("LOCO training set: n_train=", nrow(train_data), ", min_class=", min(table(train_data$Group)), "\n")
  set.seed(123)
  rf_model <- train(Group ~ ., data = train_data, method = "rf",
                    trControl = ctrl_smote, tuneGrid = tune_grid,
                    ntree = 2000, metric = "ROC")
  
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
