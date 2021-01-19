source('helpers.R')
library(data.table)
library(ROCR)

# loop through validation results and stor plot in list
train_file_name <- 'data/train_valid_model/'
test_file_name <- 'data/train_test_model/'

train_names <- list.files(train_file_name)
train_names <- train_names[!grepl('best', train_names)]

test_names <- list.files(test_file_name)
test_names <- test_names[!grepl('best', test_names)]

valid_list <- list()
valid_test_list <- list()
for(i in 1:length(train_names)){
  train_file <- train_names[i]
  if(train_file %in% test_names){
    
    # valid result
    valid <- read.csv(paste0(train_file_name, train_file))
    valid$age_label <- factor(valid$age_label, levels=c('positive', 'negative'))
    auc_valid <- round(pROC::auc(valid$age_label,valid$preds ), 3)

    valid_null <- valid[valid$cancer_diagnosis=='Unaffected',]
    valid_cancer <- valid[valid$cancer_diagnosis!='Unaffected',]
    
    # test result
    test <- read.csv(paste0(test_file_name, train_file))
    test$age_label <- factor(test$age_label, levels=c('positive', 'negative'))
    auc_test <- round(pROC::auc(test$age_label,test$preds ), 3)
    test_null <- test[test$cancer_diagnosis=='Unaffected',]
    test_cancer <- test[test$cancer_diagnosis!='Unaffected',]
    
    # get validation plots
    valid_null_plot=conmat_plot(data = valid_null,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title =paste0(train_file, 'valid auc=', auc_valid ), data_type = 'null', text_name = 'age')
    valid_cancer_plot=conmat_plot(data = valid_cancer,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title = paste0(train_file, 'valid auc=', auc_valid ), data_type = 'cancer', text_name = 'age')
    
    # get test plots
    test_null_plot=conmat_plot(data = test_null,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title =paste0(train_file, 'test auc=', auc_test ), data_type = 'null', text_name = 'age')
    test_cancer_plot=conmat_plot(data = test_cancer,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title = paste0(train_file, 'test auc=', auc_test ), data_type = 'cancer', text_name = 'age')
    
    # combine valid and test data
    valid_plot_data <- list(valid_null_plot, valid_cancer_plot)
    test_plot_data <- list(test_null_plot, test_cancer_plot)
    
    valid_test_list[[i]] <- list(valid_plot_data, test_plot_data)
    
  
   
  } else {
    temp <- read.csv(paste0(train_file_name, train_file))
    temp$age_label <- factor(temp$age_label, levels=c('positive', 'negative'))
    auc_value <- round(pROC::auc(temp$age_label,temp$preds ), 3)
    
    temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
    temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]
    
    null_plot=conmat_plot(data = temp_null,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title =paste0(train_file, ' auc=', auc_value ), data_type = 'null', text_name = 'age')
    cancer_plot=conmat_plot(data = temp_cancer,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title = paste0(train_file, ' auc=', auc_value ), data_type = 'cancer', text_name = 'age')
    
    plot_data <- list(null_plot, cancer_plot)
    valid_list[[i]] <- plot_data
  }
  
}

# Another option: create pdf where each page is a separate plot.
pdf("plots_valid")
for (i in 1:length(valid_list)) {
  print(valid_list[[i]])
}
dev.off()

pdf("plots_test")
for (i in 1:length(valid_test_list)) {
  print(valid_test_list[[i]])
}
dev.off()


# add in sensitivity (true positive) and specificity (false positive)
# from this should be able to get false negative


source('helpers.R')
# new data
# read in good results
# logit_array_first_no_cov_pc_removal_keep_cancer_lfs_first_control_gender_scale_data_no_draw
temp <- plot_results(folder_name = 'data/train_valid_model/', file_name ='logit_array_first_no_cov_pc_removal_keep_cancer_lfs_first_control_gender_scale_data_no_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

temp <- plot_results(folder_name = 'data/train_test_model/', file_name ='logit_array_first_no_cov_pc_removal_keep_cancer_lfs_first_control_gender_scale_data_no_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])


# rf_array_first_with_cov_pc_removal_remove_cancer_lfs_first_control_gender_no_scale_no_draw.csv
temp <- plot_results(folder_name = 'data/train_valid_model/', file_name ='rf_array_first_with_cov_pc_removal_remove_cancer_lfs_first_control_gender_no_scale_no_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

temp <- plot_results(folder_name = 'data/train_test_model/', file_name ='rf_array_first_with_cov_pc_removal_remove_cancer_lfs_first_control_gender_no_scale_no_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

# svm_array_first_no_cov_combat_keep_cancer_lfs_first_no_gender_no_scale_control_draw.csv
temp <- plot_results(folder_name = 'data/train_valid_model/', file_name ='svm_array_first_no_cov_combat_keep_cancer_lfs_first_no_gender_no_scale_control_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

temp <- plot_results(folder_name = 'data/train_test_model/', file_name ='svm_array_first_no_cov_combat_keep_cancer_lfs_first_no_gender_no_scale_control_draw.csv', cutoff_value = .25, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

# xgb_array_first_no_cov_combat_keep_cancer_lfs_first_control_gender_scale_data_control_draw.csv
temp <- plot_results(folder_name = 'data/train_valid_model/', file_name ='xgb_array_first_no_cov_combat_keep_cancer_lfs_first_control_gender_scale_data_control_draw.csv', cutoff_value = .125, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

temp <- plot_results(folder_name = 'data/train_test_model/', file_name ='xgb_array_first_no_cov_combat_keep_cancer_lfs_first_control_gender_scale_data_control_draw.csv', cutoff_value = .125, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])


# old data
#svm_array_first_no_cov_combat_control_gender_scale_data.csv
# svm_array_first_with_cov_pc_removal_control_gender_no_scale

temp <- plot_results(folder_name = 'data/old/train_valid_25_model/', file_name ='svm_array_first_with_cov_pc_removal_control_gender_no_scale.csv', cutoff_value = .125, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

temp <- plot_results(folder_name = 'data/old/train_test_25_model/', file_name ='svm_array_first_with_cov_pc_removal_control_gender_no_scale.csv', cutoff_value = .125, plot_var = 'age')
print(temp[[1]])
print(temp[[2]])

