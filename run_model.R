

library(doMC)

run_model <- function(array_first, with_cov, plate_correction, control_gender, control_cancerdraw, remove_age_pc, remove_cancer,remove_cancer_var,remove_cancer_first,method_name,scale_data, read_file, output_file){
  
  if(is.null(scale_data)){
    scale_data_save <- 'no_scale'
  } else {
    scale_data_save <- 'center_scale'
  }
  if(array_first){
    array_save <- 'array_first'
  } else {
    array_save <- 'plate_first'
  }
  
  if(with_cov){
    cov_save <- 'with_cov'
  } else {
    cov_save <- 'no_cov'
  }

  # make string to read in correct file
  read_string <- paste0('data/',read_file, '/', array_save, '_', cov_save, '_', plate_correction, '.rda')
  
  # load train_valid_25 data
  dat <- readRDS(read_string)
  new_train <- dat[[1]]
  new_test <- dat[[2]]
  rm(dat)

  # load gene probes
  gene_probes <- read_csv('data/all_gene_cpg_loc.csv')
  gene_region <- paste('Body', collapse = '|')
  gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]
  gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])
  
  # replace features with gene_probes
  cg_names <- names(new_train)[grepl('^cg', names(new_train))]
  int_probes <- intersect(cg_names, gene_probes)
  new_train <- replace_features(temp_data = new_train, new_feats = int_probes)
  new_test <- replace_features(temp_data = new_test, new_feats = int_probes)
  
  # subset by lfs probes 
  new_train <- replace_features(temp_data = new_train, new_feats = keep_features)
  new_test <- replace_features(temp_data = new_test, new_feats = keep_features)
  
  if(remove_cancer_first){
    remove_cancer_first_save <- 'cancer_first'
    if(remove_cancer){
      registerDoMC(cores=1)
      bh_data <- remove_cancer_signature(new_train, new_test, beta_thresh = 0.05)
      new_train <- bh_data[[1]]
      new_test <- bh_data[[2]]
      remove_cancer_save <- 'remove_cancer'
    } else {
      remove_cancer_save <- 'keep_cancer'
    }
    
    if(remove_age_pc){
      new_train <- remove_pc(new_train, pc_loc = 1)
      new_test <- remove_pc(new_test, pc_loc = 1)
      remove_age_pc_save <- 'remove_age_pc'
    } else {
      remove_age_pc_save <- 'keep_age_pc'
      
    }
  } else {
    remove_cancer_first_save <- 'age_first'
    
    if(remove_age_pc){
      new_train <- remove_pc(new_train, pc_loc = 1)
      new_test <- remove_pc(new_test, pc_loc = 1)
      remove_age_pc_save <- 'remove_age_pc'
    } else {
      remove_age_pc_save <- 'keep_age_pc'
      
    }
    
    if(remove_cancer){
      bh_data <- remove_cancer_signature(new_train, new_test, beta_thresh = 0.05)
      new_train <- bh_data[[1]]
      new_test <- bh_data[[2]]
      remove_cancer_save <- 'remove_cancer'
      message('FINISHED REMOVING CANCER')
    } else {
      remove_cancer_save <- 'keep_cancer'
    }
  }
 
  # get data ready for models
  new_train$age_label <- factor(new_train$age_label, levels = c('positive', 'negative'))
  new_test$age_label <- factor(new_test$age_label, levels = c('positive', 'negative'))
  
  y_train <- new_train$age_label
  y_test <- new_test$age_label
  
  # get gender and cancer draw
  new_train <- cbind(as.data.frame(class.ind(new_train$gender)), 
                     new_train)
  new_test <- cbind(as.data.frame(class.ind(new_test$gender)), 
                    new_test)
  new_train$gender <- NULL
  new_test$gender <- NULL
  
  new_train <- cbind(as.data.frame(class.ind(new_train$cancer_atdraw)), 
                     new_train)
  new_test <- cbind(as.data.frame(class.ind(new_test$cancer_atdraw)), 
                    new_test)
  new_train$cancer_atdraw <- NULL
  new_test$cancer_atdraw <- NULL
  
  feat_names <- names(new_train)[grepl('^cg', names(new_train))]
  
  if(control_gender){
    feat_names <- c('F', 'M', feat_names)
    control_gender_save <- 'control_gender'
  } else {
    control_gender_save <- 'no_gender'
  }
  
  if(control_cancerdraw){
    feat_names <- c('No', 'Yes', feat_names)
    control_cancerdraw_save <- 'control_cancerdraw'
  }else {
    control_cancerdraw_save <- 'no_cancerdraw'
  }
  
  pheno_names <- names(new_train)[!grepl('^cg', names(new_train))]
  train_mat <- new_train[, feat_names]
  test_mat <- new_test[, feat_names]
  test_pheno <- new_test[, pheno_names]
  
  # create string for saving data
  output_string <- paste0(output_file, '/', array_save, '_', cov_save, '_', plate_correction, '_', remove_age_pc_save, '_', remove_cancer_save, '_', remove_cancer_first_save,'_',remove_cancer_var, '_' ,control_gender_save, '_', control_cancerdraw_save,'_',method_name,'_',scale_data_save, '.rda')
  
  model_output_string <- paste0(model_output_file, '/', array_save, '_', cov_save, '_', plate_correction, '_', remove_age_pc_save, '_', remove_cancer_save, '_', remove_cancer_first_save,'_',remove_cancer_var, '_' ,control_gender_save, '_', control_cancerdraw_save,'_',method_name,'_',scale_data_save, '.rda')
  registerDoMC(cores=8)
  
  message('STARTING ', method_name)
   mod_info <- run_classifier(training_data = train_mat, train_outcome = y_train,testing_data = test_mat,method_name = method_name,test_data = test_pheno,scale_data=scale_data)
  test_results <- mod_info[[1]]
  best_model <- mod_info[[2]]
  
  saveRDS(best_model, file=model_output_string)
  saveRDS(test_results, file=output_string)
}



source('helpers.R')
load(file = 'data/lfs_probes_0.05.rda')
# "array_first_with_cov_combat_keep_age_pc_keep_cancer_cancer_first_cancer_diagnosis_no_gender_control_cancerdraw_xgbTree.rda"
method_name <- c('rf', 'xgbTree', 'svmLinear', 'svmRadial', 'bayesglm', 'gamboost')
scale_data <- c(c('center', 'scale'))
array_first = c(TRUE)
with_cov = c(TRUE)
plate_correction = c('combat','pc_removal')
control_gender=c(FALSE)
control_cancerdraw=c(FALSE)
remove_age_pc =c(FALSE)# TRUE  (next one)
remove_cancer =c(TRUE)# FALSE (next one)
remove_cancer_first=c(TRUE)
remove_cancer_var = c('cancer_diagnosis')
read_file = 'train_valid_25'
output_file = 'train_valid_25_model'
model_output_file = 'train_valid_25_model_best'

for(i in array_first){
  for(j in with_cov){
    for(l in plate_correction){
      for(m in control_gender){
        for(n in control_cancerdraw){
          for(o in remove_age_pc){
            for(p in remove_cancer){
              for(q in remove_cancer_first){
                for(r in remove_cancer_var){
                  for(s in method_name){
                    for(t in scale_data){
                      run_model(array_first = i, 
                                with_cov = j, 
                                plate_correction = l, 
                                control_gender = m, 
                                control_cancerdraw = n, 
                                remove_age_pc = o, 
                                remove_cancer = p, 
                                remove_cancer_var=r,
                                remove_cancer_first = q,
                                method_name = s,
                                scale_data = scale_data,
                                read_file = read_file, 
                                output_file = output_file)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
# # # # #
source('helpers.R')
library(data.table)
the_file_name <- 'train_valid_25_model/'

file_names <- list.files(the_file_name)
file_names <- file_names[!grepl('best|old', file_names)]
dat_list <- list()
for(i in 1:length(file_names)){
  this_file <- file_names[i]
  temp <- read.csv(paste0(the_file_name, this_file))
  temp$age_label <- factor(temp$age_label, levels=c('positive', 'negative'))
  auc_value <- pROC::auc(temp$age_label, temp$preds)

  temp_results <- data_frame(key=this_file, value=auc_value)
  dat_list[[i]] <- temp_results
  print(i)
}

# combine results
temp <- do.call('rbind', dat_list)

temp <- temp[order(temp$value, decreasing = TRUE),]

# ggplot(temp, aes(key, value)) + geom_bar(stat='identity')
# get name of top results
# top_file <- temp$key[temp$value ==max(temp$value)]
top_file <- temp$key[7]
top_file
the_file_name <- 'train_valid_25_model/'

# top_file <- 'array_first_with_cov_combat_keep_age_pc_keep_cancer_cancer_first_cancer_diagnosis_control_gender_control_cancerdraw_xgbTree_center_scale.rda'


# top_file <- temp$key[3]
temp <- read.csv(paste0(the_file_name, top_file))
library(ROCR)
# top_file <- temp$key[44]
temp$age_label <- factor(temp$age_label, levels=c('positive', 'negative'))
auc_value <- pROC::auc(temp$age_label,temp$preds )
temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]

conmat_plot(data = temp_null,predict = 'preds', actual = 'age_label', cutoff = 0.5, get_plot = TRUE,other_title = '', data_type = 'null', text_name = 'age')
conmat_plot(data = temp_cancer,predict = 'preds', actual = 'age_label', cutoff = 0.5, get_plot = TRUE,other_title = '', data_type = 'cancer', text_name = 'age')

