

library(doMC)

run_model <- function(array_first, with_cov, plate_correction,remove_cancer, remove_cancer_first,read_file, output_file){
  
 
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
  read_string <- paste0(read_file, '/', array_save, '_', cov_save, '_', plate_correction, '.rda')
  
  message('loading ', read_string)
  # load train_valid_25 data
  dat <- readRDS(read_string)
  new_train <- dat[[1]]
  new_test <- dat[[2]]
  rm(dat)
  
  # load gene probes
  gene_probes <- read_csv('data/probe_data/all_gene_cpg_loc.csv')
  gene_region <- paste('Body', collapse = '|')
  gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]
  gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])
  
  # replace features with gene_probes
  cg_names <- names(new_train)[grepl('^cg', names(new_train))]
  int_probes <- intersect(cg_names, gene_probes)
  new_train <- replace_features(temp_data = new_train, new_feats = int_probes)
  new_test <- replace_features(temp_data = new_test, new_feats = int_probes)
  
  if(remove_cancer_first){
    cancer_first_save <- 'cancer_first'
    if(remove_cancer){
      bh_data <- remove_cancer_signature(new_train, new_test, beta_thresh = 0.05)
      new_train <- bh_data[[1]]
      new_test <- bh_data[[2]]
      remove_cancer_save <- 'remove_cancer'
      message('FINISHED REMOVING CANCER')
    } else {
      remove_cancer_save <- 'keep_cancer'
    }
    
    # subset by lfs probes 
    new_train <- replace_features(temp_data = new_train, new_feats = keep_features)
    new_test <- replace_features(temp_data = new_test, new_feats = keep_features)
    
  } else {
    cancer_first_save <- 'lfs_first'
    # subset by lfs probes 
    new_train <- replace_features(temp_data = new_train, new_feats = keep_features)
    new_test <- replace_features(temp_data = new_test, new_feats = keep_features)
    
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
  
  # create string for saving data
  train_string <- paste0(output_file, '/','full_train_', array_save, '_', cov_save, '_', plate_correction, '_', remove_cancer_save,'_',cancer_first_save, '.csv')
  test_string <- paste0(output_file, '/','test_', array_save, '_', cov_save, '_', plate_correction, '_', remove_cancer_save,'_',cancer_first_save, '.csv')
  
  
  # save training data
  write.csv(new_test, file=test_string)
  message('finished saving ', test_string)
  
  # save testing data
  write.csv(new_train, file=train_string)
  message('finished saving ', train_string)
  
}

#NEXT MAKE THE CODE BELOW FIR THE CODE ABOVE (NOT THAT MANY OPTIONS)

source('helpers.R')
load(file = 'data/probe_data/lfs_probes_0.05.rda')
# "array_first_with_cov_combat_keep_age_pc_keep_cancer_cancer_first_cancer_diagnosis_no_gender_control_cancerdraw_xgbTree.rda"
array_first = c(TRUE)
with_cov = c(TRUE)
plate_correction = c('pc_removal')
remove_cancer =c(FALSE)# FALSE (next one)
remove_cancer_first = c(FALSE)
read_file = 'data/train_test'
output_file = 'data/train_test'

for(i in array_first){
  for(j in with_cov){
    for(l in plate_correction){
      for(p in remove_cancer){
        for(q in remove_cancer_first){
          run_model(array_first = i, 
                    with_cov = j, 
                    plate_correction = l, 
                    remove_cancer = p, 
                    remove_cancer_first = q,
                    read_file = read_file, 
                    output_file = output_file)
        }
      }
    }
  }
}



