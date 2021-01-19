
library(tidyverse)
library(bapred)
library(sva)

get_combat_addon <- function(training_data, testing_data){
  pheno_names <- names(training_data)[!grepl('^cg', names(training_data))]
  train_pheno <- training_data[, pheno_names]
  test_pheno <- testing_data[, pheno_names]
  
  training_data$type <- 'train'
  testing_data$type <- 'test'
  temp_dat <- rbind(training_data, testing_data)
  temp_dat$Project <- as.numeric(temp_dat$Project)
  temp_dat$Project <- factor(temp_dat$Project, levels = c(1,2,3,4,5,6))
  train_data <- temp_dat %>% filter(type == 'train')
  test_data <- temp_dat %>% filter(type != 'train')
  train_data$type <- test_data$type <- NULL
  batch_train <- train_data$Project
  batch_test <- test_data$Project
  cg_names <- names(train_data)[grepl('^cg', names(train_data))]
  
  train_mat <- as.matrix(train_data[,cg_names])
  test_mat <- as.matrix(test_data[,cg_names])
  params <- combatba(x=train_mat, batch = batch_train)
  new_train <- as.data.frame(params$xadj)
  new_test <- as.data.frame(combatbaaddon(params, x = test_mat, batch = batch_test))
  train_out <- as.data.frame(cbind(train_pheno,new_train))
  test_out <- as.data.frame(cbind(test_pheno,new_test))
  
  return(list(train_out, test_out))
}
get_fsva <- function(training_data, testing_data){
  cg_names <- names(training_data)[grepl('^cg', names(training_data))]
  pheno_names <- names(training_data)[!grepl('^cg', names(training_data))]
  train_pheno <- training_data[, pheno_names]
  train_data <- as.matrix(t(training_data[,cg_names]))
  
  # test data
  test_pheno <- testing_data[, pheno_names]
  test_data <- as.matrix(t(testing_data[,cg_names]))
  
  # should 
  trainMod = model.matrix(~age_label,data=train_pheno)
  trainMod0 = model.matrix(~1,data=train_pheno)
  trainSv = sva(train_data,trainMod,trainMod0)
  
  fsvaobj = fsva(train_data,trainMod,trainSv,test_data)
  
  # get adjusted training and test (QUESTION WE USE ADJUSTED TRAINING?)
  new_train <- as.data.frame(t(fsvaobj[[1]]))
  new_test <- as.data.frame(t(fsvaobj[[2]]))
  train_out <- as.data.frame(cbind(train_pheno,new_train))
  test_out <- as.data.frame(cbind(test_pheno,new_test))
  
  return(list(train_out, test_out))
}
# registerDoParallel(6)
get_pc_location <- function(training_dat){
  pca_data <- as.data.frame(training_dat)
  
  # get other clinical data
  column_names <- names(pca_data)[!grepl('^cg', names(pca_data))]
  
  # get features sites
  cg_sites <- names(pca_data)[grepl('^cg', names(pca_data))]
  
  # run pca
  data_length <- ncol(pca_data)
  cg_start <- which(grepl('^cg', names(pca_data)))[1]
  pca <- prcomp(pca_data[,cg_start:data_length])
  
  # get pca dataframe with results and factor to color
  pca_results <- data.frame(pca$x[, 1:10],
                            pca_data[, column_names])
  
  # anova to get pc most correlated with project
  aov_pc1 <- summary(aov(PC1 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc2 <- summary(aov(PC2 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc3 <- summary(aov(PC3 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc4 <- summary(aov(PC4 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc5 <- summary(aov(PC5 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc6 <- summary(aov(PC6 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc7 <- summary(aov(PC7 ~ Project, data = pca_results))[[1]][5][1,]
  aov_pc8 <- summary(aov(PC8 ~ Project, data = pca_results))[[1]][5][1,]
  min_aov <- min(aov_pc1, aov_pc2, aov_pc3, aov_pc4, aov_pc5,aov_pc6, aov_pc7, aov_pc8)
  if(aov_pc1==min_aov){
    remove_pc <- 1
  } else if(aov_pc2==min_aov){
    remove_pc <- 2
  }else if(aov_pc3==min_aov){
    remove_pc <- 3
  }else if(aov_pc4==min_aov){
    remove_pc <- 4
  }else if(aov_pc5==min_aov){
    remove_pc <- 5
  }else if(aov_pc6==min_aov){
    remove_pc <- 6
  }else if(aov_pc7==min_aov){
    remove_pc <- 7
  }else if(aov_pc8==min_aov){
    remove_pc <- 8
  }
  return(remove_pc)
}
get_matching_replicates <- function(data) {
  cat("[ Getting technical replicates ]","\n")
  data <- data[!is.na(data$tm_donor),]
  data_450k <- data[data$array == "450",]
  data_850k <- data[data$array == "850",]
  tech_450k <- data_450k[data_450k$tm_donor %in% data_850k$tm_donor,]
  tech_450k <- tech_450k[!duplicated(tech_450k$tm_donor),]
  tech_ids <- apply(tech_450k,1, get_850replicate_indices,data_850k)
  tech_850k <- data_850k[data_850k$ids %in% tech_ids,]
  tech_450k <- tech_450k[tech_450k$tm_donor %in% tech_850k$tm_donor,]
  tech_450k <- tech_450k[match(tech_850k$tm_donor,tech_450k$tm_donor),]
  tech_all <- list(tech_450k,tech_850k)
  return(tech_all)
}


remove_pc <- function(temp_dat, pc_loc){
  # remove pc from training and test 
  clin_dat <- temp_dat[,!grepl('^cg', names(temp_dat))]
  feat_matrix <- temp_dat[,grepl('^cg', names(temp_dat))]
  mu = colMeans(feat_matrix)
  Xpca <- prcomp(feat_matrix)
  nComp = nrow(temp_dat)
  Xhat = Xpca$x[,-pc_loc] %*% t(Xpca$rotation[,-pc_loc])
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  
  final_dat <- as.data.frame(Xhat)
  final_dat <-as.data.frame(cbind(clin_dat, final_dat))
  return(final_dat)
  
}
get_850replicate_indices <- function(data_450k,data_850k) {
  age_450 <- as.numeric(data_450k["agesamplecollection"])
  replicate_850 <- data_850k[data_850k$tm_donor == data_450k["tm_donor"],]
  replicate_850$agediff <- replicate_850$agesamplecollection - age_450
  min <- replicate_850$ids[which.min(replicate_850$agediff)]
  return(min)
}
# remove batches with only on observation
remove_factor_level_rare <- function(temp_data){
  tt <- table(temp_data$Project)
  rare_levels <- names(tt)[tt<2]
  dd <- subset(temp_data,!Project %in% rare_levels)
  return(dd)
}


temp_450 = rep_450
temp_850 = rep_850
full_data = all_450
include_cov = TRUE
linear_transform <- function(temp_450, 
                             temp_850, 
                             full_data,
                             include_cov) {
  
  probe_model <- list()
  probe_450_result <- list()
  # get cg start 
  cg_start <- which(grepl('^cg', names(full_data)))[1]
  
  for (i in cg_start:ncol(full_data)) {
    if(include_cov){
      probe_450 <- as.data.frame(temp_450[, i])
      gen_450 <- as.data.frame(temp_450$gender)
      # cancer_450 <- as.data.frame(ifelse(temp_450$cancer_diagnosis=='Unaffected', 1, 0))
      age_450 <- as.data.frame(temp_450$agesamplecollection)
      
      probe_850<- as.data.frame(temp_850[, i])
      model_data <- data.frame(probe_850= probe_850, probe_450 = probe_450,age_450 = age_450,gen_450=gen_450)
      names(model_data) <- c('probe_850', 'probe_450', 'age_450', 'gen_450')
      probe_model[[i]] <- lm(probe_850 ~ probe_450+age_450+gen_450, data = model_data)
      
      probe_full <- as.data.frame(full_data[, i])
      gen_full <- as.data.frame(full_data$gender)
      # cancer_full <- as.data.frame(ifelse(full_data$cancer_diagnosis=='Unaffected', 1, 0))
      age_full <- as.data.frame(full_data$agesamplecollection)
      
      model_data_new <- data.frame(probe_full = probe_full, age_450 = age_full, 
                                   gen_450=gen_full)
      names(model_data_new) <- c('probe_450','age_450')
      probe_450_result[[i]] <- predict(probe_model[[i]], 
                                       newdata = model_data_new, 
                                       type = 'response')
      
      # print(i)
    } else {
      probe_450 <- as.data.frame(temp_450[, i])
      probe_850<- as.data.frame(temp_850[, i])
      model_data <- data.frame(probe_850= probe_850, probe_450 = probe_450)
      names(model_data) <- c('probe_850', 'probe_450')
      model_data$probe_850 <- as.numeric(as.character(model_data$probe_850))
      model_data$probe_450 <- as.numeric(as.character(model_data$probe_450))
      
      probe_model[[i]] <- lm(probe_850 ~ probe_450, data = model_data)
      probe_full <- as.numeric(full_data[, i])
      model_data_new <- data.frame(probe_full = probe_full)
      names(model_data_new) <- 'probe_450'
      probe_450_result[[i]] <- predict(probe_model[[i]], 
                                       newdata = model_data_new, 
                                       type = 'response')
      # print(i)
    }
    
  }
  
  # transpose results
  temp <- do.call(rbind, probe_450_result)
  transform_450 <- as.data.frame(t(temp))
  
  # add cg sites
  colnames(transform_450) <- colnames(full_data)[cg_start:ncol(full_data)]
  transform_450 <- as.data.frame(transform_450)
  
  # add clinical variables
  transform_450 <- as.data.frame(cbind(ids = full_data$ids, 
                                       Project = full_data$Project, 
                                       SentrixID = full_data$SentrixID, 
                                       p53 = full_data$p53, 
                                       tm_donor = full_data$tm_donor,
                                       tissue_type = full_data$tissue_type,
                                       cancer_diagnosis = full_data$cancer_diagnosis, 
                                       ageofonset = full_data$ageofonset, 
                                       agesamplecollection = full_data$agesamplecollection,
                                       gender = full_data$gender,
                                       array = full_data$array,
                                       cancer_atdraw=full_data$cancer_atdraw,
                                       dataset =full_data$dataset,
                                       transform_450))
  
  return(transform_450)
  
}

remove_pc_plate <- function(training_data, testing_data){
  
  named_cols <- names(training_data)
  testing_data <- testing_data[, named_cols]
  ind <- 14
  p_train <- prcomp(training_data[(ind+1):length(training_data)],scale=T)
  p_train_clin <- cbind(training_data[1:ind],p_train$x)
  
  # Find pcs correlated with age of sample collection
  pval <- list()
  for (i in 1:dim(p_train$x)[1]) {
    aov_pc <- aov(p_train$x[, i] ~ p_train_clin$Project)
    pval[[i]] <- summary(aov(p_train$x[, i] ~ p_train_clin$Project))[[1]][[1,"Pr(>F)"]]
  }
  
  pca_corr <- data.frame(pval=do.call('rbind', pval))
  pca_corr$rank <- rank(pca_corr$pval)
  pca_corr$pc <- paste0("PC",row.names(pca_corr))
  
  maxpc <- pca_corr$pc[pca_corr$rank == 1] ; maxpc_p <- pca_corr$pval[pca_corr$rank == 1]
  maxpc2 <- pca_corr$pc[pca_corr$rank == 2] ; maxpc_p2 <- pca_corr$pval[pca_corr$rank == 2]
  cat(paste0("[ Max age PC ] :",maxpc,"\t","[ Max age PC Association with Batch ] :",maxpc_p,"\n"))
  cat(paste0("[ Max age PC 2 ] :",maxpc2,"\t","[ Max age PC 2 Association with Batch ] :",maxpc_p2,"\n"))
  
  cat("[ Remove PC most correlated with Batch Effect ]",'\n')
  Xhat <- p_train$x[, !(colnames(p_train$x) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
  beta_adj <- scale(Xhat, center = -(colMeans(training_data[(ind+1):length(training_data)])), scale = T)
  new_train_data <- cbind(training_data[1:ind],beta_adj)

 
  ## Predict projection onto PC space in test data
  cat("[ Predict PCA in Test Data ]","\n")
  test_pred <-predict(p_train,testing_data[(ind+1):length(testing_data)])
  Xhat_pred <- test_pred[, !(colnames(test_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
  beta_adj_pred <- scale(Xhat_pred, center = -(colMeans(testing_data[(ind+1):length(testing_data)])), scale = T)
  new_test_data <- cbind(testing_data[1:ind],beta_adj_pred)

  return(list(new_train_data, new_test_data))


  
}


get_age_label <- function(temp_data, train_or_test){
  
  if(train_or_test == 'train'){
    t_null <- temp_data%>% filter(cancer_diagnosis=='Unaffected')
    t_cancer <- temp_data %>% filter(cancer_diagnosis!='Unaffected')
    # t_null <- t_null %>% filter(agesamplecollection >72)
    t_null$age_label <- 'negative'
    t_cancer$age_label <- ifelse(t_cancer$ageofonset <72, 'positive', 'negative')
    temp_data <- as.data.frame(rbind(t_cancer,t_null))
    cg_names <- names(temp_data)[grepl('^cg', names(temp_data))]
    pheno_names <-  names(temp_data)[!grepl('^cg', names(temp_data))]
    temp_data <- temp_data[,c(pheno_names,cg_names)]
  } else {
    t_null <- temp_data%>% filter(cancer_diagnosis=='Unaffected')
    t_cancer <- temp_data %>% filter(cancer_diagnosis!='Unaffected')
    # t_null <- t_null %>% filter(agesamplecollection >72)
    t_null$age_label <- 'negative'
    t_cancer$age_label <- ifelse(t_cancer$ageofonset <72, 'positive', 'negative')
    temp_data <- as.data.frame(rbind(t_cancer,t_null))
    cg_names <- names(temp_data)[grepl('^cg', names(temp_data))]
    pheno_names <-  names(temp_data)[!grepl('^cg', names(temp_data))]
    temp_data1 <- temp_data[,c(pheno_names,cg_names)]  
  }
  return(temp_data)
  
}
train_name <- '/hpf/largeprojects/agoldenb/ben/Projects/new_predict_age/new_data/train_25.rda' #full_train_25.rda
test_name <- '/hpf/largeprojects/agoldenb/ben/Projects/new_predict_age/new_data/valid_25.rda' # test_25.rda

train_name <- 'new_data/train_25.rda'
test_name <- 'new_data/valid_25.rda'


preprocess_data <- function(train_name, test_name, array_first, with_cov, plate_correction){# read in data
  # temp_train <- readRDS(file = '/hpf/largeprojects/agoldenb/ben/Projects/new_predict_age/train_data.rda')
  temp_train <- readRDS(file = train_name)
  temp_test <- readRDS(file = test_name)
 

  if(array_first){
    # combine train and valid/test for array correction
    temp_train$dataset <- 'train'
    temp_test$dataset <- 'test'
    temp_train <- rbind(temp_train, temp_test)
    cg_names <- names(temp_train)[grepl('^cg', names(temp_train))]
    clin_names <- names(temp_train)[!grepl('^cg', names(temp_train))]
    temp_train <- temp_train[c(clin_names,cg_names)]
    rep_data <- get_matching_replicates(temp_train)
    rep_450 <- rep_data[[1]]
    rep_450_tm_donor <- unique(rep_450$tm_donor)
    rep_850 <- rep_data[[2]]
    all_450 <- temp_train[temp_train$array == '450',]
    all_450 <- all_450[!all_450$tm_donor %in% rep_450_tm_donor,]
    all_850 <- temp_train[temp_train$array =='850',]
    # run function to get new all_450 and then rbind with 850
    new_450 <- linear_transform(temp_450 = rep_450, temp_850 = rep_850, full_data = all_450, include_cov = with_cov)
    # save.image('temp_array_first_w_cov.RData')
    transformed_train <- rbind(all_850, new_450)
    temp_test <- transformed_train[transformed_train$dataset =='test',]
    transformed_train <- transformed_train[transformed_train$dataset =='train',]
    
    transformed_train <- get_age_label(transformed_train, train_or_test = 'train')
    temp_test <- get_age_label(temp_test, train_or_test = 'test')
    # save(transformed_train, temp_test, file ='temp_array_first_w_cov.rda')
    if(plate_correction == 'fsva'){
      # get fsva data
      fsva_results <- get_fsva(training_data = transformed_train, 
                               testing_data = temp_test)
      new_train <- fsva_results[[1]]
      new_test <- fsva_results[[2]]
      # save(new_train, new_test, file ='temp_array_first_w_cov_fsva.rda')
      
    } else if(plate_correction=='combat'){
      # do combat addon
      # transformed_train <- remove_factor_level_rare(transformed_train)
      # temp_test <- remove_factor_level_rare(temp_test)
      # 
      combat_results <- get_combat_addon(training_data = transformed_train,
                                         testing_data = temp_test)
      new_train <- combat_results[[1]]
      new_test <- combat_results[[2]]

    } else if(plate_correction=='pc_removal') {
      adj_pc_data <- remove_pc_plate(training_dat = transformed_train, testing_dat = temp_test)
      new_train <- adj_pc_data[[1]]
      new_test <- adj_pc_data[[2]]
    }
    
  } else {
    temp_train$dataset <- 'train'
    temp_test$dataset <- 'test'
    temp_train <- rbind(temp_train, temp_test)
    cg_names <- names(temp_train)[grepl('^cg', names(temp_train))]
    clin_names <- names(temp_train)[!grepl('^cg', names(temp_train))]
    temp_train <- temp_train[c(clin_names,cg_names)]
    
    rep_data <- get_matching_replicates(temp_train)
    rep_450 <- rep_data[[1]]
    rep_450_tm_donor <- unique(rep_450$tm_donor)
    rep_850 <- rep_data[[2]]
    all_450 <- temp_train[temp_train$array == '450',]
    all_450 <- all_450[!all_450$tm_donor %in% rep_450_tm_donor,]
    all_850 <- temp_train[temp_train$array =='850',]
    # run function to get new all_450 and then rbind with 850
    temp_train <- rbind(all_850, all_450)
    temp_test <- temp_train[temp_train$dataset == 'test',]
    temp_train <- temp_train[temp_train$dataset == 'train',]
    
    temp_train <- get_age_label(temp_train, train_or_test = 'train')
    temp_test <- get_age_label(temp_test, train_or_test = 'test')
    
    
    if(plate_correction == 'fsva'){
      # get fsva data
      fsva_results <- get_fsva(training_data = temp_train, 
                               testing_data = temp_test)
      new_train <- fsva_results[[1]]
      new_test <- fsva_results[[2]]
      
    } else if(plate_correction=='combat'){
      remove_factor_level_rare <- function(temp_data){
        tt <- table(temp_data$Project)
        rare_levels <- names(tt)[tt<2]
        dd <- subset(temp_data,!Project %in% rare_levels)
        return(dd)
      }
      
      temp_train <- remove_factor_level_rare(temp_train)
      temp_test <- remove_factor_level_rare(temp_test)
      # do combat addon
      combat_results <- get_combat_addon(training_data = temp_train,
                                         testing_data = temp_test)
      new_train <- combat_results[[1]]
      new_test <- combat_results[[2]]
      
    } else if(plate_correction=='pc_removal') {
      adj_pc_data <- remove_pc_plate(training_dat = temp_train, testing_dat = temp_test)
      new_train <- adj_pc_data[[1]]
      new_test <- adj_pc_data[[2]]
    }
   
    all_450 <- new_train[new_train$array == '450',]
    all_850 <- new_train[new_train$array =='850',]
    all_450$age_label <- NULL
    all_850$age_label <- NULL
    # run function to get new all_450 and then rbind with 850
    new_450 <- linear_transform(temp_450 = rep_450, temp_850 = rep_850, full_data = all_450, include_cov = with_cov)
    new_train <- rbind(all_850, new_450)
    new_train <- get_age_label(new_train, train_or_test = 'train')
    
  }
  return(list(new_train, new_test))
}

# create parameters to loop through
array_first <- c(FALSE)
with_cov <- c(FALSE)
plate_correction <- c('combat')

for(i in 1:length(array_first)){
  array_param <- array_first[[i]]
  if(array_param){
    array_save <- 'array_first'
  } else {
    array_save <- 'plate_first'
  }
  for(j in 1:length(with_cov)){
    with_cov_param <- with_cov[[j]]
    if(with_cov_param){
      with_cov_save <- 'with_cov'
    } else {
      with_cov_save <- 'no_cov'
    }
    for(k in 1:length(plate_correction)){
      plate_param <- plate_correction[[k]]
     
     temp <-  preprocess_data(train_name = train_name, test_name = test_name, array_first = array_param, with_cov = with_cov_param, plate_correction = plate_param)
      message(array_save,' ', with_cov_save, ' ',plate_param)
      saveRDS(temp, file = paste0('/hpf/largeprojects/agoldenb/ben/Projects/new_predict_age/train_valid/', array_save, '_', with_cov_save, '_', plate_param, '.rda'))
    }
  } 
}

# 
