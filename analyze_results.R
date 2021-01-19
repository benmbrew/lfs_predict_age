source('helpers.R')
library(data.table)
library(ROCR)

# loop through validation results and stor plot in list
the_file_name <- 'data/train_valid_model/'

file_names <- list.files(the_file_name)
file_names <- file_names[!grepl('best', file_names)]
dat_list <- list()
for(i in 1:length(file_names)){
  this_file <- file_names[i]
  temp <- read.csv(paste0(the_file_name, this_file))
  temp$age_label <- factor(temp$age_label, levels=c('positive', 'negative'))
  print(nrow(temp))
  print(this_file)
  auc_value <- pROC::auc(temp$age_label, temp$preds)
  model_name <- ifelse(grepl('xgb', this_file), 'xgb', 
                       ifelse(grepl('rf', this_file), 'rf',
                              ifelse(grepl('svm', this_file), 'svm', 
                                     ifelse(grepl('logit', this_file), 'logit', 
                                            ifelse(grepl('enet', this_file), 'enet',
                                                   ifelse(grepl('ada', this_file), 'ada', 'mlp'))))))
  temp_results <- data_frame(key=this_file, value=auc_value, model_name=model_name)
  dat_list[[i]] <- temp_results
  print(i)
}

# combine results
temp <- do.call('rbind', dat_list)
temp <- temp[!is.na(temp$key),]

temp <- temp[order(temp$value, decreasing = TRUE),]
temp$value <- as.numeric(temp$value)
ggplot(temp, aes(x= reorder(key, -value), y = value, fill = model_name)) +geom_bar(stat='identity') +
  theme(axis.text.x= element_blank()) + 
  labs(x = 'models', y = 'auc') +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))


file_names <- temp$key
plot_list <-list()
for(i in 1:length(file_names)){
  this_file = file_names[i]
  temp <- read.csv(paste0(the_file_name, this_file))
  #
  temp$age_label <- factor(temp$age_label, levels=c('positive', 'negative'))
  auc_value <- round(pROC::auc(temp$age_label,temp$preds ), 3)
  temp_null <- temp[temp$cancer_diagnosis=='Unaffected',]
  temp_cancer <- temp[temp$cancer_diagnosis!='Unaffected',]
  
  null_plot=conmat_plot(data = temp_null,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title =paste0(this_file, ' auc=', auc_value ), data_type = 'null', text_name = 'age')
  cancer_plot=conmat_plot(data = temp_cancer,predict = 'preds', actual = 'age_label', cutoff = 0.125, get_plot = TRUE,other_title = paste0(this_file, ' auc=', auc_value ), data_type = 'cancer', text_name = 'age')
  
  plot_data <- list(null_plot, cancer_plot)
  plot_list[[i]] <- plot_data
  
}

# Another option: create pdf where each page is a separate plot.
pdf("plots_valid")
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()


temp$model_name <- unlist(lapply(strsplit(temp$key, '_'), function(x) x[1]))


temp <- temp %>% group_by(model_name) %>% slice_max(order_by = value, n =10)
