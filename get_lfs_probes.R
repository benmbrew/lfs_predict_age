# load libraries
library(sva)
library(limma)
library(caret)
library(randomForest)
library(tidyverse)
library(readr)
library(ggrepel)
library(bapred)
library(multcomp)
library(bumphunter)
library(caret)
library(pROC)
source('../predict_age/all_functions.R')

# 6056-6147
# LFS in tmdonor from NIH 
# all sickkids are all 2,3,4 digits. 5 and 6 
# use array type
# read in noob data with no corrections
dat <- readRDS('../../Data/noob/Noob_beta.rds')
dat <- dat[!is.na(dat$ids),]
dat <- dat[!dat$Meth %in% c("Problem"),]
dat <- dat[!dat$SentrixID %in% c("203836210070_R08C01", "203836210068_R02C01", "203836210104_R01C01","5760666022_R06C02","203836210087_R07C01"),]

# read in recent data with no outliers and subset data
train <- readRDS('../../Data/noob/NoobCorrected_beta2_baPredComBat_TrainingSet.rds')
test <- readRDS('../../Data/noob/NoobCorrected_beta2_baPredComBat_TestSet.rds')
tm_donors <- union(test$tm_donor, train$tm_donor)

# remove outliers 
dat <- dat %>% filter(tm_donor %in% tm_donors)

# get 450 data
dat_450 <- dat %>% filter(array == '450')

method = 'noob'
cases_wt_450 <- readRDS(paste0('../../Data/', method,'/cases_wt_450_beta.rda'))
con_wt_450 <- readRDS(paste0('../../Data/', method,'/controls_wt_450_beta.rda'))
rm(test,train, dat)
# combine cases, controls, valid and pca

cases_wt_450 <- clean_dat(cases_wt_450, tech = '450k',cases_or_controls = 'cases', mut_or_wt = 'wt')
con_wt_450 <- clean_dat(con_wt_450, tech = '450k',cases_or_controls = 'con', mut_or_wt = 'wt')
con_wt_450 <- con_wt_450[!grepl('3077', con_wt_450$tm_donor),]

wt_450 <- rbind(cases_wt_450, con_wt_450)
rm(cases_wt_450, con_wt_450)

# add identifier for data
wt_450$data <- 'old_data'
dat_450$data <- 'new_data'

# homegenize column names
names(dat_450)[4] <- 'p53_germline'
wt_cgs <- names(wt_450)[grepl('^cg', names(wt_450))]
mut_cgs <- names(dat_450)[grepl('^cg', names(dat_450))]
int_cgs <- intersect(wt_cgs, mut_cgs)

dat_450 <- dat_450[, c('p53_germline','data', int_cgs)]
wt_450 <- wt_450[, c('p53_germline', 'data',int_cgs)]

# get pca
g_ranges <- readRDS('../../Data/g_ranges.rda')
g_ranges$probe <- rownames(g_ranges)
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]
names(g_ranges)[1] <- 'chr'

# get pca
get_pca(pca_data = rbind(dat_450, wt_450), column_name = 'data', show_variance = FALSE, pc_x = 1, pc_y = 2, main_title = '')

names(dat_450)[2] <- 'tech'
names(wt_450)[2] <- 'tech'
# combat to correct for data
temp <- run_combat(temp_data = rbind(dat_450, wt_450), type ='tech' )

get_pca(pca_data = temp, column_name = 'tech', show_variance = FALSE, pc_x = 1, pc_y = 2, main_title = '')

# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = wt_450, 
                        dat_2 = dat_450, 
                        bump = 'lfs', 
                        boot_num = 50, 
                        beta_thresh = 0.05,
                        methyl_type = 'beta',
                        g_ranges = g_ranges)

colnames(bh_feats)[1] <- 'chr'
bh_feats1 <- bh_feats %>% filter(p.value<=0.05)
keep_features <- inner_join(bh_feats, g_ranges)$probe
save(keep_features, file='lfs_probes_0.05.rda')

# cases
cases_450_small <- join_new_features(cases_450, new_features = bh_feats)

# lfs probes 
lfs_bump_probes <- colnames(cases_450_small)[grepl('^cg', colnames(cases_450_small))]
