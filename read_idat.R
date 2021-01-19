# get helper function
library(minfi)

methyl_data_path <- 'data/idat_files/batch5/batch5/EPIC_Malkin&Goldenberg _SUB14749-51/'
clin_data_path <- 'data/clinical_data/clin.csv'
# read in clinical data
clin <- read.csv(clin_data_path)

idat_raw <- readRDS('data/idat_files/idat_raw.rda')

# read in idat files 
# idat_raw <- read.metharray.exp(methyl_data_path, recursive = T)

# use noob preprocessing 
Mset <- preprocessNoob(idat_raw, dyeMethod = 'single')

# dat <- getBeta(Mset)
dat <- wateRmelon::BMIQ(Mset)

saveRDS(dat, file = 'data/temp_noob.rda')

# use functional normilization
# Mset <- preprocessFunnorm(data)
# dat <- getBeta(Mset)

