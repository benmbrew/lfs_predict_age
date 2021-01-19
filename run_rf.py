#!/hpf/tools/centos6/python/3.7.6_benbrew/bin/python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-mdepth', '--m_depth', type=int, help='max_depth_rf', default=2)
parser.add_argument('-nest', '--n_est', type=float, help='column sample by tree', default=0.3)

args = parser.parse_args()
m_depth= args.m_depth
n_est = args.n_est
n_est = 50
m_depth = 4

import numpy as np
import pandas as pd
import os
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

###############################
# ---- STEP 1: LOAD DATA ---- #
#dir_base = '/hpf/largeprojects/agoldenb/ben/Projects/nsqip/NSQIP_codes'
dir_base = os.getcwd()
dir_data =os.path.join(dir_base, 'data')
fn_X = 'X_imputed.csv'
fn_Y = 'y_agg.csv'
dat_X = pd.read_csv(os.path.join(dir_data, fn_X))
dat_Y = pd.read_csv(os.path.join(dir_data, fn_Y))