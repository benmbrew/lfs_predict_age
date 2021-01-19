#!/hpf/tools/centos6/python/3.7.6_benbrew/bin/python3

#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('-mdepth', '--m_depth', type=int, help='max_depth_rf', default=2)
#parser.add_argument('-nest', '--n_est', type=float, help='column sample by tree', default=0.3)

#args = parser.parse_args()
#m_depth= args.m_depth
#n_est = args.n_est


# load libraries
import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
##############################
# ---- STEP 1: LOAD DATA ---- #
#dir_base = '/hpf/largeprojects/agoldenb/ben/Projects/nsqip/NSQIP_codes'
dir_base = os.getcwd()
dir_data =os.path.join(dir_base, 'data/train_test_big/')
dir_output =os.path.join(dir_base, 'data/train_test_big_model/')
valid_output =os.path.join(dir_base, 'data/train_valid_big_model/')

# get file names
train_file = 'full_Train_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) +'.csv'
test_file = 'test_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '.csv'

# read in data
dat_train = pd.read_csv(os.path.join(dir_data, train_file))
dat_test = pd.read_csv(os.path.join(dir_data, test_file))

# remove duplicate ids
dat_train = dat_train.sort_values(by=['ageofonset']).drop_duplicates('ids', keep='first')

# get test data clinical data
clin_test = dat_test[
    ['ids', 'Project', 'SentrixID', 'p53', 'tm_donor', 'tissue_type', 'cancer_diagnosis',
     'ageofonset', 'agesamplecollection', 'gender', 'array', 'cancer_atdraw', 'dataset',
     'age_label']]

# remove unneeded columns
dat_train.drop(['Unnamed: 0', 'ids', 'Project', 'SentrixID', 'p53', 'tm_donor', 'tissue_type',
                'cancer_diagnosis', 'ageofonset', 'agesamplecollection', 'array',
                'dataset'], axis=1, inplace=True)
dat_test.drop(['Unnamed: 0', 'ids', 'Project', 'SentrixID', 'p53', 'tm_donor', 'tissue_type',
               'cancer_diagnosis', 'ageofonset', 'agesamplecollection', 'array',
               'dataset'], axis=1, inplace=True)

# recode outcome
dat_train.replace('positive', 1, inplace=True)
dat_train.replace('negative', 0, inplace=True)
dat_test.replace('positive', 1, inplace=True)
dat_test.replace('negative', 0, inplace=True)

# get outcome variable
y_train = dat_train.age_label
y_test = dat_test.age_label

# remove outcome from training and testing data
dat_train.drop(['age_label'], inplace=True, axis=1)
dat_test.drop(['age_label'], inplace=True, axis=1)


if m == 'scale_data':
    # get index for numeric columns
    num_cols_train = dat_train.columns[
        dat_train.dtypes.apply(lambda c: np.issubdtype(c, np.number))]
    num_cols_test = dat_test.columns[
        dat_test.dtypes.apply(lambda c: np.issubdtype(c, np.number))]

    if model_name == 'rf' or model_name == 'xgb' or model_name == 'svm':
        min_max_scaler = preprocessing.MinMaxScaler()
    else:
        min_max_scaler = preprocessing.StandardScaler()

    dat_train[num_cols_train] = min_max_scaler.fit_transform(dat_train[num_cols_train])
    dat_test[num_cols_test] = min_max_scaler.fit_transform(dat_test[num_cols_test])

# condition for controlling cancer at draw
if n == 'no_draw':
    dat_train.drop(['cancer_atdraw'], inplace=True, axis=1)
    dat_test.drop(['cancer_atdraw'], inplace=True, axis=1)

# make gender a dummy variable in condition
if l == 'no_gender':
    dat_train.drop(['gender'], inplace=True, axis=1)
    dat_test.drop(['gender'], inplace=True, axis=1)


# create dummy variables for gender and cancer at draw
dat_train = pd.get_dummies(dat_train)
dat_test = pd.get_dummies(dat_test)

# get model parameters and run model

#################################################
# Random forest
if model_name == 'rf':
    param_file = 'best_params_rf_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))

    # get best parameters for random forest
    bootstrap = param_data.iloc[0]['0']
    max_depth =float(param_data.iloc[1]['0'])
    n_estimators = int(param_data.iloc[2]['0'])

    # Instantiate the grid search model
    clf = RandomForestClassifier(bootstrap=bootstrap, max_depth=max_depth, n_estimators=n_estimators)

    rf_mod = clf.fit(dat_train, y_train.values.ravel())
    rf_preds = rf_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = rf_preds.tolist()
    clin_test['real'] = y_test.tolist()

    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'rf_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished random forest ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n))

#################################################
# Xgboost
#'alpha': 0, 'colsample_bytree': 1, 'eta': 0.6, 'max_depth': 8, 'min_child_weight': 0, 'subsample': 0.8}
if model_name == 'xgb':
    param_file = 'best_params_xgb_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
            l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))

    # get best parameters
    alpha = param_data.iloc[0]['0']
    colsample_bytree = param_data.iloc[1]['0']
    eta = param_data.iloc[2]['0']
    max_depth = int(param_data.iloc[3]['0'])
    min_child_weight = param_data.iloc[4]['0']
    subsample = param_data.iloc[5]['0']

    clf = XGBClassifier(alpha= alpha, colsample_bytree= colsample_bytree, eta= eta, max_depth=max_depth,
                        min_child_weight=min_child_weight, subsample= subsample )

    xgb_mod = clf.fit(dat_train, y_train.values.ravel())
    xgb_preds = xgb_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = xgb_preds.tolist()
    clin_test['real'] = y_test.tolist()

    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'xgb_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                      l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished xgb ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
        r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

#################################################
# svm
if model_name == 'svm':
    param_file = 'best_params_svm_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))

    # get best parameters
    c = float(param_data.iloc[0]['0'])
    gamma = float(param_data.iloc[1]['0'])
    kernel = param_data.iloc[2]['0']

    clf = svm.SVC(probability=True, C=c, gamma=gamma,  kernel=kernel)

    svm_mod = clf.fit(dat_train, y_train.values.ravel())
    svm_preds = svm_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = svm_preds.tolist()
    clin_test['real'] = y_test.tolist()

    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'svm_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                      l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished svm ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
        r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))


#################################################
# logitstic
if model_name == 'logit':
    param_file = 'best_params_logit_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))

    # get best parameers
    c = float(param_data.iloc[0]['0'])
    penalty = param_data.iloc[1]['0']

    clf = LogisticRegression(max_iter=200, solver='liblinear', penalty=penalty, C=c)
    logit_mod = clf.fit(dat_train, y_train.values.ravel())
    logit_preds = logit_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = logit_preds.tolist()
    clin_test['real'] = y_test.tolist()


    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'logit_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                      l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished logit ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
        r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))


#################################################
# adaboost
if model_name == 'ada':
    param_file = 'best_params_ada_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))

    # get best parameter
    n_estimators = param_data.iloc[0]['0']
    clf = estimator=AdaBoostClassifier(n_estimators=n_estimators)
    ada_mod = clf.fit(dat_train, y_train.values.ravel())
    ada_preds = ada_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = ada_preds.tolist()
    clin_test['real'] = y_test.tolist()


    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'ada_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                      l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished ada ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
        r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))


#################################################
# mlp
if model_name == 'mlp':
    param_file = 'best_params_mlp_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '_' + str(
        l) + '_' + str(m) + '_' + str(n) + '.csv'

    # load param file
    param_data = pd.read_csv(os.path.join(valid_output, param_file))
    # get best parameters
    hidden_layer_sizes = param_data.iloc[2]['0']
    if hidden_layer_sizes == '(50, 50, 50)':
        hls = tuple([50,50,50])
    if hidden_layer_sizes == '(50, 100, 50)':
        hls = tuple([50,100,50])
    if hidden_layer_sizes == '(100,)':
        hls = tuple([100,])
    activation =param_data.iloc[0]['0']
    solver=param_data.iloc[4]['0']
    alpha =float(param_data.iloc[1]['0'])
    learning_rate =param_data.iloc[3]['0']

    clf=MLPClassifier(max_iter=500, hidden_layer_sizes=hls, activation=activation, solver=solver, alpha=alpha, learning_rate=learning_rate)

    mlp_mod = clf.fit(dat_train, y_train.values.ravel())
    mlp_preds = mlp_mod.predict_proba(dat_test)[:, 1]

    # join preds, test y, and test clin
    clin_test['preds'] = mlp_preds.tolist()
    clin_test['real'] = y_test.tolist()

    # save best model and clin test
    clin_test.to_csv(os.path.join(dir_output,
                                  'mlp_' + str(i) + '_' + str(j) + '_' + str(
                                      k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                      l) + '_' + str(
                                      m) + '_' + str(n) + '.csv'), index=False)

    print('finished mlp ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
        r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

