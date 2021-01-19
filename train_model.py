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
dir_data =os.path.join(dir_base, 'data/train_valid_big/')
dir_output =os.path.join(dir_base, 'data/train_valid_big_model/')

#svm_plate_first_no_cov_pc_removal_control_gender_no_scale.csv
# set parameters to loop through
array_first = ['array_first']
with_cov = ['with_cov', 'no_cov']
plate_correction = ['pc_removal', 'combat']
remove_cancer = ['keep_cancer']
remove_cancer_first = ['lfs_first']
control_cancer = ['control_draw', 'no_draw']
control_gender = ['control_gender', 'no_gender']
scale_data = ['scale_data', 'no_scale']


for i in array_first:
    for j in with_cov:
        for k in plate_correction:
            for r in remove_cancer:
                for s in remove_cancer_first:
                    for l in control_gender:
                        for n in control_cancer:
                            for m in scale_data:
                                # get file names
                                train_file = 'train_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) +'.csv'
                                test_file = 'valid_' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(r) + '_' + str(s) + '.csv'

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

                                    #min_max_scaler = preprocessing.MinMaxScaler()
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

                                # ##################
                                # Random forest
                                param_grid = {
                                    'bootstrap': [True],
                                    'max_depth': [80,  110],
                                    'n_estimators': [50, 200]
                                }
                                # Instantiate the grid search model
                                clf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=param_grid, cv=5, n_jobs=8,
                                                   verbose=1, scoring='roc_auc')
                                rf_mod = clf.fit(dat_train, y_train.values.ravel())
                                rf_preds = rf_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = rf_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = rf_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'rf_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' +str(r) + '_' + str(s) + '_' + str(l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                              'best_params_rf_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' +str(r) + '_' + str(s) + '_' + str(l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)


                                print('finished random forest '+str(i) + '_' + str(j) + '_' + str(k)+'_'  + str(r) +'_' +  str(s) +'_' + str(l) + '_' + str(m) + '_' + str(n))

                                #################
                                #Extreme gradient boosting
                                param_grid = {
                                    'alpha': [0],
                                    'max_depth': [8],
                                    'min_child_weight': [0],
                                    'eta': [0.3],
                                    'subsample': [1],
                                    'colsample_bytree': [1]
                                }

                                clf = GridSearchCV(estimator=XGBClassifier(), param_grid=param_grid, cv=5, n_jobs=8, verbose=1, scoring='roc_auc')
                                xgb_mod = clf.fit(dat_train, y_train.values.ravel())
                                xgb_preds = xgb_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = xgb_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = xgb_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'xgb_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                  l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                                'best_params_xgb_' + str(i) + '_' + str(j) + '_' + str(
                                                                    k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                    l) + '_' + str(
                                                                    m) + '_' + str(n) + '.csv'), index=False)

                                print('finished xgb ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                    r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

                                ##################
                                # Support vector machines
                                param_grid = {'C': [1, 5],
                                             'gamma': [0.01, 0.001],
                                             'kernel': ['rbf']}
                                #
                                clf = GridSearchCV(estimator=svm.SVC(probability=True), param_grid=param_grid, cv=5, n_jobs=8,
                                                    verbose=1, scoring='roc_auc')
                                #clf = svm.SVC(probability=True, C=10, gamma=0.0001)

                                svm_mod = clf.fit(dat_train, y_train.values.ravel())
                                svm_preds = svm_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = svm_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = svm_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'svm_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                  l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                                'best_params_svm_' + str(i) + '_' + str(j) + '_' + str(
                                                                    k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                    l) + '_' + str(
                                                                    m) + '_' + str(n) + '.csv'), index=False)

                                print('finished svm ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                    r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

                                ##################
                                # Elastic net
                                # param_grid = {"penalty": ['elasticnet'],
                                #               "l1_ratio": np.arange(0.0, 1.0, 0.1)
                                #               }
                                # clf = GridSearchCV(estimator=LogisticRegression(solver='saga', max_iter=500), param_grid=param_grid,
                                #                    cv=5, n_jobs=8, verbose=1)
                                # enet_mod = clf.fit(dat_train, y_train.values.ravel())
                                # enet_preds = enet_mod.predict_proba(dat_test)[:, 1]
                                #
                                # # join preds, test y, and test clin
                                # clin_test['preds'] = enet_preds.tolist()
                                # clin_test['real'] = y_test.tolist()
                                #
                                # # get best parameters
                                # best_params = enet_mod.best_params_
                                # best_params = pd.DataFrame.from_dict(best_params, orient="index")
                                #
                                # # save best model and clin test
                                # clin_test.to_csv(os.path.join(dir_output,
                                #                               'enet_' + str(i) + '_' + str(j) + '_' + str(
                                #                                   k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                #                                   l) + '_' + str(
                                #                                   m) + '_' + str(n) + '.csv'), index=False)
                                # best_params.to_csv(os.path.join(dir_output,
                                #                                 'best_params_enet_' + str(i) + '_' + str(j) + '_' + str(
                                #                                     k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                #                                     l) + '_' + str(
                                #                                     m) + '_' + str(n) + '.csv'), index=False)
                                #
                                # print('finished enet ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                #     r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

                                ##################
                                # Logistic (lasso, ridge)
                                param_grid = {"penalty": ['l2', 'l1'],
                                              "C": np.logspace(-3, 3, 7)
                                              }
                                clf = GridSearchCV(estimator=LogisticRegression(max_iter=200, solver='liblinear'), param_grid=param_grid, cv=5,
                                                   n_jobs=8, verbose=1)
                                logit_mod = clf.fit(dat_train, y_train.values.ravel())
                                logit_preds = logit_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = logit_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = logit_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'logit_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                  l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                                'best_params_logit_' + str(i) + '_' + str(j) + '_' + str(
                                                                    k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                    l) + '_' + str(
                                                                    m) + '_' + str(n) + '.csv'), index=False)

                                print('finished logit ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                    r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

                                ##################
                                # adaboost
                                param_grid = {"n_estimators": [30, 50, 100, 200]}
                                clf = GridSearchCV(estimator=AdaBoostClassifier(), param_grid=param_grid,
                                                   cv=5,
                                                   n_jobs=8, verbose=1)
                                ada_mod = clf.fit(dat_train, y_train.values.ravel())
                                ada_preds = ada_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = ada_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = ada_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'ada_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                  l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                                'best_params_ada_' + str(i) + '_' + str(
                                                                    j) + '_' + str(
                                                                    k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                    l) + '_' + str(
                                                                    m) + '_' + str(n) + '.csv'), index=False)

                                print('finished ada ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                    r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))

                                ##################
                                # mlp classfier
                                param_grid = {
                                    'hidden_layer_sizes': [(50, 50, 50), (50, 100, 50), (100,)],
                                    'activation': ['tanh', 'relu'],
                                    'solver': ['sgd', 'adam'],
                                    'alpha': [0.0001, 0.05],
                                    'learning_rate': ['constant', 'adaptive'],
                                }
                                clf = GridSearchCV(estimator=MLPClassifier(max_iter=500), param_grid=param_grid,
                                                   cv=5,
                                                   n_jobs=8, verbose=1)
                                mlp_mod = clf.fit(dat_train, y_train.values.ravel())
                                mlp_preds = mlp_mod.predict_proba(dat_test)[:, 1]

                                # join preds, test y, and test clin
                                clin_test['preds'] = mlp_preds.tolist()
                                clin_test['real'] = y_test.tolist()

                                # get best parameters
                                best_params = mlp_mod.best_params_
                                best_params = pd.DataFrame.from_dict(best_params, orient="index")

                                # save best model and clin test
                                clin_test.to_csv(os.path.join(dir_output,
                                                              'mlp_' + str(i) + '_' + str(j) + '_' + str(
                                                                  k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                  l) + '_' + str(
                                                                  m) + '_' + str(n) + '.csv'), index=False)
                                best_params.to_csv(os.path.join(dir_output,
                                                                'best_params_mlp_' + str(i) + '_' + str(
                                                                    j) + '_' + str(
                                                                    k) + '_' + str(r) + '_' + str(s) + '_' + str(
                                                                    l) + '_' + str(
                                                                    m) + '_' + str(n) + '.csv'), index=False)

                                print('finished mlp ' + str(i) + '_' + str(j) + '_' + str(k) + '_' + str(
                                    r) + '_' + str(s) + '_' + str(l) + '_' + str(m) + '_' + str(n))


