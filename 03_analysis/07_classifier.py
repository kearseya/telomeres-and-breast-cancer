#!/bin/python3

import pandas as pd
import numpy as np
import statistics
from itertools import chain

## Inherently multiclass
from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import ExtraTreeClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC #(setting multi_class=”crammer_singer”)
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.linear_model import LogisticRegression #(setting multi_class=”multinomial”)
from sklearn.linear_model import LogisticRegressionCV #(setting multi_class=”multinomial”)
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neighbors import RadiusNeighborsClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.semi_supervised import LabelPropagation
from sklearn.semi_supervised import LabelSpreading
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import AdaBoostClassifier

#Work well
#[LogisticRegression, MLPClassifier, RidgeClassifier, GaussianNB]

from sklearn.datasets import make_classification
from sklearn.impute import SimpleImputer
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import uniform
from sklearn.model_selection import cross_validate
from sklearn. model_selection import cross_val_score
from sklearn.model_selection import LeaveOneOut

train = pd.read_csv("/path/to/project/data/train_clusters.csv")#, index_col = "Clinical index")
other = pd.read_csv("/path/to/project/data/BC_clinical_and_TL_data.csv")#, index_col = "Clinical index")

train_tumour_db = train["tumor_db"].tolist()
other = other[~other["Tumor_DB_id"].isin(train_tumour_db)]
other = other.rename(columns={"Mean": "tumor_stela"})

variables = ["tumor_stela", "AGE AT OPERATION", "NPI SCORE"]

train = train[list(chain(*[variables, ["cl_members", "Clinical index", "chromothripsis_p0.04_n8", "num_total_mut", "Signature.5", "APOBEC.bin", "diff"]]))]
train.loc[(train["chromothripsis_p0.04_n8"].isnull()) & (train["tumor_stela"] <= 3.81), "chromothripsis_p0.04_n8"] = 1
train.loc[(train["chromothripsis_p0.04_n8"].isnull()) & (train["tumor_stela"] > 3.81), "chromothripsis_p0.04_n8"] = 0

other = other[list(chain(*[variables, ["Clinical index"]]))]

impute = True
print(train[train[variables].isna().any(axis=1)][variables])
if impute == True:
    #Tumor stela can be imputed with age+npi as mainly ttl missing, only one NPI missing
    #imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp = IterativeImputer()
    imp_train = imp.fit(train[variables])
    imp_other = imp.fit(other[variables])
    train[variables] = imp_train.transform(train[variables])
    train["diff"] = imp_train.transform(train[variables])
    other[variables] = imp_other.transform(other[variables])
    other = other.reset_index()

if impute == False:
    train = train.dropna(subset=variables)
    other = other.dropna(subset=variables)

print(train)
print(other)
# view whole imputed file
# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#     print(other)




def ensemble_manual_test_new_var(train, other, method_list, variable):
    method_dict = {
    "LogisticRegression": LogisticRegression(max_iter=1000000),
    "MLPClassifier": MLPClassifier(max_iter=1000000),
    "RidgeClassifier": RidgeClassifier(),
    "GaussianNB": GaussianNB(),
    "RandomForestClassifier": RandomForestClassifier()
    }

    def fitting(train, other, include_new, method, variable, n_cl):
        if n_cl == 2:
            train = train.copy()
            train["cl_members"] = train.cl_members.str[:1]
        clf = method_dict[method]
        loo = LeaveOneOut()
        if include_new == False:
            clf.fit(train[variables], train["cl_members"])
            scores = cross_val_score(clf, train[variables], train["cl_members"], cv=loo.split(train))
        if include_new == True:
            clf.fit(train[list(chain(*[variables, [variable]]))], train["cl_members"])
            scores = cross_val_score(clf, train[list(chain(*[variables, [variable]]))], train["cl_members"], cv=loo.split(train))
        return(scores)

    include_new = False
    nochromo_2 = []
    nochromo_3 = []
    for i, method in enumerate(method_list):
        print("without: ", i+1, "/", len(method_list), " ", method, end='\r')
        nochromo_2.append(fitting(train, other, include_new, method, variable, 2))
        nochromo_3.append(fitting(train, other, include_new, method, variable, 3))
    print("")
    include_new = True
    wchromo_2 = []
    wchromo_3 = []
    for i, method in enumerate(method_list):
        print("with:    ", i+1, "/", len(method_list), " ", method, end='\r')
        wchromo_2.append(fitting(train, other, include_new, method, variable, 2))
        wchromo_3.append(fitting(train, other, include_new, method, variable, 3))
    print("")
    print("")

    print("Without", variable, ": ")
    print("2 clusters: ", np.mean(nochromo_2))
    print(np.mean(nochromo_2, axis=1))
    print("3 clusters: ", np.mean(nochromo_3))
    print(np.mean(nochromo_3, axis=1))
    print("")
    print("With", variable,":")
    print("2 clusters: ", np.mean(wchromo_2))
    print(np.mean(wchromo_2, axis=1))
    print("3 clusters: ", np.mean(wchromo_3))
    print(np.mean(wchromo_3, axis=1))


#ensemble_manual_test_new_var(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], "chromothripsis_p0.04_n8")
#ensemble_manual_test_new_var(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], "diff")
#ensemble_manual_test_new_var(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], "Signature.5")
#ensemble_manual_test_new_var(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], "num_total_mut")




def model_manual_test_mode_new_var(train, other, iterations, method, variable):
    method_dict = {
    "LogisticRegression": LogisticRegression(max_iter=1000000),
    "MLPClassifier": MLPClassifier(max_iter=1000000),
    "RidgeClassifier": RidgeClassifier(),
    "GaussianNB": GaussianNB(),
    "RandomForestClassifier": RandomForestClassifier()
    }

    def fitting(train, other, include_new, method, variable, n_cl):
        if n_cl == 2:
            train = train.copy()
            train["cl_members"] = train.cl_members.str[:1]
        clf = method_dict[method]
        loo = LeaveOneOut()
        if include_new == False:
            clf.fit(train[variables], train["cl_members"])
            scores = cross_val_score(clf, train[variables], train["cl_members"], cv=loo.split(train))
        if include_new == True:
            clf.fit(train[list(chain(*[variables, [variable]]))], train["cl_members"])
            scores = cross_val_score(clf, train[list(chain(*[variables, [variable]]))], train["cl_members"], cv=loo.split(train))
        return(scores)

    include_new = False
    nochromo_2 = []
    nochromo_3 = []
    for i in range(iterations):
        print("without: ", i+1, "/", iterations, " ", method, end='\r')
        nochromo_2.append(fitting(train, other, include_new, method, variable, 2))
        nochromo_3.append(fitting(train, other, include_new, method, variable, 3))
    print("")
    include_new = True
    wchromo_2 = []
    wchromo_3 = []
    for i in range(iterations):
        print("with:    ", i+1, "/", iterations, " ", method, end='\r')
        wchromo_2.append(fitting(train, other, include_new, method, variable, 2))
        wchromo_3.append(fitting(train, other, include_new, method, variable, 3))
    print("")
    print("")
    print("Without", variable, ": ")
    print("2 clusters: ", np.mean(nochromo_2))
    print(np.mean(nochromo_2, axis=1))
    print("3 clusters: ", np.mean(nochromo_3))
    print(np.mean(nochromo_3, axis=1))
    print("")
    print("With", variable,":")
    print("2 clusters: ", np.mean(wchromo_2))
    print(np.mean(wchromo_2, axis=1))
    print("3 clusters: ", np.mean(wchromo_3))
    print(np.mean(wchromo_3, axis=1))



# model_manual_test_mode_new_var(train, other, 10, "RandomForestClassifier", "chromothripsis_p0.04_n8")
# model_manual_test_mode_new_var(train, other, 10, "MLPClassifier", "chromothripsis_p0.04_n8")
#
# model_manual_test_mode_new_var(train, other, 10, "RandomForestClassifier", "diff")
# model_manual_test_mode_new_var(train, other, 10, "MLPClassifier", "diff")
#
# model_manual_test_mode_new_var(train, other, 10, "RandomForestClassifier", "num_total_mut")
# model_manual_test_mode_new_var(train, other, 10, "MLPClassifier", "num_total_mut")
#
# model_manual_test_mode_new_var(train, other, 10, "RandomForestClassifier", "Signature.5")
# model_manual_test_mode_new_var(train, other, 10, "MLPClassifier", "Signature.5")
#
# model_manual_test_mode_new_var(train, other, 10, "RandomForestClassifier", "APOBEC.bin")
# model_manual_test_mode_new_var(train, other, 10, "MLPClassifier", "APOBEC.bin")



def manual_mode_pred(train, other, method_list, iterations, n_cl):
    pred_table_mode = other.copy()
    raw_scores = {}
    avg_scores = []
    ind_pred_name = []

    method_dict = {
    "LogisticRegression": LogisticRegression(max_iter=1000000),
    "MLPClassifier": MLPClassifier(max_iter=1000000),
    "RidgeClassifier": RidgeClassifier(),
    "GaussianNB": GaussianNB(),
    "RandomForestClassifier": RandomForestClassifier()
    }

    def fitting(train, other, method, n_cl):
        pred_cl = other.copy()
        if n_cl == 2:
            train = train.copy()
            train["cl_members"] = train.cl_members.str[:1]
        clf = method_dict[method]
        loo = LeaveOneOut()
        clf.fit(train[variables], train["cl_members"])
        pred_cl["pred_cl"] = clf.predict(other[variables])
        scores = cross_val_score(clf, train[variables], train["cl_members"], cv=loo.split(train))
        return(pred_cl, scores)

    for i, method in enumerate(method_list):
        if method in ["MLPClassifier", "RandomForestClassifier"]:
            tmp_pred_table_mode = other.copy()
            tmp_raw_scores = []
            tmp_ind_pred_name = []
            for j in range(iterations):
                print("Fitting: ", i+1,"/",len(method_list),"  (", j+1, "/", iterations, "  ", method,")", end='\r')
                clf = fitting(train, other, method, n_cl)
                tmp_name = "_".join([str(method), str(j)])
                tmp_ind_pred_name.append(tmp_name)
                tmp_pred_table_mode[tmp_name] = clf[0]["pred_cl"]
                tmp_raw_scores.append(clf[1])
                if j == (iterations-1):
                    name = "_".join(["pred", str(method)])
                    ind_pred_name.append(name)
                    pred_table_mode[name] = tmp_pred_table_mode[tmp_ind_pred_name].mode(axis='columns')[0]
                    raw_scores[name] = np.mean(tmp_raw_scores)
        else:
            print("Fitting: ", i+1, "/", len(method_list), " ", method, end='\r')
            clf = fitting(train, other, method, n_cl)
            name = "_".join(["pred", str(method)])
            ind_pred_name.append(name)
            pred_table_mode[name] = clf[0]["pred_cl"]
            raw_scores[name] = clf[1]

    print("")
    pred_table_mode['mode_pred'] = pred_table_mode[ind_pred_name].mode(axis='columns')[0]
    for s in raw_scores:
        avg_scores.append(np.mean(raw_scores[s]))
        pred_table_mode = pred_table_mode.rename(columns={s: s+"_"+str(round(np.mean(raw_scores[s]), 3))})
    print(avg_scores)
    print("Mean score: ", np.mean(avg_scores))
    print(pred_table_mode)
    return(pred_table_mode)


#mode_pred_3 = manual_mode_pred(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], iterations = 10, n_cl = 3)
#mode_pred_3.to_csv("/path/to/project/data/clf_mode_pred_cl_n3.csv")

#mode_pred_2 = manual_mode_pred(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], iterations = 10, n_cl = 2)
#mode_pred_2.to_csv("/path/to/project/data/clf_mode_pred_cl_n2_.csv")






## use if multiple iterations not needed
def ensemble_manual_predict(train, other, method_list, n_cl=3):
    if n_cl == 2:
        train = train.copy()
        train["cl_members"] = train.cl_members.str[:1]
    method_dict ={
    "LogisticRegression": LogisticRegression(max_iter=1000000),
    "MLPClassifier": MLPClassifier(max_iter=1000000),
    "RidgeClassifier": RidgeClassifier(),
    "GaussianNB": GaussianNB(),
    "RandomForestClassifier": RandomForestClassifier()
    }
    pred_table_mode = other.copy()
    raw_scores = {}
    avg_scores = []
    ind_pred_name = []
    def fitting(train, other, method):
        clf = method_dict[method]
        loo = LeaveOneOut()
        clf.fit(train[variables], train["cl_members"])
        scores = cross_val_score(clf, train[variables], train["cl_members"], cv=loo.split(train))
        pred_cl = other.copy()
        pred_cl["pred_cl"] = clf.predict(other[variables])
        return(pred_cl, scores)
    for i, method in enumerate(method_list):
        print("Fitting: ", i+1, "/", len(method_list), " ", method, end='\r')
        clf = fitting(train, other, method)
        name = "_".join(["pred", str(method)])
        ind_pred_name.append(name)
        pred_table_mode[name] = clf[0]["pred_cl"]
        raw_scores[name] = clf[1]
    pred_table_mode['mode_pred'] = pred_table_mode[ind_pred_name].mode(axis='columns')[0]
    for s in raw_scores:
        avg_scores.append(np.mean(raw_scores[s]))
        pred_table_mode = pred_table_mode.rename(columns={s: s+"_"+str(round(np.mean(raw_scores[s]), 3))})
    print(avg_scores)
    print("Mean score: ", np.mean(avg_scores))
    print(pred_table_mode)
    return(pred_table_mode)


#mode_pred_3 = ensemble_manual_predict(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"])
#mode_pred_3.to_csv("/path/to/project/data/clf_mode_pred_cl_n3_ni.csv")

#mode_pred_2 = ensemble_manual_predict(train, other, ["LogisticRegression", "RidgeClassifier", "GaussianNB", "MLPClassifier", "RandomForestClassifier"], n_cl=2)
#mode_pred_2.to_csv("/path/to/project/data/clf_mode_pred_cl_n2_ni.csv")
