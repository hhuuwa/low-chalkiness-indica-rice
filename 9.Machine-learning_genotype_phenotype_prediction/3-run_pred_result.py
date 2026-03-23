import os, pickle, logging, pickle, joblib, sys, warnings, time
warnings.simplefilter('ignore')
from scipy import stats
import numpy as np
import pandas as pd
import lightgbm as lgb
import catboost as cb
from sklearn import ensemble, gaussian_process, linear_model, svm, metrics, model_selection
from sklearn.tree import DecisionTreeRegressor, ExtraTreeRegressor
from sklearn.neural_network import MLPRegressor
from scipy.stats import pearsonr, spearmanr

from mpire import WorkerPool

FILE_PATH = "/data2/zhoujb/project/hhf/250718/rawData/inputData4gene2K/"
OUT_PATH = "/data2/zhoujb/project/hhf/250718/Result/testRes4gene2K/"

def getPredRes(file_name, model_name):
    
    file_pre, target_col = file_name.split(".")[0].split("-")
    raw_data = pd.read_table(os.path.join(FILE_PATH, file_name), sep="\t", index_col=0)
    raw_data = raw_data.dropna(subset=target_col)
    feat_col = [x for x in raw_data.columns if x.startswith("fea_")]

    kf = model_selection.KFold(n_splits=5, shuffle=True, random_state=0)
    y_test_final, y_pred_final = [], []
    for i, (train_index, test_index) in enumerate(kf.split(raw_data)):
        print("Round {} start...".format(i))
        start_time = time.time()
    
        data_train = raw_data.iloc[train_index].copy()
        data_test = raw_data.iloc[test_index].copy()

        X_train = data_train[feat_col].copy()
        y_train = data_train[target_col].values.ravel()

        X_test = data_test[feat_col].copy()
        y_test = data_test[target_col].values.ravel()

        if model_name == "SVR":
            clf_model = svm.SVR()
            clf_model.fit(X_train, y_train)
        elif model_name == "GP":
            clf_model = gaussian_process.GaussianProcessRegressor(random_state=0)
            clf_model.fit(X_train, y_train)
        elif model_name == "RF":
            clf_model = ensemble.RandomForestRegressor(n_jobs=16, random_state=0)
            clf_model.fit(X_train, y_train)
        elif model_name == "LB":
            clf_model = lgb.LGBMRegressor(boosting_type="gbdt", random_state=0, n_jobs=16)
            clf_model.fit(X_train, y_train, callbacks=[lgb.log_evaluation(period=100)])
        elif model_name == "CB":
            clf_model = cb.CatBoostRegressor(random_state=0, thread_count=16, loss_function='RMSE')
            clf_model.fit(X_train, y_train, verbose=0, plot=False)
        elif model_name == "EN":
            clf_model = linear_model.ElasticNet(random_state=0)
            clf_model.fit(X_train, y_train)
        elif model_name == "Ridge":
            clf_model = linear_model.Ridge(random_state=0)
            clf_model.fit(X_train, y_train)
        elif model_name == "BR":
            clf_model = linear_model.BayesianRidge()
            clf_model.fit(X_train, y_train)
        elif model_name == "LR":
            clf_model = linear_model.LinearRegression()
            clf_model.fit(X_train, y_train)
        elif model_name == "ET":
            clf_model = ExtraTreeRegressor(random_state=0)
            clf_model.fit(X_train, y_train)
        elif model_name == "MLP":
            clf_model = MLPRegressor(random_state=0)
            clf_model.fit(X_train, y_train)
            
        y_pred = clf_model.predict(X_test)

        y_test_final.append(y_test)
        y_pred_final.append(y_pred)

        end_time = time.time()
        print("Round {}:{}".format(i, end_time-start_time))

    res_df = pd.DataFrame(columns=["Model", "xName", "yName", "Times", "Score", "Type"])
    for i in range(len(y_test_final)):
        score_pear = pearsonr(y_test_final[i], y_pred_final[i])[0]
        score_rmse = metrics.mean_squared_error(y_test_final[i], y_pred_final[i], squared=False)
        score_nrmse = score_rmse / np.std(y_test_final[i])

        res_df.loc[len(res_df)] = [model_name, file_pre, target_col, i+1, score_pear, "R"]
        res_df.loc[len(res_df)] = [model_name, file_pre, target_col, i+1, score_rmse, "RMSE"]
        res_df.loc[len(res_df)] = [model_name, file_pre, target_col, i+1, score_nrmse, "NRMSE"]

    with open(os.path.join(OUT_PATH, "{}-{}.pickle".format(model_name, file_name.split(".")[0])), "wb") as out_f:
        pickle.dump(res_df, out_f)

    print("{}-{} DONE".format(model_name, file_name.split(".")[0]))

    return

if __name__ == "__main__":
    file_list = os.listdir(FILE_PATH)
    model_list = ["SVR", "GP", "RF", "LB", "CB", "EN", "Ridge", "BR", "LR", "ET", "MLP"]
    param_list = []
    for i in model_list:
        for j in file_list:
            param_list.append((j, i))
    
with WorkerPool(n_jobs=4) as pool:
    pool.map(getPredRes, param_list, progress_bar=False)

    