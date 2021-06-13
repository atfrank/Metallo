from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import numpy as np
import random
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, roc_curve, auc, accuracy_score
import joblib
from collections import Counter
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import KMeans
import sys
import timeit
import matplotlib.pyplot as plt

def get_trainset(simhome, rna):
    simfile = "%s/%s_result.txt" %(simhome, rna.upper())
    with open(simfile, "r") as ifile:
        for line in ifile:
            if line.startswith("train"):
                trainset = eval(next(ifile, '').strip())
    trainset = [t.lower() for t in trainset]
    return trainset

def get_full_trainset(home):
    pdblist = "%s/data/pdblist.txt" %home
    return pd.read_csv(pdblist, delim_whitespace=True).columns.values

def resample(X_train, y_train, ratio=1):
    # downsample training set to balance
    print("resampling...")
    # sampler = ClusterCentroids(estimator=KMeans(n_clusters=5, random_state=0, verbose=10, n_jobs=20))
    sampler = RandomOverSampler(sampling_strategy=ratio, random_state=0)
    X_resampled, y_resampled = sampler.fit_sample(X_train, y_train)
    X_train, y_train = X_resampled, y_resampled

    print("After sampling data size in each category: ", sorted(Counter(y_train).items()))
    return X_train, y_train

def load_feature(feathome, rna, drop_dup=True, drop_O=False, drop_M=False):
    try:
        rna_feat = pd.read_csv("%s/%s.txt" %(feathome, rna), delim_whitespace=True)
    except:
        return
    # drop repeated rows
    if drop_dup:
        rna_feat = rna_feat.drop_duplicates(subset=rna_feat.columns[3:], keep="last")
    # drop Oxygen columns - there is no O in RNAs
    if drop_O:
        O_cols = ["0.5,O"] + [str(2**i)+",O" for i in range(7)]
        rna_feat = rna_feat.drop(O_cols, axis=1)
    if drop_M:
        rna_feat = rna_feat[rna_feat["resname"] != "M"]
    return rna_feat.iloc[:,3:], rna_feat['resname'].isin(["L", "M"]), rna_feat.iloc[:,:3] # X, y, labels


def get_feature_names(feathome, rna):
    try:
        rna_feat = pd.read_csv("%s/%s.txt" %(feathome, rna), delim_whitespace=True)
    except:
        return
    # drop Oxygen columns
    O_cols = ["0.5,O"] + [str(2**i)+",O" for i in range(7)]
    rna_feat = rna_feat.drop(O_cols, axis=1)
    return rna_feat.columns[3:]


def load_train(feathome, trainset, sample=0, ratio=1, verbose=1, add_dummy=True):
    X_train = []
    y_train = []
    random.seed(0)
    if verbose: print("loading data...")
    # load trainset data
    for rna in trainset:
        if verbose: print('loading training RNA:', rna)
        try:
            X, y, _ = load_feature(feathome, rna, drop_dup=True)
        except:
            continue
        X_train.extend(X.values)
        y_train.extend(y.values)
    print("data size in each category: ", sorted(Counter(y_train).items()))
    if add_dummy:
        n_dummy = int(0.01 * len(X_train))
        X_train.extend([X_train[0] * 0] * n_dummy)
        y_train.extend([False] * n_dummy)
        c = list(zip(X_train, y_train))
        random.shuffle(c)
        X_train, y_train = zip(*c)
    if sample:
        X_train, y_train = resample(X_train, y_train, ratio=ratio)
    return X_train, y_train


def initialize_model(n=100, depth=5, n_jobs=8, class_weight="balanced",
                     max_features="sqrt", min_samples_leaf=10, model="rf",
                     normalize=1, gridsearch=0):
    # bulid model and train

    if (gridsearch):
        # gridsearch parameters
        param_grid = {"n_estimators": [100, 500, 1000], "max_depth": [50, 100, 500], \
                "min_samples_leaf": [10, 50, 100], "max_features": [None, "sqrt", int(len(X_train[0])/3)]}

        # build a RF classifier
        classifier = GridSearchCV(
                RandomForestClassifier(class_weight='balanced', random_state=0, oob_score=True),
                cv = 5, param_grid = param_grid, return_train_score=True, verbose=10, n_jobs=n_jobs)

    else:
        if model == "rf":
            # classifier without gridsearch
            classifier = RandomForestClassifier(n_estimators=n,
                                                class_weight=class_weight,
                                                max_depth=depth,
                                                max_features=max_features,
                                                min_samples_leaf=10,
                                                random_state=0,
                                                oob_score=True,
                                                verbose=0,
                                                n_jobs=n_jobs)
        else:
            classifier = LogisticRegression(max_iter=n)

    if (normalize):
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
        scaler = StandardScaler()
        clf = Pipeline([("scaler", scaler), ("classifier", classifier)])
    else:
        clf = classifier

    return clf


def save_model(modelhome, clf, gridsearch=0, label=""):
    # save model
    model_fn = "%s/%s_model.pkl" %(modelhome, label)
    if (gridsearch):
        print(clf.best_params_)
        joblib.dump(clf.best_estimator_, model_fn)
    else:
        joblib.dump(clf, model_fn, compress=3)


def inference(clf, testset, feathome, outhome, drop_M=True):
    y_test = []
    y_pred = []
    test_df = pd.DataFrame()
    for rna in testset:
        try:
            X, y0, names = load_feature(feathome, rna, drop_dup=False, drop_M=drop_M)
        except:
            continue
        # label rna_feat with PDBID
        y1 = clf.predict_proba(X)
        np.savetxt("%s/%s.txt" %(outhome, rna), y1[:,1], fmt="%.3f")
        y_test.extend(y0.values)
        y_pred.extend(y1[:,1])

    return y_test, y_pred


def make_train_report(y_train_pred, y_train, reshome, label="train", threshold=0.5):
    # predict
    y_label = [y > threshold for y in y_train_pred]

    # Compute ROC curve and AUC
    roc = pd.DataFrame()
    roc['fpr'], roc['tpr'], roc['thresholds'] = roc_curve(y_train, y_train_pred)
    roc_auc = auc(roc['fpr'], roc['tpr'])
    roc.to_csv('%s/%s_roc.csv' %(reshome, label), index=None)
    print('train auc: %.3f' %roc_auc)

    # report performance
    cfm = confusion_matrix(y_train, y_label)
    print(cfm)
    print("f1 score: %.3f" %f1_score(y_train, y_label))
    # save statistics
    # metrics.to_csv('%s/cfm.csv' %reshome, index=['0','1'])
    return roc_auc


def make_report(reshome, y_train_pred, y_train, y_test, y_pred, label="", threshold = 0.5):
    # predict
    y_label = [y > threshold for y in y_pred]

    fpr, tpr, _ = roc_curve(y_train, y_train_pred)
    print("train auc: %.3f" %(auc( fpr, tpr ) ))

    # Compute ROC curve and AUC
    roc = pd.DataFrame()
    roc['fpr'], roc['tpr'], roc['thresholds'] = roc_curve(y_test, y_pred)
    roc_auc = auc(roc['fpr'], roc['tpr'])
    roc.to_csv('%s/%s_roc.csv' %(reshome, label), index=None)
    print('test auc: %.3f' %roc_auc)

    # report performance
    # print("Accuracy on training set: %.5f" %clf.score(X_train,y_train))
    # print("Accuracy on test set: %.5f" %accuracy_score(y_test, y_label))
    cfm = confusion_matrix(y_test, y_label)
    print(cfm)
    print("f1 score: %.3f" %f1_score(y_test, y_label))
    # save statistics
    metrics = pd.DataFrame()
    metrics['0'] = cfm[:,0]
    metrics['1'] = cfm[:,1]
    # metrics.to_csv('%s/cfm.csv' %reshome, index=['0','1'])
    return roc, roc_auc, cfm


def plot_roc(roc, roc_auc, home, penalty=50):
    fpr, tpr, thresholds = roc['fpr'], roc['tpr'], roc['thresholds']
    # calculate optimal threshold
    optimal_idx = np.argmax(tpr - penalty*fpr)
    optimal_threshold = thresholds[optimal_idx]
    print('ROC optimal threshold (optimal [fpr, tpr]): ', optimal_threshold, [fpr[optimal_idx], tpr[optimal_idx]])

    # plot ROC curve
    #plt.subplot(311)
    plt.plot(fpr, tpr, color='darkorange',
                     lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.scatter(fpr[optimal_idx], tpr[optimal_idx], color='navy',
                label='optimal threshold = %0.3f' %optimal_threshold)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.savefig('%s/plot/roc.png' %home)
    return optimal_threshold
