import os, sys
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
import joblib
from sklearn.metrics import roc_curve, auc
from sklearn.tree import export_graphviz

import matplotlib as mp
import matplotlib.pyplot as plt
import pylab as pl

FONTSIZE = 20
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : FONTSIZE}

mp.rc('font', **font)
train_data = pd.read_csv('../../ENHEP/datasets/hzz_vbf_ggf_train.csv')
test_data  = pd.read_csv('../../ENHEP/datasets/hzz_vbf_ggf_test.csv')
print(train_data[:10])
print("len(train_data): %10d" % len(train_data))


def plotData(data, xmin=0, xmax=8, ymin=0, ymax=2000, N=1000, ftsize=FONTSIZE):
    
    # set size of figure
    plt.figure(figsize=(8, 8));
    
    # get axis info
    axes = plt.gca()
    # set axes' limits
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax)
    
    # annotate axes
    plt.xlabel(r'$|\Delta\eta|_{jj}$', fontsize=24)
    plt.ylabel(r'$m_{jj}$ (GeV)', fontsize=24)
    
    # split into sig and bkg for the purposes of plotting
    sig = data[data.target > 0.5][:N]
    bkg = data[data.target < 0.5][:N]

    plt.scatter(sig.detajj, sig.massjj, marker='o',
                s=50, c='blue', alpha=0.3, label='VBF')
    pl.legend(loc='upper left') # activate legend
    
    plt.scatter(bkg.detajj, bkg.massjj, marker='*',
                s=100, c='red',  alpha=0.3, label='ggF')
    pl.legend(loc='upper left') # activate legend
    plt.tight_layout()
    plt.savefig('higgs_vbf_ggf_variables.png')
    plt.show()
    
plotData(train_data)

def standardize_data(train_data, test_data, inputs):
    scaler  = StandardScaler()
    scaler.fit(train_data[inputs])
    
    X_train = scaler.transform(train_data[inputs])
    X_test  = scaler.transform(test_data[inputs])
    y_train = np.array(train_data['target'])
    y_test  = np.array(test_data['target'])

    return (X_train, X_test, y_train, y_test, scaler)
    
inputs = ['detajj', 'massjj']
X_train, X_test, y_train, y_test, scaler = standardize_data(train_data,
                                                            test_data,
                                                            inputs)
print(type(X_train), type(y_train))

dt  = DecisionTreeClassifier(max_depth=6,
                             min_samples_leaf=1)

ada = AdaBoostClassifier(dt,
                         n_estimators=100,
                         random_state=1,
                         algorithm="SAMME")

# possible combinations of params
params = {'base_estimator__min_samples_leaf': [100, 200, 300, 400, 500],
          'n_estimators' : [100, 200, 300, 400, 500]}
print ("============TRAINING===========")
ntrain = 5000
rcv = RandomizedSearchCV(ada, params, n_iter=20, verbose=1)
rcv.fit(X_train[:ntrain], y_train[:ntrain])

print ("Best set of parameters: %s" % rcv.best_params_)

adabest = rcv.best_estimator_
adabest.fit(X_train, y_train)

print("Training set score: %f" % adabest.score(X_train, y_train))
print("Test set score:     %f" % adabest.score(X_test,  y_test))

filename = 'higgs_vbf_ggf_bdt.pkl'
print( "save to %s" % filename)
joblib.dump([adabest, scaler], filename)
