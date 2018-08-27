import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
import time
import itertools
from sklearn.model_selection import KFold
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold

import multiprocessing
from multiprocessing import Pool
############################################################################## Data ##############################################################################


data=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTrain0710.csv",header=0,sep=",")


train_x=data.drop(["ID","Target"],axis=1)
train_y=data["Target"]

IRFTest=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTest0710.csv",header=0,sep=",")

IRFtest_x=IRFTest.drop(["ID","Target"],axis=1)
IRFtest_y=IRFTest["Target"]

############################################################################## Sklearn ##############################################################################


def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
warnings.filterwarnings("ignore")
start=time.clock()


def RunModel(ParameterList):
    accuracyList=[]
    for i in ParameterList:
        trainingProcess=str(ParameterList.index(i) + 1) + "/" + str(len(ParameterList))
        print("++++++++++++++++++++++++++++++ %s ++++++++++++++++++++++++++++++++++"%(trainingProcess))
        print(i)
        single=time.clock()
        accuracyList.append(Adaboost(i))
        print("acc: "+str(Adaboost(i)))
        print("currentBest_Acc: "+str(max(accuracyList))+" #"+str(accuracyList.index(max(accuracyList))))
        print("Time"+str(time.clock()-single))
        print("++++++++++++++++++++++++++++++ Finished ++++++++++++++++++++++++++++++++++")
    print("Accuracy:"+" "+str(accuracyList))
    print("++++++++++++++++++++++++++++++Summary++++++++++++++++++++++++++++++++++")
    print("Best Accuracy:"+" "+str(max(accuracyList)))
    print("Number of Experiments done:"+" "+str(len(ParameterList)))
    print("ParameterList:"+" "+str(ParameterList))
    print("Best Parameter Combination:"+" "+str(ParameterList[accuracyList.index(max(accuracyList))]))
    print("++++++++++++++++++++++++++++++Done++++++++++++++++++++++++++++++++++")
    return accuracyList


def Adaboost(parameter):
    AdaBo = AdaBoostClassifier(DecisionTreeClassifier(criterion=parameter[0], max_depth=parameter[1],max_features=parameter[2]), n_estimators=parameter[3],
                               learning_rate=parameter[4],random_state=5)

    if (parameter[5] == True):

        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = pd.DataFrame(min_max_scaler.fit_transform(train_x))
        kfold=StratifiedKFold(n_splits=10, random_state=3, shuffle=True)
        scores = cross_val_score(AdaBo, X_train_minmax, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()

    else:
        kfold = StratifiedKFold(n_splits=10, random_state=3, shuffle=True)
        scores = cross_val_score(AdaBo, train_x, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()



def BestModel(parameter):
    AdaBo = AdaBoostClassifier(
        DecisionTreeClassifier(criterion=parameter[0], max_depth=parameter[1], max_features=parameter[2]),
        n_estimators=parameter[3],
        learning_rate=parameter[4], random_state=5)
    if(parameter[3]==True):
        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = min_max_scaler.fit_transform(train_x)
        AdaBo.fit(X_train_minmax,train_y)
        p = AdaBo.predict(IRFtest_x)
    else:
        AdaBo.fit(train_x,train_y)
        p = AdaBo.predict(IRFtest_x)
    accu=accuracy_score(IRFtest_y, p)
    matrix = pd.crosstab(IRFtest_y, p)
    return accu,matrix



criterion=["entropy","gini"]
max_depth=[None,50,100]
max_features=["auto"]
n_estimators=[20,200]
learning_rate=[0.0001,0.001]
scaling=[True,False]

def AdaboostParameters(criterion,max_depth,max_features,n_estimators,learning_rate,scaling):
    return [[a,b,c,d,e,f] for a,b,c,d,e,f in itertools.product(criterion,max_depth,max_features,n_estimators,learning_rate,scaling)]

print("                                                                                                     ")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# print("================================== Model : Neural Network ===========================================")
AdaP=AdaboostParameters(criterion,max_depth,max_features,n_estimators,learning_rate,scaling)
print(AdaP)
print(AdaP[32])
# Adaboostacc=RunModel(AdaP)
# print("totalTime:")
# print(time.clock()-start)
#
# print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("================================== Best Model ===========================================")
bestParameter=['gini', 50, 'auto', 20, 0.0001, True]
# bestParameter=AdaP[Adaboostacc.index(max(Adaboostacc))]
# bestParameter=[10, 'entropy', None, True]
print("Use Best Parameter: "+ str(bestParameter))
test_accu,matrix=BestModel(bestParameter)
print("accuracy: "+ str(test_accu))
print("confusion matrix:")
print(matrix)

print("================================== Best Model ===========================================")

# ================================== Best Model ===========================================
# Use Best Parameter: ['gini', 50, 'auto', 20, 0.0001, True]
# accuracy: 0.9181890389197777
# confusion matrix:
# col_0   DTA  DTC  DTH  DTM  DTT  NonTIR
# Target
# DTA     186    0    6    1    0       7
# DTC       7  177    2    0    2       9
# DTH       3    0  185    4    7       1
# DTM       3    1    7  188    0       1
# DTT       5    1   12    2  180       0
# NonTIR    6    1   10    3    2     240
# ================================== Best Model ===========================================



