import pandas as pd
import os
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
import time
import itertools
from sklearn.model_selection import KFold
import multiprocessing
from multiprocessing import Pool
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold

############################################################################## Data ##############################################################################

data=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTrain0710.csv",header=0,sep=",")


train_x=data.drop(["ID","Target"],axis=1)
train_y=data["Target"]

IRFTest=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTest0710.csv",header=0,sep=",")

IRFtest_x=IRFTest.drop(["ID","Target"],axis=1)
IRFtest_y=IRFTest["Target"]

############################################################################## Data ##############################################################################


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
        accuracyList.append(NN(i))
        print("acc: "+str(NN(i)))
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



def NNParameter(hiddenlayer,activation,solver,learning_rate,learning_rate_init,scaling):
    return [[a,b,c,d,e,f] for a,b,c,d,e,f in itertools.product(hiddenlayer,activation,solver,learning_rate,learning_rate_init,scaling)]

def NN(parameter):
    dnn = MLPClassifier(hidden_layer_sizes=parameter[0],activation=parameter[1],solver=parameter[2], alpha=1e-5, learning_rate=parameter[3],learning_rate_init=parameter[4],validation_fraction=0.1,random_state = 5,batch_size=50,early_stopping=True)
    if (parameter[5] == True):

        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = pd.DataFrame(min_max_scaler.fit_transform(train_x))
        kfold=StratifiedKFold(n_splits=10, random_state=3, shuffle=True)
        scores = cross_val_score(dnn, X_train_minmax, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()

    else:
        kfold = StratifiedKFold(n_splits=10, random_state=3, shuffle=True)
        scores = cross_val_score(dnn, train_x, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()

def BestModel(parameter):
    dnn = MLPClassifier(hidden_layer_sizes=parameter[0],activation=parameter[1],solver=parameter[2], alpha=1e-5, learning_rate=parameter[3],learning_rate_init=parameter[4],validation_fraction=0.1,random_state = 5,batch_size=50,early_stopping=True)
    if(parameter[5]==True):
        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = min_max_scaler.fit_transform(train_x)
        dnn.fit(X_train_minmax,train_y)
        p = dnn.predict(IRFtest_x)

    else:
        dnn.fit(train_x,train_y)
        p = dnn.predict(IRFtest_x)
    accu=accuracy_score(IRFtest_y, p)
    matrix = pd.crosstab(IRFtest_y, p)
    return accu,matrix


hiddenlayer=[(1,50),(2,200),(50,50),(100,100)]
activation=["tanh","relu"]
solver=["lbfgs", "adam"]
learning_rate=["constant","adaptive"]
learning_rate_init=[0.0001]
scaling=[True,False]

# hiddenlayer=[(2,2)]
# activation=["relu"]
# solver=["adam"]
# learning_rate=["adaptive"]
# learning_rate_init=[0.01]


print("                                                                                                     ")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("================================== Model : Neural Network ===========================================")
nnP=NNParameter(hiddenlayer,activation,solver,learning_rate,learning_rate_init,scaling)
NNacc=RunModel(nnP)
print("totalTime:")
print(time.clock()-start)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("================================== Best Model ===========================================")
# bestParameter=[(50, 50), 'tanh', 'adam', 'constant', 0.0001, True]
bestParameter=nnP[NNacc.index(max(NNacc))]
print("Use Best Parameter: "+ str(bestParameter))

accu,matrix=BestModel(bestParameter)
print("accuracy: "+ str(accu))
print("confusion matrix:")
print(matrix)
