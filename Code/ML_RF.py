import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
import time
import itertools
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score


############################################################################## Data ##############################################################################

data=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTrain0710.csv",header=0,sep=",")


train_x=data.drop(["ID","Target"],axis=1)
train_y=data["Target"]

IRFTest=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeTest0710.csv",header=0,sep=",")

IRFtest_x=IRFTest.drop(["ID","Target"],axis=1)
IRFtest_y=IRFTest["Target"]


############################################################################## Sklearn ##############################################################################
#
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
        accuracyList.append(RandomForest(i))
        print("acc: "+str(RandomForest(i)))
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


def RandomForest(parameter):
    randomforest = RandomForestClassifier(n_estimators=parameter[0], criterion=parameter[1], oob_score="True",
                                         max_depth=parameter[2], bootstrap="True", random_state=5)
    if (parameter[3] == True):

        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = pd.DataFrame(min_max_scaler.fit_transform(train_x))
        kfold=StratifiedKFold(n_splits=10, random_state=5, shuffle=True)
        scores = cross_val_score(randomforest, X_train_minmax, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()

    else:
        kfold = StratifiedKFold(n_splits=10, random_state=5, shuffle=True)
        scores = cross_val_score(randomforest, train_x, train_y, cv=kfold, scoring="accuracy")
        return scores.mean()

def RFParameter(n_estimator,criterion,max_depth,scalling):
    return [[a,b,c,d] for a,b,c,d in itertools.product(n_estimator,criterion,max_depth,scalling)]



def BestModel(Parameter):
    randomforest = RandomForestClassifier(n_estimators=Parameter[0], criterion=Parameter[1], oob_score="True",
                                         max_depth=Parameter[2], bootstrap="True", random_state=5)
    if(Parameter[3]==True):
        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = min_max_scaler.fit_transform(train_x)
        randomforest.fit(X_train_minmax,train_y)
    else:
        randomforest.fit(train_x,train_y)
    p = randomforest.predict(IRFtest_x)
    accu=accuracy_score(IRFtest_y, p)
    matrix = pd.crosstab(IRFtest_y, p)
    return accu,matrix


os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData")

allData=pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/MaizeData/MaizeAll0710.csv",header=0,sep=",")
allData_x=data.drop(["ID","Target"],axis=1)
allData_y=data["Target"]






def Prediction(Parameter,chr):
    pre = pd.read_table("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/SorghumData/predition/chr%s.csv"%(chr),
                        header=0, sep=",")
    pre_x = pre.drop(["ID"], axis=1)
    randomforest = RandomForestClassifier(n_estimators=Parameter[0], criterion=Parameter[1], oob_score="True",
                                         max_depth=Parameter[2], bootstrap="True", random_state=5)
    if(Parameter[3]==True):
        min_max_scaler = preprocessing.MinMaxScaler()
        X_train_minmax = min_max_scaler.fit_transform(allData_x)
        randomforest.fit(X_train_minmax,allData_y)
    else:
        randomforest.fit(allData_x,allData_y)
    predi = randomforest.predict(pre_x)
    p = pd.DataFrame(predi, columns=['predection'])
    r = pd.DataFrame(pre['ID']).join(p)
    r.to_csv('chr%s_pre.csv'%(chr),index=None)



n_estimator=[10,50,200]
criterion=["entropy","gini"]
max_depth=[None,2,10,100]
data_scale=[True,False]

print("                                                                                                     ")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("================================== Model : Random Forest ===========================================")
RFP=RFParameter(n_estimator,criterion,max_depth,data_scale)
RFacc=RunModel(RFP)
print("totalTime:")
print(time.clock()-start)


print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("================================== Best Model ===========================================")
bestParameter=RFP[RFacc.index(max(RFacc))]
print("Use Best Parameter: "+ str(bestParameter))

accu,matrix=BestModel(bestParameter)
print("accuracy: "+ str(accu))
print("confusion matrix:")
print(matrix)
