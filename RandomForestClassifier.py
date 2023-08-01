# Random Forest Classifier is a machine learning approach that is used to predict the phenotype (diease or normal) based on the count values/expression levels.
# One dataset is required to test and train the model and an external independent dataset is required to test the model.

from sklearn.metrics import recall_score, accuracy_score, f1_score, precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import pandas as pd

clf = RandomForestClassifier(n_estimators=10) #initialise model

def read_dataframe(file, value):
    
    df = pd.read_csv(file, delimiter= ';', index_col= 0)
    df = df.T
    df['y'] = value
    
    return(df) 

PDAC = read_dataframe('/Users/khadija/Desktop/ML/biomarker_data/tumour_biomarker.csv', 1)#assign pdac to 1, this is the phenotype

Normal = read_dataframe('/Users/khadija/Desktop/ML/biomarker_data/normal_biomarker.csv', 0) #assign normal to 0

def merge(Normal,PDAC):
    com = [Normal,PDAC]
    combine = pd.concat(com)
    return (combine)

df = merge(Normal,PDAC) # combine

X = df.iloc[:, :-1] #Deletes last column 

y = df['y'] #Select last column

y.head(2) # y values (0 and 1)

X.head(2) # counts and genes

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, shuffle=True) #Split train and test (shuffle values)  

clf_trained_model = clf.fit(X_train, y_train) #Model Training on training data (70% of data)

y_predicted = clf_trained_model.predict(X_test) #Test model on test data (30% data)

recall_score(y_test, y_predicted) #Recall Specificity Score

accuracy_score(y_test, y_predicted) #Accuracy Score 

f1_score(y_test, y_predicted) #F1-Score

precision_score(y_test, y_predicted) #Precision Score

#test the model on an independent count matrix
#Prepare the dataframe as the intail steps, however, do not split the frame.

PDAC_D1 = read_dataframe('/Users/khadija/Desktop/dataset_1/ML/D1_updwn/d1_tumour.csv', 1)#assign pdac to 1, this is the phenotype
Normal_D1 = read_dataframe('/Users/khadija/Desktop/dataset_1/ML/D1_updwn/d1_normal.csv', 0) #assign normal to 0

df_1 = merge(Normal_D1,PDAC_D1) # combine

X_1 = df_1.iloc[:, :-1] #Deletes last column 

y_1 = df_1['y'] #Select last column

y_predicted = clf_trained_model.predict(X_1) #Test model on test data (100% data of the new dataset)

recall_score(y_1, y_predicted) #Recall Specificity Score #take note of the y_1, it is the phenetype for the external dataset

accuracy_score(y_1, y_predicted) #Accuracy Score

f1_score(y_1, y_predicted) #F1-Score

precision_score(y_1, y_predicted) #Precision Score

