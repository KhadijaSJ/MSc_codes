# Random Forest Classifier is a machine learning approach that is used to predict the phenotype (diease or normal) based on the count values/expression levels.
# One dataset is required to test and train the model and an external independent dataset is required to test the model.

from sklearn.metrics import recall_score, accuracy_score, f1_score, precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import pandas as pd

