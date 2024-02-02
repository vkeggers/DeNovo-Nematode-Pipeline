"""
Created on Thu Feb 01 14:19:06 2024
@author: vkeggers and virallyDanny
@description: Compares the EDTA output GFF3 file against the output GFF file of a gene
              annotation software (such as Braker3) to detect EDTA-predicted TEs found
              inside predicted gene sequences. These sequences are filtered from the
              EDTA annotation GFF file and outputted to a rejected_TEs.gff file. The
              non-conflicting TE prediction are outputted to a filtered_TEs.gff file.

@Usage: EDTA_Gene_Annotation_Comparison.py --EDTA [EDTAoutput.gff] --gene [gene_annotation_output.gff]
"""

import pandas as pd
import numpy as np
import xgboost
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split, KFold, cross_val_score
import matplotlib.pyplot as plt
import argparse
import warnings
warnings.filterwarnings('ignore')

# Create an argument parser object to parse the value of the --Genus argument to the genus variable
# parser = argparse.ArgumentParser(description='XGBoost model to predict the origin of a TE')
# parser.add_argument('--Genus', type=str, help='Genus of the TE')
# args = parser.parse_args()
# genus = args.Genus

genus = 'Oscheius'

# Load SIDR stats tsv file to a pandas dataframe titled data
stats = pd.read_csv('./data/SIDRstats.tsv', sep='\t', header=0)

# Make a new dataframe titled trainingDF that only contains rows which do not equal 'No hits found' in the 'Origin' column
trainingDF = stats[stats['Origin'] != 'No hits found']

# Make a new dataframe titled orig where if the 'Origin' column equals the 'Genus' string variable. If it is True then change the Origin value to True, else False
orig = trainingDF[['Origin']]
orig['Origin'] = trainingDF['Origin'].str.contains(genus).astype(int)

# Make a variable titled "Train" that samples trainingDF with a random sampling 1/3 of the data in the dataframe
Train = trainingDF.sample(frac=1/3)

# Make a variable titled stat1_test which takes all the rows not in Train and keeps all columns except the 'contig' and 'Origin' columns
stat1_test = trainingDF[~trainingDF.index.isin(Train.index)].drop(['contig', 'Origin'], axis=1)

stat1_train = trainingDF[trainingDF.index.isin(Train.index)].drop(['contig', 'Origin'], axis=1)

origin_test = orig[~orig.index.isin(Train.index)]

origin_train = orig[orig.index.isin(Train.index)]

####################
### WORKS SO FAR ###
####################

model = XGBClassifier(n_estimator=5000, max_depth=2, random_state=3, objective='binary:logistic', distribution='bernoulli')

model.fit(stat1_train, origin_train)

kfold = KFold(n_splits=5, random_state=7, shuffle=True)

results = cross_val_score(model, stat1_train, origin_train, cv=kfold)

print("Accuracy: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))

# Print SIDR predictions to stats dataframe
stats['SIDR_predictions'] = model.predict(stats.drop(['contig', 'Origin'], axis=1))

# Save the stats dataframe to a new tsv file
stats.to_csv('Full_SIDR_predictions.tsv', sep='\t', index=False)

# Save the stats dataframe to a new tsv file if the SIDR_predictions column equals 1
stats[stats['SIDR_predictions'] == 1].to_csv('Kept_SIDR_predictions.tsv', sep='\t', index=False)
stats[stats['SIDR_predictions'] == 0].to_csv('Removed_SIDR_predictions.tsv', sep='\t', index=False)


# Feature Importance
feature_importance = model.feature_importances_features = stats.columns

# Plot the top 7 features
xgboost.plot_importance(model, max_num_features=12)

# Prediction Report
y_pred = model.predict(stat1_test)
report = classification_report(origin_test, y_pred)

# Show the plot
plt.show()