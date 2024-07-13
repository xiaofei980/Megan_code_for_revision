# Import packages
import pandas as pd
import os 
from scipy.spatial.distance import cdist
from datetime import datetime 

# Get current directory
current_directory = os.getcwd()
date = datetime.now().strftime("%Y%m%d")
cutoff = 25
# Find datafile

path = f"{current_directory}/Results/Figures_output/20240626_allROIs_25px_Iris_neighbours_distanceCalculation_allCells.csv"
# path = f"{current_directory}/Data/distances/dists_cellpairs_25px_20240229.csv"
dists = pd.read_csv(path)

# Convert datafile of distances to relative frequency of source and target cell per source_id and per patient
count_df = dists.groupby(['ROI_ID','source_ID','source_cluster', 'target_cluster']).size().reset_index(name='count')
total_counts_scid = count_df.groupby(['ROI_ID','source_ID','source_cluster'])['count'].sum().reset_index(name='total_scid')


count_df = pd.merge(count_df, total_counts_scid, on=['ROI_ID','source_ID' ,'source_cluster'])
totalscid = count_df['total_scid']
# print(count_df.head())
count_df['relative_frequency'] = (count_df['count'] / count_df['total_scid'])

# Convert dataframe to correct format for clustering
sc = dists['source_cluster'].unique()
freq = count_df.pivot(columns=['target_cluster'], index=['ROI_ID','source_ID','source_cluster'], values=['relative_frequency'])
freq.reset_index(inplace=True)
freq.columns  = freq.columns.droplevel(0)
freq.columns = ['ROI_ID','source_ID','source_cluster',*freq.columns[3:]]
freq.fillna(0, inplace=True)

print('Dataset converted, below is a preview:\n',freq.head())

# Save file
output_path = f"{current_directory}/Results/Figures_output/"

# totalscid.to_csv(f"{output_path}frequencies_cells_{cutoff}px_{date}.csv", index = False)
freq.to_csv(f"{output_path}frequencies_{cutoff}px_{date}.csv", index = False)


