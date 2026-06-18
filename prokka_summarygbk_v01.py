#!/usr/bin/env python3
# prokka_summarygbk_v01.py
# Author : Tassy Joseph-S. Bazile, Bioinformatian ApHL Bioinformatics Fellow 2021-2023  
#sys.path.append('/usr/sammy/') to get append the path of a module
import os
import pandas as pd
import numpy as np
import re
import xlsxwriter
from Bio import SeqIO
from collections import Counter
import argparse

parser = argparse.ArgumentParser(prog = 'prokka_summarygbk_v01.py [Options]')
parser.add_argument('--pkanres', type=str, help =' provide the path to prokka output directory ', required=True)
args = parser.parse_args()
pk_dir = args.pkanres

"""
Business Logic:
The module creates the summary report of the prokka annotation.

Usage:
python prokka_summarygbk_v01.py --pgapres

"""
features_set = set()
# Getting report as a dataframe
def pgap_report(s,path_to_gbk):

    # SeqIO.parse(*.gbk) used for a mutiple-record file; SeqIO.read(*.gbk) for a single record file
    for seq_record in SeqIO.parse(path_to_gbk, "genbank"):
        
        """
        Everything below is for each single record
        From record, getting features
        feature_Rec = seq_record.features
        From features, getting types
        """
        feature_Rec = seq_record.features
        all_feature_types = [feature.type for feature in feature_Rec]
        
        feature_type = set(all_feature_types)
        
        feature_counts = Counter(all_feature_types) # a dictionary mapping feature and number "for each record"
        
        # Looking at the features
        feature_counts.keys()
        # feature_counts.keys(): dict_keys(['source', 'gene', 'CDS', 'tRNA']) ... for each record 
        #Removing source       
        del feature_counts['source']
        # List of Features
        features_keylist = list(feature_counts.keys()) # list of features per record
        [features_set.add(it) for it in features_keylist] # features for all records

        # Adding the dictionary feature_count to dataframe
        df_features = pd.DataFrame(list(feature_counts.items())) # creating df from dictionary
        df_featuresT =  df_features.transpose()
        header_row = df_featuresT.iloc[0] # converting row into column header
        df_features2 = pd.DataFrame(df_featuresT.values[1:], columns =  header_row)  
        df_features2.insert(0,"SampleID", s, True)  # adding at column the position zero
        return df_features2

# parsing df to csv and txt files   
def report_csv(df_features2):
    df_features2.to_csv('prokka_annotation_report.txt', sep=' ', index=True, header=True)
    df_features2.to_csv('prokka_annotation_report.csv', sep=' ', index=True, header=True)
    with pd.ExcelWriter('prokka_annotation_report.xlsx', engine='xlsxwriter') as writer: # writer is to add more than one sheet in the workbook
        df_features2.to_excel(writer, sheet_name='Annotation_Report', startrow=1)
        #writer.save()
    return
                       
def main():
    pk_resPath = pk_dir
    sampledir_list = os.listdir(pk_dir)
    sampledir_list.sort()
    df_list = []
    for s in sampledir_list:  
        saved_path = pk_resPath +'/'+ s 
        rep_file = s.split("_")[0]+'.gbk'
        path_to_gbk = os.path.join(saved_path, rep_file)
        feature_df = pgap_report(s, path_to_gbk)
        df_list.append(feature_df)
        feature_dfn = pd.concat(df_list)
        feature_dfn.index = np.arange(1, len(feature_dfn) + 1)
    #print(feature_dfn)
        #print(s)
    report_csv(feature_dfn)
if __name__ == "__main__":
    main()

