#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:07:14 2022

@author: chingchinglam

linking the pre-clustering process scripts


"""




import pandas as pd 
import traceback

from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict


########### linking the pre-clustering process scripts ########### 

from isolate_key_dihedral_v5 import isolate_dihedral
from cal_dihedral_v2 import dihedral_df
from dihedral_parameter_v2 import  remove_fixed_dihedrals, remove_fixed_bond_df
from correcting_dihedral_v1 import correction_by_gap

########## 


def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]


def get_cluster_df(path):
    
    sdf_path=path
    
    ## pipeline to extract the descriptors 
    
    select_dihedral = isolate_dihedral(sdf_path)
    raw_df = dihedral_df(sdf_path, select_dihedral[0], select_dihedral[1])
    tobe_filtered_dihedral = remove_fixed_dihedrals(raw_df, select_dihedral[0], select_dihedral[1])[2]
    dihedral_df_final=remove_fixed_bond_df(raw_df,tobe_filtered_dihedral)[0]
    dihedral_df_final1=correction_by_gap(dihedral_df_final)[0]
    
    descriptor=dihedral_df_final1.to_numpy()
    
    ## create n_ls
    cluster_ls=list(range(1, len(descriptor)+1, 1))
    
    ## analyse and append result into a large df 
    
    nor_n_cluster=[i/cluster_ls[-1] for i in cluster_ls]
    cluster_result_sort_ls =[]

    for k in cluster_ls:
        model = AgglomerativeClustering(n_clusters=k).fit(descriptor)
        
        
        cluster_result_sort=list_duplicates(model.labels_)
        cluster_result_sort_ls.append(cluster_result_sort)
        
    
    result_dict ={'n':cluster_ls, 'nor_n': nor_n_cluster, 'clusters': cluster_result_sort_ls}
    result_df = pd.DataFrame(result_dict)
    
    name=path.split('/')[-1]
    result_df['name'] = name
    first_column = result_df.pop('name')
    result_df.insert(0, 'name', first_column)
    
    return result_df


def get_cluster_df_multi(mol_path_ls):

    ## perform get_cluster_df() for multiple files in the same directory 

    result_df_ls=[]
    for i in mol_path_ls:
        try:
            result_df_ls.append(get_cluster_df(i))
    
        except Exception:
            print('Error: '+i)
            traceback.print_exc()
    

    result_df=pd.concat(result_df_ls)
    
    
    
    return result_df





