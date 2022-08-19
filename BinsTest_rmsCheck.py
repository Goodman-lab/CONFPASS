#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:38:19 2022

@author: chingchinglam
"""

import pandas as pd
import itertools
from collections import defaultdict
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit import Chem

#import numpy as np
#import traceback
#import os



def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]



def get_raw_bins2(dft_df, bin_size=0.1):
    
    ## put conformers into bins according to the delG with the neighbours 
    ## cut off for the energy difference == 0.1 kcal mol-1
    
    delG_dft_df1=dft_df.sort_values('delG_kcal')
    
    sort_delG=delG_dft_df1['delG_kcal'].tolist()
    sort_idx=delG_dft_df1['idx'].tolist()
    
    
    sort_delG_v=sort_delG[0]
    sort_ls=[[sort_idx[0]]]
    sort_ls_cur=0
    for idx in range(1, len(sort_idx)):
        if sort_delG[idx]-sort_delG_v <= bin_size: 
            sort_ls[sort_ls_cur].append(sort_idx[idx])
        else:
            sort_ls.append([sort_idx[idx]])
                 
        sort_delG_v=sort_delG[idx]
        sort_ls_cur=len(sort_ls)-1   
        
    sort_ls2=[[idx, sort_ls[idx]] for idx in range(0,len(sort_ls))]
    
    return sort_ls2
        


def get_binary_combination(stuff): 
    
    binary_ls=[]
    for subset in itertools.combinations(stuff, 2):            
        binary_ls.append(subset)
                
    return binary_ls


def get_same_conformer_rdkit(test_ls,mols, mol_name):
    
   
    ## generate binary_combinations for conformers within the same bins by energy 
    bi_test_ls=get_binary_combination(test_ls)
    
    ## perform rms calculations for pairs of conformers 
    bi_test_ls_select=[]
    for i in bi_test_ls:
        align_rmsd=GetBestRMS(mols[i[0]],mols[i[1]])

        if align_rmsd <= 0.005:
            bi_test_ls_select.append([i]+[align_rmsd])
            
    return bi_test_ls_select


def put_item2list(cluster_ls,new_item):

    new_cluster_ls=cluster_ls
    new_TF=[]

    for i in new_cluster_ls:
        if new_item[0] not in i and new_item[1] not in i: 
            new_TF.append('F')
        
        elif new_item[0] in i and new_item[1] not in i: 
            new_TF.append('T1')
        
        elif new_item[0] not in i and new_item[1]  in i: 
            new_TF.append('T0')
            
        elif new_item[0] in i and new_item[1] in i: 
            new_TF.append('T')
            

    if 'T1' in new_TF: 
        append_idx=new_TF.index('T1')
        new_cluster_ls[append_idx].append(new_item[1]) 
    
    elif 'T0' in new_TF: 
        append_idx=new_TF.index('T0')
        new_cluster_ls[append_idx].append(new_item[0]) 
    
    elif 'T' in new_TF:
        pass

    else: 
        new_cluster_ls.append(new_item)

    return new_cluster_ls


def put_list2cluster(repeat_ls):

    ## resort the conformers into bins according to the rms calculation result 
    
    select_idx=[list(i[0]) for i in repeat_ls]
    
    new_cluster_ls=[]
    try:
        new_cluster_ls=[select_idx[0]]
        
    except IndexError:
        pass

    for i in select_idx[1:]:
        put_item2list(new_cluster_ls,i)
           
    return  new_cluster_ls


################################

class rmsCheck:
    
    def __init__(self, sdf_file_path, delG_dft_df, mol_name):
        
        self.dft_df=delG_dft_df
        self.mol_name=mol_name
        suppl = Chem.SDMolSupplier(sdf_file_path)
        self.mols = [mol for mol in suppl]

        
        self.delG_dft_ls=self.dft_df['delG_kcal'].tolist()
         
        
    def get_rmsCheck_cluster(self):
        
        ## use the following line for complete search
        #self.raw_cluster=[[0,list(range(0,len(self.dft_df)))]]
        
        self.raw_cluster=get_raw_bins2(self.dft_df, bin_size=0.1)
        
        repeat_conformers=[]
        for idx in range(0,len(self.raw_cluster)):
            test_ls=self.raw_cluster[idx][1]
            repeat_conformers.append(get_same_conformer_rdkit(test_ls,self.mols,self.mol_name))
            

        sort_repeat_conformers=[]
        for ls in repeat_conformers: 
            sort_repeat_conformers.append(put_list2cluster(ls))
    
        
        ## assemble the raw rmsCheck_cluster list
            
        conformers_with_repeats =  [x for y in sort_repeat_conformers for x in y]
        conformers_with_repeats2= [x for y in conformers_with_repeats for x in y]

        raw_cluster_idx=[i[1] for i in self.raw_cluster]
        complete_idx_ls=[x for y in raw_cluster_idx for x in y]
        conformers_with_NOrepeats = set(complete_idx_ls)-set(conformers_with_repeats2)
        
        ## raw rmsCheck_cluster list
        rmsCheck_cluster=conformers_with_repeats
        for i in conformers_with_NOrepeats:
            rmsCheck_cluster.append([i])
            
        ## sort the rmsCheck_cluster according to delG(DFT) value from lowest to highest
            
        rmsCheck_delG_ls=[]
        for i in rmsCheck_cluster:
            rmsCheck_delG_ls.append(sorted([[self.delG_dft_ls[j],j] for j in i]))
            
        ## rmsCheck_delG_ls=[[[0.0, 1], [0.975148986, 0], [0.9764040040000012, 4]],...
 
        rmsCheck_delG_idx_ls=sorted(rmsCheck_delG_ls)
        
        ## get the final sorted rmsCheck_cluster_ls 
        rmsCheck_idx_cluster_ls0 = []
        for i in rmsCheck_delG_idx_ls:
            rmsCheck_idx_cluster_ls0.append([j[1] for j in i])
        
        
        self.rmsCheck_idx_cluster_ls=[]
        for sub in rmsCheck_idx_cluster_ls0:
            sub_ls=[]
            for i in sub:
                sub_ls.append(i)
            self.rmsCheck_idx_cluster_ls.append(sub_ls)
            
            
        cluster_idx=[]
        for idx in range(0,len(self.rmsCheck_idx_cluster_ls)):
            cluster_idx_sub=[]
            for i in self.rmsCheck_idx_cluster_ls[idx]:
                cluster_idx_sub.append([i, idx])
            
            cluster_idx.append(cluster_idx_sub)
    
    
        unpack_cluster_idx_ls = sorted([x for y in cluster_idx for x in y])
        self.cluster_idx_dict= { unpack_cluster_idx_ls[idx][0]: unpack_cluster_idx_ls[idx][1] for idx in range(0, len(unpack_cluster_idx_ls))}

        self.total_conformers=len(self.cluster_idx_dict)
        self.total_clusters=len(self.rmsCheck_idx_cluster_ls)
        self.repeat_conformer_ratio= (self.total_conformers-self.total_clusters)/self.total_conformers
        
        
        ## return the rms calculation result 
        
        self.rmsCheck_cluster_df=pd.DataFrame({'name':[self.mol_name], 'conf_no': [self.total_conformers], 
                                               'clus_no': [self.total_clusters], 'ratio': [self.repeat_conformer_ratio],
                                               'clusters':[ self.rmsCheck_idx_cluster_ls], 'dict': [self.cluster_idx_dict]})
            
            
            
        