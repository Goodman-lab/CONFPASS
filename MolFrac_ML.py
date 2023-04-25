#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 18:03:50 2022

@author: chingchinglam
"""


import pandas as pd
import numpy as np
import math
#import ast

#import sys
#import json

from collections import defaultdict



def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]




def get_delG_df(single_p_ls,csv_df,ideal_dict,temperature=298.15):
    
    
    ## get delG df for a selection of reopt conformers 
    
    #print(single_p_ls)
    
    p_ls_clno=[ideal_dict.get(i) for i in single_p_ls]
    p_ls_clno_dup=[i[1] for i in list_duplicates(p_ls_clno)]
    p_ls_single_idx=[i[0] for i in p_ls_clno_dup]
   

    filter_single_p_ls=[single_p_ls[idx] for idx in p_ls_single_idx]
    #print(filter_single_p_ls)

    select_conf=csv_df.loc[csv_df['idx'].isin(filter_single_p_ls)] 
       
    
    
    
    min_G=min(select_conf['G_hatree'])

    select_conf_df=select_conf.copy()

    select_conf_df['delG_kJ']=(select_conf_df['G_hatree'] - min_G)*2625.5
    mol_fraction = [math.exp(-((g/(temperature*8.314)))*1000) for g in select_conf_df['delG_kJ'].tolist()]
    select_conf_df['mol_fraction']=mol_fraction
    
    return select_conf_df
        


def get_n0(mol_ls, cutoff=1):
    
    ## given the result, newly introduced mol fraction list
    ## find n0 - the n value at which molfrac no longer deviates from x
    
    n_ls=list(range(1,len(mol_ls)+1))
    
    min_ls_re=mol_ls[::-1]
    #n_ls_re=n_ls[::-1]
    
    n0_re=len(n_ls)
    
    for idx in range(0,len(min_ls_re)):
        if min_ls_re[idx] > cutoff:
            n0_re=idx
            break
    
    n0=len(n_ls)-n0_re
    
    nor_n0=n0/len(n_ls)
    
    return n0, nor_n0



def get_new_mol_data(p_ls,csv_df,ideal_dict,T=298.15):
    
    ## go through from reopt conf = 1 to n 
    
    #print(p_ls)
    
    sel_df_ls=[]
    for single_p_ls in p_ls:
        sel_df=get_delG_df(single_p_ls,csv_df,ideal_dict,temperature=T)
        sel_df_ls.append(sel_df)
        
    del_mol_ls=[1]
    
    for idx in range(1,len(sel_df_ls)):
        
        ## take into the consideration of conformers with repeat structure 
        
        try:
            del_idx=list(set(sel_df_ls[idx]['idx'].tolist())-set(sel_df_ls[idx-1]['idx'].tolist()))[0]
            
            curr_df=sel_df_ls[idx]
            mol_frac=curr_df[curr_df['idx']==del_idx]['mol_fraction'].tolist()[0]
            #print(curr_df[curr_df['idx']==del_idx])
            
        except IndexError:
            mol_frac=0
            
        del_mol_ls.append(mol_frac)
        
    
    ## calculate the descriptors 
    
    
    _,n0_molfac = get_n0(del_mol_ls, cutoff=1)
    
    len_del_mol_ls = len(del_mol_ls)
    
    
    
    del_mol_ls_0 = len([i for i in del_mol_ls if i == 0])/len_del_mol_ls
    del_mol_ls_1 = len([i for i in del_mol_ls if i == 1])/len_del_mol_ls
    del_mol_ls_frac1=len([i for i in del_mol_ls if i <= 0.1])/len_del_mol_ls
    del_mol_ls_frac2=len([i for i in del_mol_ls if i <= 0.2])/len_del_mol_ls
    del_mol_ls_frac5=len([i for i in del_mol_ls if i <= 0.5])/len_del_mol_ls
    
    del_mol_ls_per = del_mol_ls[round(len_del_mol_ls*0.6):]
    len_del_mol_ls_per=len(del_mol_ls_per)
    
    if len_del_mol_ls_per == 0:
        del_mol_ls_per_0 = 0      
        del_mol_ls_per_1 = 0
        del_mol_ls_per_frac1=0
        del_mol_ls_per_frac2=0
        del_mol_ls_per_frac5=0 
        
    else:
        del_mol_ls_per_0 = len([i for i in del_mol_ls_per if i == 0])/len_del_mol_ls_per         
        del_mol_ls_per_1 = len([i for i in del_mol_ls_per if i == 1])/len_del_mol_ls_per 
        del_mol_ls_per_frac1=len([i for i in del_mol_ls_per if i <= 0.1])/len_del_mol_ls_per 
        del_mol_ls_per_frac2=len([i for i in del_mol_ls_per if i <= 0.2])/len_del_mol_ls_per 
        del_mol_ls_per_frac5=len([i for i in del_mol_ls_per if i <= 0.5])/len_del_mol_ls_per 
        
        
    
    
    
    result_df=pd.DataFrame({'n0_molfac':[n0_molfac], 'f0':[del_mol_ls_0], 'f1':[del_mol_ls_1],
                            'ffrac1':[del_mol_ls_frac1],'ffrac2':[del_mol_ls_frac2],
                            'ffrac5':[del_mol_ls_frac5],'p0':[del_mol_ls_per_0], 
                            'p1':[del_mol_ls_per_1],'pfrac1':[del_mol_ls_per_frac1],
                            'pfrac2':[del_mol_ls_per_frac2],
                            'pfrac5':[del_mol_ls_per_frac5],'conf_no':[len_del_mol_ls] })
    
    
    return result_df


def get_just_norn(p_ls,csv_df,ideal_dict, T=298.15):
    
    ## just calculate n0_molfac
    
    ## go through from reopt conf = 1 to n 
    
    sel_df_ls=[]
    for single_p_ls in p_ls:
        sel_df=get_delG_df(single_p_ls,csv_df,ideal_dict,temperature=T)
        sel_df_ls.append(sel_df)
        
    del_mol_ls=[1]
    
    for idx in range(1,len(sel_df_ls)):
        
        ## take into the consideration of conformers with repeat structure 
        
        try:
            del_idx=list(set(sel_df_ls[idx]['idx'].tolist())-set(sel_df_ls[idx-1]['idx'].tolist()))[0]
            
            curr_df=sel_df_ls[idx]
            mol_frac=curr_df[curr_df['idx']==del_idx]['mol_fraction'].tolist()[0]
            
        except IndexError:
            mol_frac=0
            
        del_mol_ls.append(mol_frac)
        
    
    _,n0_molfac = get_n0(del_mol_ls, cutoff=1)
    
    result_df=pd.DataFrame({'n0_molfac':[n0_molfac]})
        
    return result_df


def get_ls_given_priority_ls(priority_ls):
    
    
    analy_ls=[priority_ls]

    for i in range(1,len(priority_ls)):
        analy_ls.append(priority_ls[:-i])
    
    analy_ls.reverse()
    
    return analy_ls

##########
## execution of the analyses 
    
class MolFraTest_ml:
    
    def __init__(self, rmsCheck_df, result_df, delG_df):
        
        ## result_df = priority_ls 
        ## load data
        
        self.rmsCheck_df=rmsCheck_df
        self.result_df=result_df
        self.delG_df_ls=[delG_df]
        
        self.molname_ls=[i[0] for i in list_duplicates(self.result_df['name'].tolist())]
        self.conf_no_ls=[len(i[1]) for i in list_duplicates(self.result_df['name'].tolist())]
        
        
        
    def NewMolfrac_df(self, just_norn='No', temp=298.15):
        
        
        ## take out the required information from the data frames 
        
        dict_ls=[self.rmsCheck_df[self.rmsCheck_df['name']==mol]['dict'][0] for mol in self.molname_ls]
        
        
        other_select_ls=[]
        for mol in self.molname_ls:
            mol_other_df=self.result_df[self.result_df['name']==mol]

            other_ls = [x for x in mol_other_df['priority_ls'].tolist()][0]
                
            other_select_ls.append(get_ls_given_priority_ls(other_ls))
    
        other_df_ls=[]
        
        for idx in range(0,len(self.molname_ls)):
            

            if just_norn=='No':
                
                ## calculate all the descriptors 
                other_df=get_new_mol_data(other_select_ls[idx],self.delG_df_ls[idx],dict_ls[idx], temp)
            
            else:
                ## only calculate n0_molfac
                other_df=get_just_norn(other_select_ls[idx],self.delG_df_ls[idx],dict_ls[idx], temp)
                
            other_df['name'] = self.molname_ls[idx]
            other_df_ls.append(other_df)
            
        
        #ideal_df=pd.concat(ideal_df_ls)
        self.other_df=pd.concat(other_df_ls)
    
    
####################


def get_descriptor(rmsCheck_df, priority_df, delG_df, temp=298.15, pls_len=None):
    
    ## using class MolFraTest_ml to generate desriptor df for a priority df 
    ## only applicable for one molecule at the time 
    
    ## extract the index of the conformer based on their name in the delG_df
    
    delG_df['opt'] = [len(i.split('_')) for i in delG_df['name'].tolist()]
    
    opt_len_ls = delG_df['opt'].unique().tolist()
    
    reopt_idx_ls=[]
    for i in delG_df['name'].tolist():
        if len(opt_len_ls)>1:
            len_no=len(i.split('_'))
            if len_no == np.amax(opt_len_ls):
                reopt_idx_ls.append(int(i.split('_')[-2]))
            else:
                reopt_idx_ls.append(int(i.split('_')[-1]))
                
        else:
            reopt_idx_ls.append(int(i.split('_')[-1]))
            
    
    reopt_idx_ls1=[i-1 for i in reopt_idx_ls]
    ## priority_df - there should only be one line of data 
    priority_ls = priority_df['priority_ls'][0]

    reopt_idx_ls_sort = [i for i in priority_ls if i in reopt_idx_ls1][:pls_len]
    
    priority_df_reopt=priority_df.copy()
    priority_df_reopt=priority_df_reopt[['idx']]
    priority_df_reopt['priority_ls'] = [reopt_idx_ls_sort]
    priority_df_reopt['nor_nc'] = [len(reopt_idx_ls)/len(priority_ls)]
    
    name_ls=[n[:-4] for n in priority_df['name'].tolist()]
    priority_df_reopt['name'] = name_ls
    
    
    ## test 2 MolFraTest_ml on the incomplete priority ls 
    test2=MolFraTest_ml(rmsCheck_df, priority_df_reopt, delG_df)
    test2.NewMolfrac_df(temp=temp)
    
    
    descriptor_df1=pd.merge(test2.other_df, priority_df_reopt,on = 'name')
    descriptor_df2=descriptor_df1[['n0_molfac','f0','f1','ffrac1','ffrac2','ffrac5','p0','p1',
                               'pfrac1','pfrac2','pfrac5','nor_nc','conf_no']]

    
    
    descriptor_ary=descriptor_df2.to_numpy()
    
    return descriptor_ary, reopt_idx_ls1



        





