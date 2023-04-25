#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:33:52 2022

@author: chingchinglam


"""

import numpy as np
import pandas as pd
from collections import defaultdict

import random


def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]


##################
## random approach 


def random_ls_generator(length,total):
    
    ## generate a list with len(list)=length with random, non-repetitive number < total 
    
    random_ls=[]
    #for i in range(0,total+1):
    while len(random_ls) < length:
        #print(random_ls)
        random_no=random.randint(0,total-1)
        if random_no not in random_ls:
            random_ls.append(random_no)
    #else:
        #break
            
    return random_ls

def get_conf_ramdom_ls(conf_no):
    
    ## n_ls here refers to the index of the conformers 
    ## generate lists of randomly chosen conformer idx
    ## get_conf_ramdom_ls(3)
    ## Out[49]: [[1], [0, 1], [2, 0, 1]]
    
    n_ls=list(range(1,conf_no+1))

    randomlist_ls=[]
    for i in n_ls:
        #if i==0:
            #randomlist = [random_ls_generator(1,conf_no)[0]]
        #else:
        randomlist = random_ls_generator(i,conf_no) 
        randomlist_ls.append(randomlist)
        
    return randomlist_ls

##################
## every nth approach 


def get_nth_prority(conf_no, nth=2):
    
    ## Prioritise every nth conformer 
    ## e.g. [0,2,4,6,8,1,3,5,7,9]
    
    complete_ls=list(range(0,conf_no))
    
    prority_ls=[]
    
    for i in range(0,nth):
        prority_ls.append(complete_ls[i::nth])
        
    prority_ls1=[x for y in prority_ls for x in y]
        
    
    return prority_ls1

##################
## pipe descend

def get_pipe_de_ls(result_df, mol_name):
    
    ## result_df -- the df that contains clustering result 
    ## name of the molecule 
    ## choose a conformer from the cluster with the lowest delG at the FF level 
    
    result_mol=result_df[result_df['name']==mol_name]
    result_mol_ls=result_mol['clusters'].tolist()
    #print(result_mol_ls)
    
    de_pipe_select=[]
    for i in result_mol_ls:
        sub_ls=[]
        for j in i: 
            sub_ls.append(j[1][j[1].index(np.amin(j[1]))]) 
        
        de_pipe_select.append(sub_ls)
    
    return de_pipe_select


def get_prority_ls(p_ls):
    #print(p_ls)
    
    prority_ls=p_ls[0]
    
    for ls in p_ls[1:]:
        for i in ls:
            if i not in prority_ls:
                prority_ls.append(i)
                
    #print(prority_ls)
            
    return prority_ls

##################
## pipe x 


def get_format_result_df(result_df, x):
    
    ## given a result df, extract information and product a rmsCheck format df using result at nor_n = x
    
    mol_name_ls=[i[0] for i in list_duplicates(result_df['name'].tolist())]
    
    result_format_df_ls=[]
    for name in mol_name_ls:
        mol_df=result_df[result_df['name']==name]
        conf_no=len(mol_df)
        clus_no=round(conf_no*x)
        result_at_clus_no=mol_df[mol_df['n']==clus_no]
        cluster_ls = [x for x in result_at_clus_no['clusters']][0]
        cluster_ls2 = [ i[1] for i in cluster_ls]

        cluster_ls3=[]
        for i in cluster_ls2:
            cluster_ls3.append([[ j, cluster_ls2.index(i)] for j in i] )
    
        #cluster_ls4=sorted([x for y in cluster_ls3 for x in y])

        #cluster_dict={cluster_ls4[idx][0]:cluster_ls4[idx][1] for idx in range(0,len(cluster_ls4))}


        result_format_df_ls.append(pd.DataFrame({'name':[name], 'conf_no':[conf_no], 'clus_no':[clus_no],
                              'clusters':[cluster_ls2]}))
    


    result_format_df=pd.concat(result_format_df_ls)
    
    return result_format_df


def get_conf_x_ls(rmsCheck_df, mol_name):
    
    ## generate lists of conformers 
    ## The ideal scenario for selecting conformers for reoptimisations based on DFT result and RMSD check 
    
    ideal_mol=rmsCheck_df[rmsCheck_df['name']==mol_name]
    ideal_mol_ls=[x for x in ideal_mol['clusters'].tolist()]
    #print(ideal_mol_ls)
    
    ## start with assembling the complete list 
    rmsCheck_idx_ls=ideal_mol_ls[0]
    
    ## isolate the repetitive structures from the list of the list 
    ## the repetitive structures are put at the end of the list 
    rmsCheck_idx_ls_repeat=[i[1:] for i in rmsCheck_idx_ls]
    rmsCheck_idx_ls_repeat2=[x for y in rmsCheck_idx_ls_repeat for x in y]
    rmsCheck_idx_ls_single=[i[0] for i in rmsCheck_idx_ls]
    rmsCheck_idx_ls_sort_complete=rmsCheck_idx_ls_single+rmsCheck_idx_ls_repeat2
    
    complete_ideal_ls=[rmsCheck_idx_ls_sort_complete]

    for i in range(1,len(rmsCheck_idx_ls_sort_complete)):
        complete_ideal_ls.append(rmsCheck_idx_ls_sort_complete[:-i])
    
    complete_ideal_ls.reverse()
    
    return complete_ideal_ls


#################
## pipe x as -- default setting -- 20% of ls1(pipe x) + 80% of ls2 (pipe de)

def get_pip_mix(ls1,ls2, per=0.2):
    
    ## 20% of ls1 + 80% of ls2 
    
    conf_no=len(ls1)
    new_ls=ls1[0:round(conf_no*per)]
    
    for i in ls2:
        if i not in new_ls:
            new_ls.append(i)
            
    return new_ls



##################################
## implementation 


class GetPriority:
    
    def __init__(self, result_df):
        
        ## load data 
    
        self.result_df=result_df
        self.molname_ls=[i[0] for i in list_duplicates(self.result_df['name'].tolist())]
        self.conf_no_ls=[len(i[1]) for i in list_duplicates(self.result_df['name'].tolist())]

    

    def priority_df(self, x_=0.8, x_de_=0.2, n_=3, method_ = 'pipe_x_as'):
    
        ## generate priority list 
    
        idx_ls=list(range(0,len(self.molname_ls)))
        
        
        if method_ == 'pipe_as':
        
            pipe_de_select_ls=[get_prority_ls(get_pipe_de_ls(self.result_df, mol)) for mol in self.molname_ls]
            
            self.priority_ls_df=pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':pipe_de_select_ls})
            
            
        
        elif method_ == 'pipe_x':
            
            x_result_df = get_format_result_df(self.result_df, x_)
            pipe_x_ls=[get_conf_x_ls(x_result_df, mol)[-1] for mol in self.molname_ls]
            self.priority_ls_df=pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':pipe_x_ls})
        
        elif method_ == 'pipe_x_as':
            
            x_result_df = get_format_result_df(self.result_df, x_)
            
            pipe_x_ls=[get_conf_x_ls(x_result_df, mol)[-1] for mol in self.molname_ls]
            #print(pipe_x_ls)
            pipe_de_select_ls=[get_prority_ls(get_pipe_de_ls(self.result_df, mol)) for mol in self.molname_ls]
            
            
            pipe_x_de_ls=[get_pip_mix(ls1,ls2, per=x_de_) for ls1,ls2 in zip(pipe_x_ls,pipe_de_select_ls)]
            self.priority_ls_df=pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':pipe_x_de_ls})
            
            
        elif method_ == 'ascend':
        
            descend_select_ls=[get_nth_prority(no, nth=1) for no in self.conf_no_ls]
            self.priority_ls_df=pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':descend_select_ls})
        
        
        elif method_ == 'nth':
            
            nth_select_ls = [get_nth_prority(no, nth=n_) for no in self.conf_no_ls]
            self.priority_ls_df = pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':nth_select_ls})
            
            
        elif method_ == 'random':
            
            random_select_ls=[get_conf_ramdom_ls(no)[-1] for no in self.conf_no_ls]
            self.priority_ls_df=pd.DataFrame({'idx':idx_ls,'name':self.molname_ls, 'priority_ls':random_select_ls})
            
        
        else:
           
           print('method not available: please choose from pipe_de, pipe_x, pipe_x_as, ascend, nth or random')
            
    
    
    







