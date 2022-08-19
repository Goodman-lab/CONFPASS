#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:22:42 2022

@author: chingchinglam
"""

import numpy as np
import pandas as pd

def find_gaps(numbers):
    
    ## find the difference between sequential value in a list 
    
    adjacent_differences = [(y - x) for (x, y) in zip(numbers[:-1], numbers[1:])]
    max_gap=np.amax(adjacent_differences)
    max_idx=adjacent_differences.index(max_gap)
    
    return max_idx,  max_gap, numbers[max_idx]


def correct_dih(dih_ls):
    
    ## given a list of dihedral angle values 
    
    ## Covert the θ to between a range of 0-360 degree
    shift_dih_ls=[i+180 for i in dih_ls]
    
    ## Find the maximum gap in data between 0-100 degree. 
    ## If this gap is greater than or equal to 15 degree, 
    ## move the data point from the left side of the gap by adding 360

    try:
        s_shift_dih_ls=[i for i in shift_dih_ls if i < 100]
        s_shift_dih_ls+=[100]
        ss_shift_dih_ls=sorted(s_shift_dih_ls)

        max_idx, max_gap, max_no= find_gaps(ss_shift_dih_ls)

        if max_gap >= 15:

            shift_dih_ls2=[]
            for i in shift_dih_ls:
                if i <= max_no:
                    shift_dih_ls2.append(i+360)
        
                else:
                    shift_dih_ls2.append(i)
        
            del_shift_dih_ls2 = [i-shift_dih_ls2[0] for i in shift_dih_ls2]
    
    
        else:
            del_shift_dih_ls2 = [i-shift_dih_ls[0] for i in shift_dih_ls]
        
        
    except ValueError:
        
        del_shift_dih_ls2 = [i-shift_dih_ls[0] for i in shift_dih_ls]
               
    
    return del_shift_dih_ls2


def shift_axis_df(df):
    
    ## shift the axis so that there are no negative values
    
    col_ls=[col for col in df.columns]
    
    nor_ls=[]
    
    for c in col_ls:
        col_as_ls = df[c].tolist()
        col_min= np.amin(col_as_ls)
        nor_col_as_ls=[]
        for v in col_as_ls:
            
            nor_col_as_ls.append(v-col_min)
        
        nor_ls.append(nor_col_as_ls)
    
    nor_dict = {col_ls[idx]:nor_ls[idx] for idx in range(0, len(col_ls))}
    nor_df = pd.DataFrame(nor_dict)
    
    return nor_df 

def normalise_dih_df(df):

    ## normalise the value in each column
    
    col_ls=[col for col in df.columns]
    
    nor_ls=[]
    
    for c in col_ls:
        col_as_ls = df[c].tolist()
        col_max = np.amax(col_as_ls)
        nor_col_as_ls=[]
        for v in col_as_ls:
            
            nor_col_as_ls.append(v/col_max)
        
        nor_ls.append(nor_col_as_ls)
    
    nor_dict = {col_ls[idx]:nor_ls[idx] for idx in range(0, len(col_ls))}
    nor_df = pd.DataFrame(nor_dict)
    
    return nor_df



def correction_by_gap(df):
    
    ## the execution function 
    
    dih_at_ls=list(df.columns)
    corr_dih_ls_ls=[]

    
    for i in dih_at_ls:
        dih_ls=df[i].tolist()
        corr_dih_ls_ls.append(correct_dih(dih_ls))

    refined_dih_df=pd.DataFrame({dih_at_ls[idx]:corr_dih_ls_ls[idx] for idx in range(0,len(dih_at_ls))})

    ## Find the relative θ compared to the θ of the global minimum at the MMFF level 
    ## Shift the axis so that there are no negative values; θ’ 

    refined_dih_df2=shift_axis_df(refined_dih_df)
    
    ## normalisation 
    refined_dih_df3=normalise_dih_df(refined_dih_df2)
    
    return refined_dih_df2, refined_dih_df3





