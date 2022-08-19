#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:20:11 2021

@author: chingchinglam

#update from v2:
#    introduce functions that process the dihedral angle df 

"""

import numpy as np
import pandas as pd


def dihedral_para(df):
    '''

    Parameters
    ----------
    df : dataframe 
         table of dihedral angle values with the rotational bond as the column index 

    Returns
    -------
    pd_outcome : dataframe 
        statistical parameters for each set of dihedral angle values 

    '''
    
    ## convert the df to array format 
    dihedral_ar=df.to_numpy()
    
    ## get number of column in the df and list the title of the column (the selected bonds)
    column_len=len(dihedral_ar[0,:])
    bonds=list(df.columns)

    ## use numpy to get statistical parameters 
    stand_dev=[]
    mean_va=[]
    median_va=[]
    max_va=[]
    min_va=[]
    range_va=[]
    for i in range(0,column_len):
        stand_dev.append(np.std(dihedral_ar[:,i]))
        mean_va.append(np.mean(dihedral_ar[:,i]))
        median_va.append(np.median(dihedral_ar[:,i]))
        max_va.append(np.amax(dihedral_ar[:,i]))
        min_va.append(np.amin(dihedral_ar[:,i]))
        range_va.append(np.amax(dihedral_ar[:,i])-np.amin(dihedral_ar[:,i]))
        
    ## compile the parameters into df
    dict_outcome={'bond':bonds,'stand_dev':stand_dev,'mean':mean_va,'median':median_va,
              'max':max_va, 'min':min_va, 'range':range_va}

    pd_outcome=pd.DataFrame(dict_outcome)

    return pd_outcome


def dihedral_para_ab(df):
    
    
    ## absoluate all values in the df 
    ab_df=np.absolute(df)
    
    return dihedral_para(ab_df)

def gen_merge_df(df):
    
    ## calculate the statistical parameters
    ## compile the relvant parameters into a df
    
    noab_df=dihedral_para(df)
    ab_df=dihedral_para_ab(df)
    noab_df2=noab_df[['bond','stand_dev','range']]
    ab_df2=ab_df[['bond','stand_dev','range']]
    ab_df3=ab_df2.rename(columns={'stand_dev':'stand_dev_ab', 'range':'range_ab'})

    merge_df= pd.merge(noab_df2, ab_df3, on='bond') 
    
    return merge_df
    
    


def list_fixed_bond(df):
    '''
    Parameters
    ----------
    df : dataframe 
         table of dihedral angle values with the rotational bond as the column index

    Returns
    -------
    remove_bond_format : list
        list of fixed bonds - to be removed 
    remove_bond : list
        list of fixed bonds - to be removed - input for remove_fixed_bond_df
    merge_df : dataframe
        df of statistical parameters 

    '''
    
    ## calculate the statistical parameters
    ## compile the relvant parameters into a df
    noab_df=dihedral_para(df)
    ab_df=dihedral_para_ab(df)
    noab_df2=noab_df[['bond','stand_dev','range']]
    ab_df2=ab_df[['bond','stand_dev','range']]
    ab_df3=ab_df2.rename(columns={'stand_dev':'stand_dev_ab', 'range':'range_ab'})

    merge_df= pd.merge(noab_df2, ab_df3, on='bond') 
    array_merge_df=merge_df.to_numpy()

    remove_bond=[]
    for i in array_merge_df:
        ## if the standard deviation is less than 5 and range is less than 10 (select)
        if i[1]<2.4 and i[2]<8.7:
            remove_bond.append(i[0])
        ## else if the range is greater than 355 and the standard deviation and the range of the 
        ## absoluate dihedral angle are less than 5 
        elif i[2]>359.2 and i[3]<1.5 and  i[4]<6.6:
            remove_bond.append(i[0])
        

    remove_bond_format=[[b.split('_')[0], b.split('_')[1]] for b in remove_bond]
    
    return remove_bond_format, remove_bond, merge_df


##################### embedded functions above 

def remove_fixed_dihedrals(df, dihedral_list, bond_list):
    '''
    
    Parameters
    ----------
    df : dataframe 
         table of dihedral angle values with the rotational bond as the column index
    dihedral_list : list
        list of dihedrals [['C 1','C 2','C 3', 'C 4'],[...]]
    bond_list : list
        list of bonds [['C 1','C 2'],['C 3', 'C 4'],[...]]

    Returns
    -------
    new_dihedral_list : list 
        list of dihedrals [['C 1','C 2','C 3', 'C 4'],[...]] 
        after removing the fixed dihedrals
    new_bond_list : list 
        list of bonds [['C 1','C 2'],['C 3', 'C 4'],[...]]
        after removing the fixed bonds
    remove_bond : list
        list of fixed bonds - to be removed - input for remove_fixed_bond_df

    '''
    
    ## conduct statistical parameter analyses
    fixed_bond=list_fixed_bond(df)
    remove_bond_format=fixed_bond[0]
    remove_bond=fixed_bond[1]
    #print(remove_bond_format)
    
    ## remove fixed bonds from the bond list 
    rm_bond_index=[]
    for b in remove_bond_format:
        if b in bond_list: 
            rm_bond_index.append(bond_list.index(b)) 
    
    
    ## remove the fixed dihedrals from the dihedral list 
    rm_dihedral=[]
    for idx in rm_bond_index: 
        rm_dihedral.append(dihedral_list[idx]) 
    
    
    ## generte the new bond and dihedral list 
    new_dihedral_list = [y for y in dihedral_list if y not in rm_dihedral]
    new_bond_list=[x for x in bond_list if x not in remove_bond_format]
    
    
    return new_dihedral_list, new_bond_list, remove_bond

def remove_fixed_bond_df(df,remove_label):
    '''

    Parameters
    ----------
    df : dataframe 
        table of dihedral angle values with the rotational bond as the column index
    remove_label : list
        list of fixed bonds - to be removed ['C 1_C 2', 'C 3_C 4', ...]

    Returns
    -------
    df_new : dataframe 
        table of dihedral angle values with the rotational bond as the column index
        the coloumns of dihedrals in the remove_label lists are removed

    '''
    
    
    leftover_label=list(set(list(df.columns)) - set(remove_label))
    ## final df (selected dihedral after the filter) 
    df_new=df[leftover_label]
    
    ## dihedrals filtered by the system
    df_remove=df[remove_label]
    
    
    
    return df_new, df_remove
    

