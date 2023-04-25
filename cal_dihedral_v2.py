#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 10:43:30 2022

@author: chingchinglam

update from version 1:
    rewrite dihedral_angle function to aviod the use of Bio.PDB.vectors and math package 
    which should also improve the cpu time consumption 

"""


import numpy as np
import pandas as pd

#from math import pi
#from Bio.PDB.vectors import Vector, calc_dihedral

def parse_line(number, string):
    '''
    extra xyz coor info from mol file

    Parameters
    ----------
    number : int
        e.g. 1
    string : str
        e.g. '0.090001 1.000401 2.000031  C'

    Returns
    -------
    lista : list
        [1, 'C', [0.090001, 1.000401, 2.000031]]

    '''
    
    stringa=string.split() 
    
    for index in range(0, 3):
        stringa[index] = float(stringa[index])

    lista = [number]+[stringa[-1]]+[stringa[0:3]]
    
    return lista


def dihedral_angle(list_A1, list_A2, list_B1, list_B2):
    '''
    calculate the dihedral angle

    Parameters
    ----------
    list_A1, list_A2, list_B1, list_B2 : list
        [0.090001, 1.000401, 2.000031] - xyz coor

    Returns
    -------
    di_angle : float

    '''
    A1 = np.array(list_A1)
    A2 = np.array(list_A2)  
    B1 = np.array(list_B1) 
    B2 = np.array(list_B2)   
    vec_a = A2 - A1 
    vec_b = B1 - A2  
    vec_c = B2 - B1 

    m = np.cross(vec_a,vec_b)
    n = np.cross(vec_b,vec_c)
    psi = np.sign(np.dot(vec_a,n))*np.arccos(np.dot(n,m)/( np.sqrt(np.dot(n,n)) *np.sqrt(np.dot(m,m)) ) )
    di_angle=psi*180/np.pi
    
    return di_angle


def dihedral_descriptor(dihedral_list, structure):
    '''
    

    Parameters
    ----------
    dihedral_list : list
        [['C 1','C 2','C 3', 'C 4'],[...]]
    structure : TYPE
        [[1, 'C', [0.090001, 1.000401, 2.000031]],
         [2, 'C', [1.090301, 5.004901, 1.030031]],
         ...]

    Returns
    -------
    dihedral_conf : list
        [134.25, 176.34, ...]

    '''
    
    
    ## extract the index from the str(sym+index) 
    dihedral_list_idx=[]
    for dihedral in dihedral_list:
        idx=[int(dihedral[no].split(' ')[1]) for no in range(0,4)]
        dihedral_list_idx.append(idx)
    #print(dihedral_list_idx)
    
    ## extract the xyz cooridinate based on index
    dihedral_list_xyz=[]
    for idx in dihedral_list_idx:
        cor_idx=[ a-1 for a in idx]
        ## index-1 from str(sym+index) to get index of the atom in the sdf/mol
        xyz_coor=[ structure[cor_idx[n]][2] for n in range(0,4)]
        dihedral_list_xyz.append(xyz_coor)
    #print(dihedral_list_xyz)

    ## calculate the all the dihedral angle value for all the dihedrals in the list 
    ## and compile into a list 
    dihedral_conf=[dihedral_angle(dihedral_list_xyz[l][0], dihedral_list_xyz[l][1], 
               dihedral_list_xyz[l][2], dihedral_list_xyz[l][3]) for l in range(0, len(dihedral_list_xyz))]
    
    #print(dihedral_list_xyz)
    
    return dihedral_conf





def dihedral_df(sdf_file, dihedral_list, bond_list):
    '''
    generate dihedral angle df
    
    Parameters
    ----------
    sdf_file : str
        e.g. 'test_01.sdf'
    dihedral_list : list
        [['C 1','C 2','C 3', 'C 4'],[...]]
    bond_list : list
        list of bonds [['C 1','C 2'],['C 3', 'C 4'],[...]]

    Returns
    -------
    df : dataframe 
         table of dihedral angle values with the rotational bond as the column index
         
    '''
    
    ## convert the bond list format to ['C1_C2', 'C3_C4', ...]

    
    ## convert the bond list format to ['C 1_C 2', 'C 3_C 4', ...]
    bond_list2=[bond_list[o][0]+'_'+bond_list[o][1] for o in range(0,len(bond_list))]
    #print(bond_list2)
    
    
    ## extract the xyz coor of mol in the sdf file
    sdf_output=[]
    
    with open(sdf_file) as output:
        for line in output:
            sdf_output+= [line]
    END_index_list =[k for k,d in enumerate(sdf_output) if d=='M  END\n']
    #print(END_index_list)
    
    #try:
        #atom_no=int(sdf_output[3].split(' ')[1])
        
        
    #except ValueError:
        
    atom_no=int(sdf_output[3][:3])
    
    #print(atom_no)
        
    #bond_no=int(sdf_output[3].split(' ')[2])
    
    len_mol=END_index_list[0]

    division1=[]
    for j in END_index_list:
        indv_mol=sdf_output[j-len_mol:j+1]
        indv_mol_atom1=indv_mol[4:4+atom_no]
        indv_mol_atom2=[j[0:35] for j in indv_mol_atom1]
        indv_mol_atom3=[ parse_line(idx, indv_mol_atom2[idx-1]) for idx in range(1,len(indv_mol_atom2)+1)]
    
        division1.append(indv_mol_atom3)

    ## calculate the dihedral angle values 
    dihedral_descriptor_list=[]
    for conformer in division1:
        dihedral_descriptor_list.append(dihedral_descriptor(dihedral_list, conformer))
    #print()
    
    ##compile into df
    dihedral_descriptor_array=np.array(dihedral_descriptor_list)
    #print(dihedral_descriptor_array)
    df = pd.DataFrame(dihedral_descriptor_array, columns=bond_list2)
    
    return df 
