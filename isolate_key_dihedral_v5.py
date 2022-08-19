#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 11:04:04 2022

@author: chingchinglam

Update from v1:
    change in the remove_FixedBond criteria 
    
update from v2: 
    include termindal dihedrals that may be linked to potential H bonding 
    
update from v3: 
    Remove C-H (CH3) and Si-H (SiH3) from termindal dihedrals 
    Reintroduce remove_terminalCCX3(df_e, bond_list) function into the pipeline
    Update remove_terminalCCX3 function: only take away C-CX3
    Update the compose_dihedral functions 

update from v4: 
    update functions to aviod using scm.plams and openbabel package 
    remove functions that are no longer in use in this script


"""

######################## import packages 
import pandas as pd
from rdkit import Chem
from collections import defaultdict


######################## supporting functions

def remove(duplicate): 
    ## remove duplicates in a list 
    
    final_list = [] 
    for num in duplicate: 
        if num not in final_list: 
            final_list.append(num) 
    return final_list 


def list_duplicates(seq):
    ## from a list of items to [item, [list of indexes at which the item appears in the list]]
    ## e.g. input: [1, 2, 3, 4, 4, 3]
    ## output: [[1, [0]], [2, [1]], [3, [2, 5]], [4, [3, 4]]]

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,len(locs)] for key,locs in tally.items() if len(locs)>0]


######## generate bond list and bond order matrix


def bond_info_array(mol):
    ## rdkit mol
    ## give a list of bond info incl. atom index and bond order
    bonds_info = [[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()] for bond in mol.GetBonds()]
    return bonds_info 

def atoms_info_array(mol):
    ## rdkit mol
    ## give a list of atom info incl. atom symbol, atomic number and if atom is in a ring (True or False)
    atoms_info = [ [ atom.GetSymbol(), atom.GetAtomicNum(), atom.IsInRing()] for atom in mol.GetAtoms()]
    return atoms_info


def get_atomlist(mol_rdkit_H):
    '''
    Parameters
    ----------
     mol_rdkit_H: Chem.MolFromMolFile(mol_file, removeHs=False)


    Returns
    -------
    atom_list, noH_atom_list : list
        [[1, 'C', 'C 1', [-0.2276, -0.2958, -1.3124], 0],
         [2, 'C', 'C 2', [0.4749, 0.4281, -0.2341], 1], ...]
    '''
    
    ## mol = rdkit mol
    ## generate atom list (atom symbol, index, str(symbol+index) and [xyz])
    coor_ls=mol_rdkit_H.GetConformer().GetPositions().tolist()
    atom_list =[]
    for atom,ls in zip(mol_rdkit_H.GetAtoms(), coor_ls):
        atom_list.append([atom.GetIdx()+1, atom.GetSymbol(), atom.GetSymbol()+' '+str(atom.GetIdx()+1), ls])
        
    ## generate a atom list without H atoms 
    ## updated to (atom symbol, index, str(symbol+index), [xyz], index in the list without H)
    
    noH_atom_list=[]
    for k in atom_list:
        if k[1] != 'H':
            noH_atom_list.append(k)
        
    for rd_idx in range(0,len(noH_atom_list)):
        noH_atom_list[rd_idx].append(rd_idx)
        
    ## return the atom list and atom list wihtout H
    
    return atom_list, noH_atom_list
    
#dict_index={noH_atom_list[idx][4]:noH_atom_list[idx][2] for idx in range(0,len(noH_atom_list))}
#inv_dict_index={noH_atom_list[idx][2]:noH_atom_list[idx][4] for idx in range(0,len(noH_atom_list))}


def get_BondOrderMatrix(mol_rdkit, mol_rdkit_H):
    '''

    Parameters
    ----------
    mol_rdkit : rdkit.Chem.rdchem.Mol
    mol_rdkit_H: Chem.MolFromMolFile(mol_file, removeHs=False)
       

    Returns
    -------
    df_e : dataframe 
        bond order matrix (no heavy atom)
    bond_list : list
        list of bonds in the mol (no heavy atom)
    dict_index : dictionary
        {0: 'C 1', 1: 'C 2' }; {index_without_H : str(sym+index(with_H))}

    '''
    
    am = Chem.GetAdjacencyMatrix(mol_rdkit)
    df = pd.DataFrame(am)
    
    ## float all the element in df
    for col in df.columns: 
        df[col] = df[col].astype('float', errors='ignore')
        
    ## replace the element (1,0) in the df to bond order with bond info list 
    change_list=bond_info_array(mol_rdkit)
    for i in range(0,len(change_list)):
        df[change_list[i][0]][change_list[i][1]]= change_list[i][2]
        df[change_list[i][1]][change_list[i][0]]= change_list[i][2]
    
  
    ## use get_atomlist function to generate noH_atom_list
    ## convert the index and column label of the df to str(symbol+index); 
    ## index = gview index / index in the mol or xyz file with H 
    noH_atom_list=get_atomlist(mol_rdkit_H)[1]
    dict_index={noH_atom_list[idx][4]:noH_atom_list[idx][2] for idx in range(0,len(noH_atom_list))}
    ## element = [noH_atom_list[i][2] for i in range(0,len(noH_atom_list))]
    
    df_e = df.rename(index=dict_index,  columns=dict_index)
    
    ## assembling the bond list without H from bond info
    ## [atom1 str(symbol+index), atom2 str(symbol+index), float(bond order)]
    ## atom to atom mapping to accommodate index differences of heave atoms in list with and without H
    
    begin_atom_list=[item[0] for item in change_list]
    end_atom_list=[item[1] for item in change_list]
    bond_order_list=[item[2] for item in change_list]
    begin_atom_list_idx = list((pd.Series(begin_atom_list)).map(dict_index))
    end_atom_list_idx = list((pd.Series(end_atom_list)).map(dict_index))

    bond_list=[[begin_atom_list_idx[no],end_atom_list_idx[no],bond_order_list[no]] for no in range(0,len(begin_atom_list))]


    return df_e, bond_list, dict_index
    


########## remove bond scripts 


def remove_FixedBond(mol_rdkit, dict_index, bond_list):
    ## to clean up the list: remove fixed bond (bond order > 1)-- sort1
    
    ###### to get bonds with bond order = 1 
    bond_list_sort1=[bond for bond in bond_list if bond[2] <= 1.0]
    
    ###### to get bonds with bond order = 1.5 but with heteroatoms (C-X or X-X) 
    
    bond_list_sort15=[bond for bond in bond_list if bond[2] == 1.5]
    bond_list_sort15_hetero=[bond for bond in bond_list_sort15 if bond[0][0] != 'C' or bond[1][0] != 'C']
    
    ###### to get X-X bonds (X= heteroatoms) with bond order = 2.0 
    
    bond_list_sort2=[bond for bond in bond_list if bond[2] == 2.0]
    bond_list_sort2h_hetero=[bond for bond in bond_list_sort2 if bond[0][0] != 'C' and bond[1][0] != 'C']
    
    ###### to get bonds with C=N (C=N must not be within a ring) with bond order = 2  
    ## get a list of ring N atoms and a list of ring C atoms 
    info_ring = mol_rdkit.GetRingInfo()
    ring_atom_list = list(info_ring.AtomRings())
    sort_ring_atom_list = [x for y in ring_atom_list for x in y]
    ring_atom_sym = list(map(dict_index.get, sort_ring_atom_list))
    ring_atom_sym_N = [i  for i in ring_atom_sym if i[0] == 'N']
    ring_atom_sym_C = [i  for i in ring_atom_sym if i[0] == 'C']
    
    ## get a list of C=N bonds
    #bond_list_sort2=[bond for bond in bond_list if bond[2] == 2.0]
    bond_list_sort2_NC = [bond for bond in bond_list_sort2 if bond[0][0] == 'N' and bond[1][0] == 'C'] + [
        bond for bond in bond_list_sort2 if bond[0][0] == 'C' and bond[1][0] == 'N']

    bond_list_sort2_NC_noring = []
    for b_ls in bond_list_sort2_NC:
        check=0
        for b in ring_atom_sym_N:
            if b in b_ls: 
                check += 1 
        for c in ring_atom_sym_C:
            if c in b_ls:
                check += 1 
        if check <= 1 :
            ## i.e. C(ring)=N, C=N(ring) and C=N are all ok; C(ring)=N(ring) is not ok  
            bond_list_sort2_NC_noring.append(b_ls)
    
    ## checking
    #print(ring_atom_sym_N)
    #print(ring_atom_sym_C)
    #print(bond_list_sort2_NC)
    #print(bond_list_sort2_NC_noring)
    
    
    ###### to get bonds with C-C(-X) with bond order = 1.5
    ## (X-O,S or N, C-X bond order = 2) 
    ## get a list of C atom within C=X bond, where X = O or N
    bond_withX = [bond for bond in bond_list_sort2 if bond[0][0] == 'O' or bond[1][0] == 'O' or bond[0][0] == 'N' or bond[1][0] == 'N']
    #bond_withX = [bond for bond in bond_list if bond[0][0] == 'O' or bond[1][0] == 'O' 
    #              or bond[0][0] == 'N' or bond[1][0] == 'N' or bond[0][0] == 'S' or bond[1][0] == 'S']
    bond_withX_at = [k[0:2] for k in bond_withX]
    bond_withX_at2 = [x for y in bond_withX_at for x in y]    
    bond_CInCX= [at for at in bond_withX_at2 if at[0]=='C']
    
    
    bond_list_sort15_CCX =[]
    for bond in bond_list_sort15:
        check_CX = 0 
        for at in bond_CInCX:
            if at in bond: 
                check_CX += 1 
        if check_CX >= 1: 
            bond_list_sort15_CCX.append(bond)
    
    ## checking
    #print(bond_withX)
    #print(bond_withX_at)
    #print(bond_CInCX)
    #print(bond_list_sort15_CCX)
    
    result = remove(bond_list_sort1 + bond_list_sort15_hetero + bond_list_sort2_NC_noring + bond_list_sort15_CCX + bond_list_sort2h_hetero) 
    
    
    
    return result

def remove_CX(bond_list):
    #and bonds that cannot be extended to become a dihedral -- sort2 , i.e. C-X, X=F,Cl,Br,I 
    bond_list_sort2=[bond for bond in bond_list if bond[0][0] != 'F' and bond[0][0] != 'Cl'
                    and bond[0][0] != 'Br' and bond[0][0] != 'I'and bond[1][0] != 'F' and 
                     bond[1][0] != 'Cl' and bond[1][0] != 'Br' and bond[1][0] != 'I']
    return bond_list_sort2


def split_terminalCX(df_e, bond_list):
    # put terminal C-C bonds in a separate list - sort3 -- C-CH3; C-NH2; C-OH; C-SH
    # for each atom in the bond, find the number of neighbouring atoms; 
    # identify the bond if number of neighbouring atoms > 1 
    
    bond_list_sort3=[]
    a1_a2_ls=[]
    
    for b in bond_list:
        a1 = len(df_e[df_e[b[0]] != 0.0]) 
        a2 = len(df_e[df_e[b[1]] != 0.0]) 

        con = False
        if a1 == 1 or a2 == 1:
            con = True 
            a1_a2_ls.append(b)
            
        
        if con == False:
            bond_list_sort3.append(b)
    
    #print(a1_a2_ls)
            
    return bond_list_sort3, a1_a2_ls


def remove_terminalCX(df_e, bond_list):
    ## remove the true terminal CX bond, e.g. C=O, CN triple bond
    bond_list_sort=[]
    
    for b in bond_list:
        a1 = len(df_e[df_e[b[0]] != 0.0]) 
        a2 = len(df_e[df_e[b[1]] != 0.0]) 
    
        con = False
        if a1 == 1 or a2 == 1:
            con = True 
        
        if con == False:
            bond_list_sort.append(b)
    
    return bond_list_sort

def remove_CH3(df_e_complete, bond_list):
    ## remove CH3 and SiH3 bonds from the bond_list based on the complete bond matrix df 
    
    bond_list_sort=[]
    
    for b in bond_list:
        a1_df = df_e_complete[df_e_complete[b[0]] != 0.0]
        a2_df = df_e_complete[df_e_complete[b[1]] != 0.0]
    
        a1_ls=list(a1_df.index)
        a2_ls=list(a2_df.index)
        a1_ls_X = [a1_ls[x_no][0:2] for x_no in range(0, len(a1_ls))]
        a2_ls_X = [a2_ls[x_no][0:2] for x_no in range(0, len(a2_ls))]
        
        a1_ls_X_noH = [ele for ele in a1_ls_X if ele == 'H ']
        a2_ls_X_noH = [ele for ele in a2_ls_X if ele == 'H ']
        
        
        toremove = False 
        if len(a1_ls_X_noH) == 3 and ( b[0][0] == 'C' or b[0][0:2] == 'Si' ): 
            toremove = True
        elif len(a2_ls_X_noH) == 3 and ( b[1][0] == 'C' or b[1][0:2] == 'Si' ): 
            toremove = True 
        
        if toremove == False:
            bond_list_sort.append(b)
        
        
    
    return bond_list_sort



def remove_3MemRing(mol_rdkit, dict_index, bond_list):
    
    ## remove the bonds within the three membered ring 
    info_ring = [ atom.IsInRingSize(3) for atom in mol_rdkit.GetAtoms()]
    
    ring_dict={True:1, False:-1}
    ring_info_no=list(map(ring_dict.get, info_ring ))
    sym_label=list(dict_index.values())

    sum_ring_no_dict={sym_label[idx]:ring_info_no[idx] for idx in range(0,len(ring_info_no))}
    
    
    ## with the above, bond within a 3-mem ring = 2 (selected); 
    ## bond connected to a 3-mem ring = 0 (not selected); 
    ## other bonds = -2 (not selected)
    
    new_bond_list=[]
    for bond in bond_list:
        bond_point=sum_ring_no_dict.get(bond[0])+sum_ring_no_dict.get(bond[1])
        if bond_point < 1: 
            new_bond_list.append(bond)
    
            
    return new_bond_list


def remove_terminalCCX3(df_e, bond_list):
    
    ## remove terminal C-CX3 bonds - sort4 -- C-CCl3
    ## for each atom in the bond, find the number of neighbouring; 
    ## collect the str(symbol+index) of the neighbours

    bond_list_sort4=[]
    for b in bond_list:
        a1_df =df_e[df_e[b[0]] != 0.0]
        a2_df =df_e[df_e[b[1]] != 0.0]
        
        a1_ls=list(a1_df.index)
        a2_ls=list(a2_df.index)
        a1_ls_X = [a1_ls[x_no][0:2] for x_no in range(0, len(a1_ls))]
        a2_ls_X = [a2_ls[x_no][0:2] for x_no in range(0, len(a2_ls))]
        #print(a2_ls_X)
        #print(a1_ls_X)

        # retain only X atoms in the list
        wanted = {'F ','Cl','I ','Br'}
        a1_ls_X_rm = [ele for ele in a1_ls_X if ele in wanted]
        a2_ls_X_rm = [ele for ele in a2_ls_X if ele in wanted]
        
        wanted_ls=[['F ','F ','F '],['Cl','Cl','Cl'], ['I ','I ','I '],['Br','Br','Br']]
        
        # remove bond if there are 3 X atoms in the list of neighbours

        X_con=False
        if a1_ls_X_rm in wanted_ls or a2_ls_X_rm in wanted_ls:
            X_con=True 
            #print(X_con)
    
        if X_con == False:
            bond_list_sort4.append(b)
            
    #print(bond_list_sort4)
    
    return bond_list_sort4 


###############################################################

########################
########################


def select_atom(df_e,df_e_complete,bond_list, mol, dict_index):
    '''
    Parameters
    ----------

    mol : rdkit.Chem.rdchem.Mol
    df_e : dataframe 
        bond order matrix (no heavy atom)
    bond_list : list
        list of bonds in the mol (no heavy atom)
    dict_index : dictionary
        {0: 'C 1', 1: 'C 2' }; {index_without_H : str(sym+index(with_H))}

    Returns
    -------
    final_bond_list : list
        proccessed bond list 

    '''
    
    
    ## mol = rdkit mol
    ## to clean up the list: remove fixed bond (bond order > 1)-- sort1 
    bond_list_sort1=remove_FixedBond(mol, dict_index,bond_list)
    
    ## and bonds that cannot be extended to become a dihedral -- sort2 , i.e. C-X, X=F,Cl,Br,I 
    bond_list_sort2=remove_CX(bond_list_sort1)
    
    ## retain terminal C-X bonds in a separate list -- C-CH3; C-NH2; C-OH; C-SH 
    split_bonds = split_terminalCX(df_e, bond_list_sort2)
   
    teriminal_CX = split_bonds[1]
    
    teriminal_CX2 = remove_terminalCX(df_e_complete, teriminal_CX )
    teriminal_CX3 = remove_CH3(df_e_complete, teriminal_CX2)
    
    final_bond_list2= [t[0:2] for t in teriminal_CX3]
    
    ## remove terminal C-X bonds in the main list -- C-CH3; C-NH2; C-OH; C-SH 
    bond_list_sort4 = split_bonds[0]
        
    ## remove terminal C-C bonds - sort3 -- C-CX3 X=F,Cl,Br,I     
    bond_list_sort5=remove_terminalCCX3(df_e, bond_list_sort4)
    
    ## remove 3-membered ring 
    bond_list_sort6=remove_3MemRing(mol, dict_index, bond_list_sort5)
    

    ## remove bonds that contain sp atom 
    #bond_list_sort5=remove_sp(mol, dict_index, bond_list_sort4)
    ## remove ring linkage bonds in fused/bridged ring systems 
    #bond_list_sort7=remove_fixring(mol, dict_index, df_e, bond_list_sort6)
    
    final_bond_list = [k[0:2] for k in bond_list_sort6]
    
    return final_bond_list, final_bond_list2



################### compose dihedral script


def CXH_prepare4_compose_dihedral(mol_rdkit_H):
    
    am = Chem.GetAdjacencyMatrix(mol_rdkit_H)
    df = pd.DataFrame(am)

    for col in df.columns: 
        df[col] = df[col].astype('float', errors='ignore')
        
    ## replace the element (1,0) in the df to bond order with bond info list 
    change_list=bond_info_array(mol_rdkit_H)
    for i in range(0,len(change_list)):
        df[change_list[i][0]][change_list[i][1]]= change_list[i][2]
        df[change_list[i][1]][change_list[i][0]]= change_list[i][2]
    
    ## get atom symbol with atom_info
    ## convert the index and column label of the df to str(symbol+index); 
    ## index = gview index / index in the mol or xyz file with H 

    atom_info=atoms_info_array(mol_rdkit_H)
    atom_sym_ls=[atom_info[i][0]+' '+str(i+1) for i in range(0,len(atom_info))]
    no_ls=list(range(0,len(atom_info)))
    dict_index_H={no_ls[idx]:atom_sym_ls[idx] for idx in range(0,len(atom_info))}
            
    df_e_complete = df.rename(index=dict_index_H,  columns=dict_index_H)
    
    return dict_index_H, df_e_complete



def atom_within_same_ring(mol_rdkit, dict_index, chosen_atom):
    
    ## use in the compose_dihedral function
    ## output: [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 
    ## each data point corresponds to an atom (1 --> within the same ring as the chosen atom; 0 --> not within the same ring)
    
    info_ring=mol_rdkit.GetRingInfo()
    ring_atom_list = list(info_ring.AtomRings())
    ring_atom_list_map=[list(map(dict_index.get, atom_ls)) for atom_ls in ring_atom_list]

    atom_in_this_list=[]
    for i in ring_atom_list_map:
        if chosen_atom in i:
            atom_in_this_list.append(ring_atom_list_map.index(i))

    atom_in_ring_point_list=[]
    for at in list(dict_index.values()):
        atom_in_ring_point=0
        for idx in atom_in_this_list:
            if at in ring_atom_list_map[idx]:
                atom_in_ring_point+=2
    
        atom_in_ring_point_list.append(atom_in_ring_point)
    
    #print(atom_in_ring_point_list)

    return atom_in_ring_point_list



def compose_dihedral(mol_rdkit, mol_rdkit_H, dict_index, df_e, final_bond_list, bond, anyH='No'):
    '''

    Parameters
    ----------
    mol_rdkit : rdkit.Chem.rdchem.Mol
    
    df_e : dataframe 
        bond order matrix (no heavy atom)
    final_bond_list : list
        shortlisted bonds
    bond : i in a for loop
        

    Returns
    -------
    dihedral : list
        shortlisted dihedrals
    
    '''
    
    # this function is written for one bond with two atoms 
    # composing the dihedral angles - part 1 - identify lists of neighbouring atoms 
    # find neighbouring atoms from the final_bond_list

    mol=mol_rdkit
    noH_atom_list=get_atomlist(mol_rdkit_H)[1]
    
    if anyH == 'Yes':
        noH_atom_list=get_atomlist(mol_rdkit_H)[0]
        
    #print(noH_atom_list)
    
    
    a1_df =df_e[df_e[bond[0]] != 0.0]
    a2_df =df_e[df_e[bond[1]] != 0.0]
    
    a1_ls=list(a1_df.index)
    a2_ls=list(a2_df.index)
    
    
    final_bond_atom_ls=[x for y in final_bond_list for x in y]
    
    ## compile neighbouring atoms that are also within the final_bond_atom_ls
    ele_0_simpl=[at for at in a1_ls if at in final_bond_atom_ls and at != bond[1]]
    ele_1_simpl=[at for at in a2_ls if at in final_bond_atom_ls and at != bond[0]]
    
    

    # when there's no neighouring atom in the final_bond_list, 
    # find neighouring atoms using the df_e bond order matrix

    if len(ele_0_simpl) == 0:
        # need to remove the atom in the bond from the list and 
        # combine sublist within the list
        
        ele_0_simpl.append(list(set(a1_ls)-set([bond[1]])))
        ele_0_simpl=[x for y in ele_0_simpl for x in y]   

    if len(ele_1_simpl) == 0:
        ele_1_simpl.append(list(set(a2_ls)-set([bond[0]])))
        ele_1_simpl=[x for y in ele_1_simpl for x in y]
    
    ## same_ele_0_simpl and same_ele_1_simpl are the final lists
    ## composing the dihedral angles - part 2 - construct the dihedral angle list 

    atom_info=atoms_info_array(mol)

    ## if the atom is not a C or H -- +2
    atom_point_list=[]
    for at in atom_info: 
        atom_point = 0
        if at[1] != 6 and at[1] != 1: 
            atom_point += 2
        if at[1] == 6:
            atom_point += 1
            
    
        atom_point_list.append(atom_point)
    
    ## the list of points for atom with in the same ring as atom1 or atom0
    atom0=atom_within_same_ring(mol_rdkit, dict_index, bond[0])
    atom1=atom_within_same_ring(mol_rdkit, dict_index, bond[1])
    
    
    ## atom-to-atom mapping to assign the corresponding point to the shortlisted atoms
    dict_atom_point_0={noH_atom_list[idx][2]: atom_point_list[idx]+atom0[idx] for idx in range(0,len(noH_atom_list)) }
    dict_atom_point_1={noH_atom_list[idx][2]: atom_point_list[idx]+atom1[idx] for idx in range(0,len(noH_atom_list)) }
    
    ele_0_point = list(map(dict_atom_point_0.get, ele_0_simpl))
    ele_1_point = list(map(dict_atom_point_1.get, ele_1_simpl))
    
    
    #print(bond)
    #print(ele_0_point)
    #print(ele_1_point)

    ## the atom with the highest point to be selected 
    ele_0_chosen=ele_0_simpl[ele_0_point.index(max(ele_0_point))]
    ele_1_chosen=ele_1_simpl[ele_1_point.index(max(ele_1_point))]

    ## compose the list of atoms in the dihedral
    dihedral = [ele_0_chosen]+bond+[ele_1_chosen]
    #print(dihedral)
    
    return dihedral



##################### embedded functions above 

def isolate_dihedral(sdf_file):
    '''
    Parameters
    ----------
    sdf_file : str -- path to the sdf file

    Returns
    -------
    dihedral_list : list -- [['O 1', 'C 2', 'N 15', 'C 16'], ...]
    final_bond_list : list -- [['C 2', 'N 15'], ...]
    dict_index_list : dict -- {0: 'O 1', 1: 'C 2', 2: 'N 3', ...}

    '''
    
    ## combine all the process starting from a sdf file
    ## get rdkit mol without H atoms
    suppl = Chem.SDMolSupplier(sdf_file)
    mol_rdkit =  [mol for mol in suppl][0]
    
    ## get rdkit mol with H atoms
    suppl_H = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol_rdkit_H =  [mol for mol in suppl_H][0]
    
    dihedral_list, final_bond_list, dict_index_list = isolate_dihedral_mol(mol_rdkit, mol_rdkit_H)
    
    
    return dihedral_list, final_bond_list, dict_index_list


def isolate_dihedral_mol(mol_rdkit, mol_rdkit_H):
    
    ## combine all the process from mol_rdkit
    dict_index_H,df_e_complete =CXH_prepare4_compose_dihedral(mol_rdkit_H)
    df_e,bond_list,dict_index_list =get_BondOrderMatrix(mol_rdkit, mol_rdkit_H)
    final_bond_list1, final_bond_list2 = select_atom(df_e, df_e_complete, bond_list, mol_rdkit, dict_index_list)
    
    
    final_bond_list = final_bond_list1 + final_bond_list2
    
    ## compose dihedral angles 
    dihedral_list = [compose_dihedral(mol_rdkit, mol_rdkit_H, dict_index_list, df_e, final_bond_list, b) for b in final_bond_list1]
    dihedral_list2 = [compose_dihedral(mol_rdkit_H, mol_rdkit_H, dict_index_H, df_e_complete, final_bond_list, b,'Yes') for b in final_bond_list2]
    
    dihedral_list += dihedral_list2
    
    return dihedral_list, final_bond_list, dict_index_list








