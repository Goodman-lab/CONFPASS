#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 09:47:50 2022

@author: chingchinglam

update from version 1:
    update functions to accommodate radicals 

"""

import os
import natsort
import pandas as pd
from rdkit import Chem
#import numpy as np 
#import math

def check(filename):
    
    ## open the output file 
    output_line=[]
    with open(filename) as output:
        for line in output:
            output_line.append(line)
        
    ## check if the structure has been optimised     
    Stationary_idx = 0
    for l in output_line: 
        if'    -- Stationary point found.' in l:
            Stationary_idx = output_line.index(l)
      
    
    #yesline_opt=output_line[Stationary_idx-6:Stationary_idx-2]
    #yesline_opt_yes=[i.split(' ')[-1] for i in yesline_opt]

    #yes_check=['YES\n', 'YES\n', 'YES\n', 'YES\n']
    #yesline_opt_test=True


    #if yesline_opt_yes != yes_check:
    #    yesline_opt_test=False
    
    opt_check=True
    ## use the following in if yesline_opt_test were to be considered
    #if Stationary_idx<1 or yesline_opt_test == False:
    if Stationary_idx<1:
        opt_check=False
        
    
    ## check the frequency calculation result     
    ind_Freq=0
    for l in output_line:
        if ' and normal coordinates:\n' in l:
            ind_Freq= output_line.index(l)

    Freq_useful=output_line[ind_Freq+3]
    Freq_sq=Freq_useful.split()
    
    try:
        Freq_no1=float(Freq_sq[2])
        Freq_no2=float(Freq_sq[3])
    except:
        Freq_no1=0
        Freq_no2=0
        
    
    
    
    freq_check=False

    if Freq_no1>0 and Freq_no2>0:
        freq_check=True

    
    ## assign: fail or pass 
    ## asign: error
    
    error = 'no error'    
    check='Fail'

    if freq_check==True and opt_check==True:
        check='Pass'
        
    elif freq_check==False and opt_check==True:
        error = 'frequency error'
        
    elif opt_check==False:
        error = 'optimisation error'
    
            
    return check, error



def thermal_correction_energy(filename):
    
    ## extract Thermal correction to Gibbs Free Energy and Thermal correction to Enthalpy
    output_line=[]
    with open(filename) as output:
        for line in output:
            output_line.append(line)
            
            
    gibbs_free_energy = 0
    enthalpies = 0

    for l in output_line:
        if'Thermal correction to Gibbs Free Energy='in l:
            gibbs_free_energy = float(l.split()[-1])
    
        if 'Thermal correction to Enthalpy=' in l:
            enthalpies = float(l.split()[-1])
        
    return gibbs_free_energy, enthalpies


def getSPE(filename):
    
    SCF_done_no = 0
    
    SCF_done = []

    with open(filename) as output:
            for line in output:               
                if 'SCF Done:' in line:
            
                    SCF_done += [line]
                    SCF_done_sp = SCF_done[0].split()
                    SCF_done_no = float(SCF_done_sp[4])
    #print(filename)
    #print(SCF_done)
    #print(SCF_done_no)
    #print()
    
    return SCF_done_no

def get_delG(path_opt, get_csv='yes'):
    
    ## for ground state conformers
    ## the spe folder need to be inside the opt file
    
    path_spe=path_opt+'/spe'
    
    opt_files_raw = [f for f in os.listdir(path_opt) if f.endswith('.out')]
    opt_files = natsort.natsorted(opt_files_raw)

    spe_files_raw = [f for f in os.listdir(path_spe) if f.endswith('.out')]
    spe_files = natsort.natsorted(spe_files_raw)


    ## opt df 

    ther_cor_G=[]
    ther_cor_H=[]
    opt_name=[]
    #for p in pass_ls:
    for p in opt_files:
        path=path_opt+'/'+p
        result=thermal_correction_energy(path)
        ther_cor_G.append(result[0])
        ther_cor_H.append(result[1])
        opt_name.append(p[:-4])
    
    
    energy={'idx':list(range(0,len(opt_name))),'name':opt_name, 'ther_cor_H':ther_cor_H,'ther_cor_G':ther_cor_G}
    df_energy=pd.DataFrame(energy)

    ## spe df 
    scf=[]
    spe_name=[]
    for s in spe_files:
        
        path=path_spe+'/'+s
        scf.append(getSPE(path))
        spe_name.append(s[:-8])

    scf_dict={'name':spe_name, 'scf':scf}
    scf_energy=pd.DataFrame(scf_dict)

    ## merge spe_df and opt_df 
    data_df= pd.merge(df_energy, scf_energy, on=['name'])

    ## calculate energy
    data_df['H_hatree']=data_df['ther_cor_H'] + data_df['scf']
    data_df['G_hatree']=data_df['ther_cor_G'] + data_df['scf']
    #data_df['G_hatree']=data_df.apply(lambda x: x['ther_cor_G'] + x['scf'], axis=1)

    min_G=min(data_df['G_hatree'])
    min_H=min(data_df['H_hatree'])

    data_df['delH_kcal']=(data_df['H_hatree'] - min_H)*627.509
    data_df['delG_kcal']=(data_df['G_hatree'] - min_G)*627.509
    
    data_df['delG_kJ']=(data_df['G_hatree'] - min_G)*2625.5

    if get_csv == 'yes':
    
        csv_name=path_opt.split('/')[-1]+'_delG.csv'
        data_df.to_csv(csv_name,index=False)
    
    return data_df




########

def xyz_str(filename):
    
    ##based on https://sites.google.com/site/rangsiman1993/comp-chem/techniques/gaussian-extract-xyz
    ## get xyz from g16 output file as a block of str with /n between lines
    
    output_line=[]
    
    with open(filename) as output:
        for line in output:
            output_line.append(line)
        
    start = 0
    end = 0


    for i in range (len(output_line)):
        if "Standard orientation:" in output_line[i]:
            start = i

    for m in range (start + 5, len(output_line)):
        if "---" in output_line[m]:
            end = m

            break

    ## Convert to Cartesian coordinates format
    ## convert atomic number to atomic symbol

    #xyz=str(len(output_line[start+5 : end]))+'\n'
    #xyz+=filename[:-4]+'\n'
    xyz=''

    for line in output_line[start+5 : end] :
        words = line.split()
        word1 = int(words[1])
        #word3 = str(words[3])
        
        
        if   word1 ==   1 : word1 = "H"
        elif word1 ==   2 : word1 = "He"
        elif word1 ==   3 : word1 = "Li"
        elif word1 ==   4 : word1 = "Be"
        elif word1 ==   5 : word1 = "B"
        elif word1 ==   6 : word1 = "C"
        elif word1 ==   7 : word1 = "N"
        elif word1 ==   8 : word1 = "O"
        elif word1 ==   9 : word1 = "F"
        elif word1 ==  10 : word1 = "Ne"
        elif word1 ==  11 : word1 = "Na"
        elif word1 ==  12 : word1 = "Mg"
        elif word1 ==  13 : word1 = "Al"
        elif word1 ==  14 : word1 = "Si"
        elif word1 ==  15 : word1 = "P"
        elif word1 ==  16 : word1 = "S"
        elif word1 ==  17 : word1 = "Cl"
        elif word1 ==  18 : word1 = "Ar"
        elif word1 ==  19 : word1 = "K"
        elif word1 ==  20 : word1 = "Ca"
        elif word1 ==  21 : word1 = "Sc"
        elif word1 ==  22 : word1 = "Ti"
        elif word1 ==  23 : word1 = "V"
        elif word1 ==  24 : word1 = "Cr"
        elif word1 ==  25 : word1 = "Mn"
        elif word1 ==  26 : word1 = "Fe"
        elif word1 ==  27 : word1 = "Co"
        elif word1 ==  28 : word1 = "Ni"
        elif word1 ==  29 : word1 = "Cu"
        elif word1 ==  30 : word1 = "Zn"
        elif word1 ==  31 : word1 = "Ga"
        elif word1 ==  32 : word1 = "Ge"
        elif word1 ==  33 : word1 = "As"
        elif word1 ==  34 : word1 = "Se"
        elif word1 ==  35 : word1 = "Br"
        elif word1 ==  36 : word1 = "Kr"
        elif word1 ==  37 : word1 = "Rb"
        elif word1 ==  38 : word1 = "Sr"
        elif word1 ==  39 : word1 = "Y"
        elif word1 ==  40 : word1 = "Zr"
        elif word1 ==  41 : word1 = "Nb"
        elif word1 ==  42 : word1 = "Mo"
        elif word1 ==  43 : word1 = "Tc"
        elif word1 ==  44 : word1 = "Ru"
        elif word1 ==  45 : word1 = "Rh"
        elif word1 ==  46 : word1 = "Pd"
        elif word1 ==  47 : word1 = "Ag"
        elif word1 ==  48 : word1 = "Cd"
        elif word1 ==  49 : word1 = "In"
        elif word1 ==  50 : word1 = "Sn"
        elif word1 ==  51 : word1 = "Sb"
        elif word1 ==  52 : word1 = "Te"
        elif word1 ==  53 : word1 = "I"
        elif word1 ==  54 : word1 = "Xe"
        elif word1 ==  55 : word1 = "Cs"
        elif word1 ==  56 : word1 = "Ba"
        elif word1 ==  57 : word1 = "La"
        elif word1 ==  58 : word1 = "Ce"
        elif word1 ==  59 : word1 = "Pr"
        elif word1 ==  60 : word1 = "Nd"
        elif word1 ==  61 : word1 = "Pm"
        elif word1 ==  62 : word1 = "Sm"
        elif word1 ==  63 : word1 = "Eu"
        elif word1 ==  64 : word1 = "Gd"
        elif word1 ==  65 : word1 = "Tb"
        elif word1 ==  66 : word1 = "Dy"
        elif word1 ==  67 : word1 = "Ho"
        elif word1 ==  68 : word1 = "Er"
        elif word1 ==  69 : word1 = "Tm"
        elif word1 ==  70 : word1 = "Yb"
        elif word1 ==  71 : word1 = "Lu"
        elif word1 ==  72 : word1 = "Hf"
        elif word1 ==  73 : word1 = "Ta"
        elif word1 ==  74 : word1 = "W"
        elif word1 ==  75 : word1 = "Re"
        elif word1 ==  76 : word1 = "Os"
        elif word1 ==  77 : word1 = "Ir"
        elif word1 ==  78 : word1 = "Pt"
        elif word1 ==  79 : word1 = "Au"
        elif word1 ==  80 : word1 = "Hg"
        elif word1 ==  81 : word1 = "Tl"
        elif word1 ==  82 : word1 = "Pb"
        elif word1 ==  83 : word1 = "Bi"
        elif word1 ==  84 : word1 = "Po"
        elif word1 ==  85 : word1 = "At"
        elif word1 ==  86 : word1 = "Rn"
        elif word1 ==  87 : word1 = "Fe"
        elif word1 ==  88 : word1 = "Ra"
        elif word1 ==  89 : word1 = "Ac"
        elif word1 ==  90 : word1 = "Th"

    ## copy from atom list.

        xyz+= word1+line[30:-1] + '\n'
    
    return xyz


def format_xyz(line):
    
    # format to sdf style (xyz + symbol part)
    
    line_sp=line.split()

    format_xyz_ls=[]
    for i in line_sp[1:]:
        if i[0] == '-':
            if float(i) > -10:
                format_xyz_ls.append(format(float(i), '.4f'))
            else: 
                format_xyz_ls.append(format(float(i), '.3f'))
        else:
            if float(i) < 10:
                format_xyz_ls.append(' '+format(float(i), '.4f'))
            else: 
                format_xyz_ls.append(' '+format(float(i), '.3f'))
        
    format_xyz_sym=' '*3+format_xyz_ls[0]+' '*3+format_xyz_ls[1]+' '*3+format_xyz_ls[2]+' '+line_sp[0]
    
    if len(format_xyz_sym) == 32:
        format_xyz_sym += ' '
    
    return format_xyz_sym


def get_xyz_lines_format(g16_path):
    
    xyz_lines=xyz_str(g16_path)
    xyz_lines_sp=xyz_lines.split('\n')
    xyz_lines_format=[format_xyz(ln) for ln in xyz_lines_sp[:-1]]
    
    return xyz_lines_format


def get_g16sdf(path, radical = False, rmAtom_ls=[]):
    
    ## corresponding folder of the path need to contain g16 out files and the sdf of the FF conformers
    ## only the first mol in the sdf will be used for the construction of the txt
    
    ## atom to be removed to get the radical structure -- rmAtom_ls
    ## atom index starts from 1
    
    ## get sdf output name 
    get_sdf_name = [f for f in os.listdir(path) if f.endswith('.sdf')]
    
    sdf_file=path+'/'+get_sdf_name[0]
    sdf_file1=path+'/'+get_sdf_name[0]
    
    ##### to accommodate radicals 
    
    if radical == True:
        suppl = Chem.SDMolSupplier(sdf_file,removeHs=False)
        mol = [x for x in suppl][0]
        em1 = Chem.EditableMol(mol)
        for at in rmAtom_ls:
            em1.RemoveAtom(at-1)
        
        m2 = em1.GetMol()
        
        file = open(sdf_file[:-4]+'_r.sdf', "w") 
        file.write(Chem.MolToMolBlock(m2)) 
        file.close()
        
        sdf_file=sdf_file1[:-4]+'_r.sdf'
        
    else:
        sdf_file=path+'/'+get_sdf_name[0]
        
    #####
    
    ## get the and the key index from the sdf file 
    output_line=[]

    
    with open(sdf_file) as output:
        for line in output:
            output_line.append(line)

    end_idx=output_line.index('M  END\n')
    #atom_no=int(output_line[3].split()[0])
    atom_no=int(output_line[3][:3])

    ## import g16 output info 
    group_name=path.split('/')[-1]
    opt_files_raw = [f for f in os.listdir(path) if f.endswith('.out')]
    opt_files = natsort.natsorted(opt_files_raw)
    
    pass_ls = []

    for f in opt_files:
        chk,err = check(path+'/'+f)
        if chk == 'Pass':
            pass_ls.append(f)
        elif chk == 'Fail':
            print(f+': '+err)

    ## get xyz in the right format 
    get_xyz_lines_format_ls=[]
    for opt in  pass_ls:
        get_xyz_lines_format_ls.append(get_xyz_lines_format(path+'/'+opt))
    

    ## aseemble the info
    start_ls1 = output_line[:2]
    start_ls1_ls = [start_ls1+[' '+pass_ls[idx][:-4]+'\n',output_line[3] ] for idx in range(0,len(pass_ls))]

    atom_info = output_line[4:4+atom_no]
    atom_info_noxyz=[l[33:] for l in atom_info]

    end_ls=output_line[4+atom_no:end_idx+1]
    end_ls+=[ '\n','$$$$\n']

    sdf_ls=[]
    for idx in range(0, len(pass_ls)):
        xyz_lines_format=get_xyz_lines_format_ls[idx]
        atom_info_xyz = [xyz_lines_format[idx]+atom_info_noxyz[idx] for idx in range(0,len(xyz_lines_format))]
        one_sdf=start_ls1_ls[idx]+atom_info_xyz+end_ls
        sdf_ls.append(one_sdf)

    
    sdf_txt_ls=[]
    for mol in sdf_ls:
        mol_txt=''
        for i in mol:
            mol_txt+=i
        sdf_txt_ls.append(mol_txt)
    
    sdf_txt=''
    for txt in sdf_txt_ls:
        sdf_txt+=txt
    
    g16_sdf=open(group_name+'_g16.sdf', "w")
    g16_sdf.write(sdf_txt)
    g16_sdf.close()
    
    ##### to accommodate radicals  -- remove _r.sdf file
    if radical == True:
        os.remove(path+'/'+get_sdf_name[0][:-4]+'_r.sdf')
        
    
    