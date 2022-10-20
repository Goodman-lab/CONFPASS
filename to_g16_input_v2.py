#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 10:34:37 2022

@author: chingchinglam

update from version 1:
    update functions to accommodate radicals 

"""



def right_spacing(coor_list):
    coor_list_space=[]
    for cline in coor_list:
        split_coor=cline.split()
        new_coor=''
        space1=18
        space2=4
        space3=4
        if len(split_coor[0]) == 2:
            space1=17
        for o in split_coor:  
            if '-' in o:
                index_list=[]
                index_list.append(split_coor.index(o))
                if 1 in index_list:
                    space1=17
                elif 2 in index_list:
                    space2=3
                elif 3 in index_list:
                    space3=3
        new_coor=str(split_coor[0])+space1*' '+str(split_coor[1])+'00'+space2*' '+str(split_coor[2])+'00'+space3*' '+str(split_coor[3])+'00\n'
        coor_list_space.append(new_coor)
    return coor_list_space



def sdf2gjfs(filename,keywords,space, priority_ls, per, radical = False, rmAtom_ls=[]):

    sdf_output=[]
    
    with open(filename) as output:
        for line in output:
            sdf_output+= [line]
        
    END_index_list =[k for k,d in enumerate(sdf_output) if d=='M  END\n']

    atom_no=int(sdf_output[3][0:3])

    len_mol=END_index_list[0]

    division1=[]
    for j in END_index_list:
        indv_mol=sdf_output[j-len_mol:j+1]
        #atom_no=int(indv_mol[3].split(' ')[1])
    
        indv_mol_atom1=indv_mol[4:4+atom_no]
        indv_mol_atom2=[j[30:35]+j[0:31] for j in indv_mol_atom1]

        division1.append(indv_mol_atom2)

    ### to accommodate radicals 
    if radical == True:
        
        division1_rm=[]
        for conf in division1:
            division1_rm.append([conf[idx] for idx in range(0,len(conf)) if idx+1 not in rmAtom_ls])
            
        division2=[right_spacing(ls) for ls in division1_rm]
        
    else:

        division2=[right_spacing(ls) for ls in division1]
    
    
    ###
    
    division3=[]
    
    for structure in division2: 
        coor_content=''
        for individual_structure in structure:
            coor_content+=str(individual_structure)
        content=keywords+coor_content+space
        division3.append(content)    
    
    
    len_priority_ls = len(priority_ls)   
    priority_ls1=[x+1 for x in priority_ls]
            
    len_per_priority_ls = round(len(priority_ls)*per)
    
    select = [y+1 for y in priority_ls]
    select1 = select[:len_per_priority_ls]
    select2 = priority_ls[:len_per_priority_ls]
    
    notselect = select[len_per_priority_ls:]
    #notselect2 = priority_ls[len_per_priority_ls:]

    division4=[]
    for o in select2:
        division4.append(division3[o])
    
    
    filename_list2=[]
    for num2 in select1:
        filename_list2.append(filename[:-4]+'_'+str(num2)+'.gjf')
               
    for c,d in zip(filename_list2,division4):
            #print(c,'\n',d)
        file = open(c, "w") 
        file.write(d) 
        file.close() 
    
    
    
    result_txt='CONFPASS - summary \n \n priority_ls: '+str(priority_ls1)+'\n' 
    result_txt+='\n To be optimised (with .gjf files generated):'+str(select1)+'\n \n waitlisted: '
    result_txt+=str(notselect)+ '\n \n Total number of conformers: '+str(len_priority_ls) 
    
    if radical == True:
        result_txt+='\n radical: '+ str(radical) 
        result_txt+=' \n  rmAtom_ls: '+str(rmAtom_ls)
    
    
    file = open(filename[:-4]+'_summary.txt', "w") 
    file.write(result_txt) 
    file.close() 
    
    #print(result_txt)
    


def sdf2gjfs_v2(filename,keywords,space,numbering, by_num='Yes', radical = False, rmAtom_ls=[]):

    sdf_output=[]
    
    with open(filename) as output:
        for line in output:
            sdf_output+= [line]
        
    END_index_list =[k for k,d in enumerate(sdf_output) if d=='M  END\n']

    atom_no=int(sdf_output[3][0:3])

    len_mol=END_index_list[0]

    division1=[]
    for j in END_index_list:
        indv_mol=sdf_output[j-len_mol:j+1]
        #atom_no=int(indv_mol[3].split(' ')[1])
    
        indv_mol_atom1=indv_mol[4:4+atom_no]
        indv_mol_atom2=[j[30:35]+j[0:31] for j in indv_mol_atom1]

        division1.append(indv_mol_atom2)
        
    
    ### to accommodate radicals 
    if radical == True:
        
        division1_rm=[]
        for conf in division1:
            division1_rm.append([conf[idx] for idx in range(0,len(conf)) if idx+1 not in rmAtom_ls])
            
        division2=[right_spacing(ls) for ls in division1_rm]
        
    else:

        division2=[right_spacing(ls) for ls in division1]
    
    ###
    
    division3=[]
    
    for structure in division2: 
        coor_content=''
        for individual_structure in structure:
            coor_content+=str(individual_structure)
        content=keywords+coor_content+space
        division3.append(content)    
    
    
    name = filename.split('/')[-1][:-4]
    
    if by_num=='No':
        
        filename_ls=[]
        num_ls=list(range(1,len(division3)+1))
        for idx in num_ls:
            filename_ls.append(name+'_'+str(idx)+'.gjf')
        
            #double check
        #for b,a in zip(filename_ls,division3):
            #print(b,'\n',a)

        #write the file
        for ifilename, gcontent in zip(filename_ls,division3):
            file = open(ifilename, "w") 
            file.write(gcontent) 
            file.close()

    else:
                    
        select=[y-1 for y in numbering]
        select1=[y+1 for y in select]

        division4=[]
        for o in select:
            division4.append(division3[o])
    
    
        filename_list2=[]
        for num2 in select1:
            filename_list2.append(name+'_'+str(num2)+'.gjf')
               
        for c,d in zip(filename_list2,division4):
            ## write the file
            #print(c,'\n',d)
            file = open(c, "w") 
            file.write(d) 
            file.close() 
            
            
            
            
            
            