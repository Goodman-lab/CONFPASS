#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 18:25:04 2022

@author: chingchinglam

14072022 - update to accommodate radicals; update pas test - more nor_nmx options
15082022 - update to change the pas test nor_nx setting and model -- 0.5 --> 0.99
14092022 - include temperature variable for pas test; check the completion of the optimisation calculation;
            pas.add_cal() function 
04022023 - introduce the plugin: pas.preparation(chk_structure=False)
19042023 - introduce repeat parameter pas.make_prediction -- conduct the prediction making process 
            for x time at reopt (n-x/m) to (n/m); repeat = x

            
"""

__version__ = "v19042023"

import numpy as np
from datetime import date
from optparse import OptionParser
import pandas as pd
import pickle

import json
import os 
import sys 
sys.path.insert(0, os.path.abspath(__file__)[:-11])


from clustering_dih_v7 import get_cluster_df_multi
from GetPriority_v3 import GetPriority
from to_g16_input_v2 import sdf2gjfs, sdf2gjfs_v2

from get_dft_output_v2 import get_g16sdf, get_delG
from BinsTest_rmsCheck import rmsCheck
from MolFrac_ML import get_descriptor

import xyz2mol as x2m

today = date.today()
date_str=today.strftime("%d%m%Y")

class conp:
    
    def __init__(self, sdf_files):
        
        self.sdf_files = sdf_files
        

    def get_priority(self, x=0.8, x_as=0.2, n=3, method = 'pipe_x_as'):

        self.methodf = method
        
        ## perform the clustering calculation 
        self.clustering_result = get_cluster_df_multi(self.sdf_files)
        
        ## generate the priority list 
        test1 = GetPriority(self.clustering_result)
        test1.priority_df(x_=x, x_de_=x_as, n_=n, method_ = method)
        
        self.priority_df = test1.priority_ls_df
        
        p_ls= self.priority_df['priority_ls'].tolist()
        
        reidx_p_ls=[]
        for ls in p_ls: 
            reidx_p_ls.append([i+1 for i in ls])
        
        name_ls=self.priority_df['name'].tolist()
        
        self.reidx_priority_df = pd.DataFrame({'name':name_ls,'priority_ls':reidx_p_ls})
        
        
    def priority2csv(self):
        
        ## apply after self.get_priority()
        ## create a .csv version of the priority list data frame 
        
        self.reidx_priority_df.to_csv(self.methodf+'_'+date_str+'.csv',index=False)
        
        print('o '+self.methodf+'_'+date_str+'.csv created')
        
    def priority2gjf(self, keywords, space, per, radical = False, rmAtom_ls=[]):
        
        ## apply after self.get_prority()
        ## generate .gjf input for g16 calculations according to the priority list
        
        if radical == True:
            
            ## can only process one radical molecule at a time when generating gjf
            if len(self.sdf_files) > 1:
                print('Error: Radical setting: Can only process one molecule at a time when generating gjf')
                pass
            
            sdf2gjfs(self.sdf_files[0], keywords, space, self.priority_df['priority_ls'][0], per, radical = True, rmAtom_ls=rmAtom_ls)
        
        else:
        
            for idx in range(0, len(self.sdf_files)):
                sdf2gjfs(self.sdf_files[idx], keywords, space, self.priority_df['priority_ls'][idx], per)
            
        
        
    
    def get_gjf(self, keywords, space, conf_idx_ls,  radical = False, rmAtom_ls=[]):
        
        ## apply after self.get_priority()
        ## create gjf files for DFT calculations based on a list of conformer indexes specified by the user
        
        if radical == True:
            
            ## can only process one radical molecule at a time when generating gjf
            if len(self.sdf_files) > 1:
                print('Error: Radical setting: Can only process one molecule at a time when generating gjf')
                pass
            
            sdf2gjfs_v2(self.sdf_files[0], keywords, space, conf_idx_ls[0], radical = True, rmAtom_ls=rmAtom_ls)
        
        else:
            for idx in range(0, len(self.sdf_files)):
                sdf2gjfs_v2(self.sdf_files[idx], keywords, space, conf_idx_ls[idx])
        

class pas:
    
    def __init__(self, path):
        
        self.path = path
        self.molname=path.split('/')[-1]
    
    def preparation(self, molf='1', ra= False, rm_ls=[], chk_structure=False, 
                    p_x=0.8, p_x_as=0.2, p_n=3, p_method = 'pipe_x_as'):
        
        ## extract energy info from g16 output files 
        self.delG_df=get_delG(self.path , get_csv='yes')
        
        origin_structure_no=len(self.delG_df)
        structure_no_postchk=len(self.delG_df)
        
        ## extract structure information from g16 output files 
        get_g16sdf(self.path, radical = ra, rmAtom_ls=rm_ls)
        self.g16sdf_name=self.molname+'_g16.sdf'
        
        ## check structures function starts here
        ### sometimes bond forming or breaking may occur upon reoptimisations at DFT leve
        ### the below plugin check the structures and separate different chemical systems into different sdf
        ### associated delG csv is also given and updated for each sdf file
        ### CONFPASS will be applied on the sdf file with the same chemical system as the sdf at the FF level
        ### or sdf file with the most conformer (if the previous is not available 
        ###  We used the script from https://github.com/jensengroup/xyz2mol/blob/master/xyz2mol.py 
        ### (which is based on the work of Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777) to convert xyz coordinates to a rdkit.Chem.mol object.
        
        self.chk_structure=chk_structure
        
        if chk_structure==True:
        
            self.delG_csv_name=self.molname+'_delG.csv'
            self.sdf_csv_path=self.path +'/'+ self.molname +'.sdf'
        
            x2m.mol_separations(self.g16sdf_name,self.delG_csv_name,self.sdf_csv_path)
            
            self.delG_df = pd.read_csv(self.delG_csv_name)
            
            structure_no_postchk=len(self.delG_df)
            
            if structure_no_postchk < origin_structure_no:
                print(str(structure_no_postchk)+'/'+str(origin_structure_no)+' pass the structure check')
            
        
        self.structure_chk_df=pd.DataFrame({'name':[self.molname],'structure_no_prechk':[origin_structure_no], 'structure_no_postchk':
                                            [structure_no_postchk]})
            
            
        ## check structures function ends here
        
        ## perform rms calculations 
        get_rmsCheck_df=rmsCheck(self.g16sdf_name, self.delG_df, self.molname)
        get_rmsCheck_df.get_rmsCheck_cluster()
        self.rmsCheck_cluster_df=get_rmsCheck_df.rmsCheck_cluster_df
        
        ## load the RF model 
        if molf == '1':
            filename = os.path.abspath(__file__)[:-11]+'LR_model_24042023_x10.sav'
            self.model = pickle.load(open(filename, 'rb'))
            
        #elif molf == '0.99':
        #    filename = os.path.abspath(__file__)[:-11]+'RF_model_09042023_x99.sav'
        #    self.model = pickle.load(open(filename, 'rb'))
            
        self.molf = molf
        
        
        ## generate the priority list 
        ## sdf file located inside the directory in the path specify 
        self.sdf_files = [self.path +'/'+ self.molname +'.sdf']
        
        ptest = conp(self.sdf_files)
        ptest.get_priority(method = p_method, x=p_x, x_as=p_x_as, n= p_n)
        
        self.priority_df=ptest.priority_df
    
    
    
    def process_prediction(self, label, prob):
        
        ## this function is used in the make_prediction function 
        ## not to be excuted by itself 
        
        completion ='complete'
        
        if label == 1:
            prob_v=prob[1]
        elif label == 0:
            prob_v=prob[0]
            
            completion ='incomplete'
        
        self.label = completion
        self.prob_v = round(prob_v,3)
        ## 0 - reoptimisation incomplete 
        ## 1 - reoptimisation completed
        
        if self.molf =='1':
            A = -6.293126587013565
            
       # elif self.molf =='0.1':
       #     A= -10.83234412
            
        
        per_conf=(1 / (1 + np.exp(A  * (prob_v -0.5))))*100
        
        
        if label == 0:
            per_conf2=100-per_conf
        else:
            per_conf2=per_conf

        self.confidence=round(per_conf2,3)
        
        
    
    def make_prediction(self, T=298.15, print_result=True, repeat=0):
        
        if repeat == 0:
            pls_len_ls=[None]
            
        else:
            pls_len_ls=[-repeat+i for i in range(repeat)][1:]+[None]
        
        
        label_ls=[]
        prob_ls=[]
        
        for pls in pls_len_ls:
        
        
        ## generate the descriptor array 
            descriptor_y, self.reopt_idx_ls =get_descriptor(self.rmsCheck_cluster_df, self.priority_df, self.delG_df,temp=T, pls_len=pls)
        
        ## make prediction 
            label= self.model.predict(descriptor_y)[0]
            label_ls.append(label)
        #print(label)
            prob = list(self.model.predict_proba(descriptor_y)[0])
            prob_ls.append(prob)
        
        
        completion_ls=[]
        confidence_ls=[]
        prob_v_ls=[]
        
        for l,p in zip(label_ls, prob_ls):
            
            self.process_prediction(l, p)
            completion_ls.append(self.label)
            confidence_ls.append(self.confidence)
            prob_v_ls.append(self.prob_v)
        
        ## calculate reopt
        
        total_conformer_no = len(self.priority_df['priority_ls'].tolist()[0])
        reoptimized_conf_no = len(self.delG_df)
        self.current_ropt= round(reoptimized_conf_no/total_conformer_no,3)
        

        if repeat != 0:
            self.reopt_ls=[round((reoptimized_conf_no-i)/total_conformer_no,3) for i in range(repeat)]
            self.confidence_ls=confidence_ls

        
        ### print out the result
        if print_result==True:
        
            print('reoptimisation: '+ self.label +'; confidence level: ' + str(self.confidence))
            print('probability: ' + str(self.prob_v))
            print('current ropt (number of optimised / total number of conformers): '+str(self.current_ropt))
            
            if repeat != 0:
                
                print('breakdown: reoptimisation: '+ str(completion_ls)+'; confidence level: '+str(confidence_ls))
                print( 'ropt list: '+str(self.reopt_ls))
        
         
        
    def add_cal(self, keywords, space,  radical = False, rmAtom_ls=[], per = 0.1):
        
        if self.chk_structure==True:
        
            delG_df = get_delG(self.path , get_csv='no')
        
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
            self.reopt_idx_ls =reopt_idx_ls1
        

        complete_priority_ls = self.priority_df['priority_ls'].tolist()[0]
        left_ls = [i for i in complete_priority_ls if i not in self.reopt_idx_ls]
        
        
        ## per = percentage of left_ls to be reoptimised in this round 
        tobe_reopt_ls = left_ls[:round(len(complete_priority_ls)*per)]
        tobe_reopt_ls2 = [i+1 for i in tobe_reopt_ls]
        
        ## create gjf files for DFT calculations based on tobe_reopt_ls
        
        if radical == True:
            sdf2gjfs_v2(self.sdf_files[0], keywords, space, tobe_reopt_ls2, radical = True, rmAtom_ls=rmAtom_ls)
        
        else:

            sdf2gjfs_v2(self.sdf_files[0], keywords, space, tobe_reopt_ls2)
            
        
        ## update // generate a new summary file 
        
        complete_priority_ls2 = [i+1 for i in complete_priority_ls]
        reopt_idx_ls2= [i+1 for i in self.reopt_idx_ls]
        notselect = [i for i in complete_priority_ls2 if i not in reopt_idx_ls2 and i not in tobe_reopt_ls2]
        
        result_txt='CONFPASS - summary \n \n priority_ls: '+str(complete_priority_ls2)+'\n' 
        result_txt+='\n optimised:'+str(reopt_idx_ls2)+'\n' 
        result_txt+='\n To be optimised (with .gjf files generated):'+str(tobe_reopt_ls2)+'\n \n waitlisted: '
        result_txt+=str(notselect)+ '\n \n Total number of conformers: '+str(len(complete_priority_ls2)) 
        
        if radical == True:
            result_txt+='\n radical: '+ str(radical) 
            result_txt+=' \n  rmAtom_ls: '+str(rmAtom_ls)
        
        
        file = open(self.molname+'_summary.txt', "w") 
        file.write(result_txt) 
        file.close()
        
        




#########################################
##
##             Execution 
##    
#########################################

    


def main():

    parser = OptionParser()

    parser.add_option('-m',dest='m', help='priority list assembling method', default='pipe_x_as',
                    choices=('pipe_x_as', 'pipe_x','pipe_as','nth','random','ascend'))

    parser.add_option('-x', help='hyperparameter for the pipe_x method, default = 0.8 i.e. select clustering result when n_clusters=total_conformer_no*x', 
                  dest='x', default=0.8)

    parser.add_option('--x_as', help='hyperparameter for the pipe_x_as method, default = 0.2 i.e. 20% (pipe x priority list) + 80%(pipe d priority list)', 
                  dest='x_as', default=0.2)

    parser.add_option('-n', help='hyperparameter for the nth method, default = 3 i.e. prioritise every 3rd conformer', 
                  dest='n', default=3)

    parser.add_option('--path', help='path for pas test', dest='path')

    parser.add_option('--pas', dest='pas', action="store_true", help='predict the completeness of the reoptimisation process', 
                  default=False)

    parser.add_option('--pas_multi', dest='pas_multi', action="store_true", help='predict the completeness of the reoptimisation process - for multiple molecules', 
                  default=False)

    parser.add_option('--csv', dest='csv', action="store_true", help='get csv of the priority list dataframe', 
                  default=False)

    parser.add_option('--p2gjf', dest='p2gjf', action="store_true", help='execute priority2gjf()', 
                  default=False)

    parser.add_option('--togjf', dest='togjf', action="store_true", help='get_gjf()', 
                  default=False)

    parser.add_option('--per', help='percentage of the conformers to be converted to gjf format in the priority list', 
                  dest='per', default=0.2)
    
    
    #parser.add_option('--mx',dest='mx', help='nor_nx=mx; mol fraction', default='1',
    #                choices=('1'))
    
    parser.add_option('-T', help='Temperature setting (required for pas test)', 
                  dest='T', default=298.15)
    
    
    ## for the radicals 
    
    parser.add_option('--radical', dest='radical', action="store_true", help='process radicals (with a different number of atoms compared to the pseudo structure)', 
                      default=False)
    
    parser.add_option('--rmatom', dest='rmatom', help='process radicals (with a different number of atoms compared to the pseudo structure) - list of atom to be removed',
                      default='[]')
    
    
    

    (options, args) = parser.parse_args()



    files = []
    if len(sys.argv) > 1:
        for elem in sys.argv[1:]:
            try:
                if '.sdf' in elem:
                    files.append(elem)
       

            except IndexError: pass
     

    if options.pas == False and options.pas_multi == False:
    
    ## part 1: generate priority list     
    
        test = conp(files)
        test.get_priority(method = options.m, x=float(options.x), x_as=float(options.x_as), n= int(options.n))
        print('method: '+options.m)
        for idx in range(0,len(test.priority_df)):
            print(test.reidx_priority_df['name'][idx])
            print(test.reidx_priority_df['priority_ls'][idx])
    
    
        if options.csv == True:
            test.priority2csv()
        
        
        if options.p2gjf == True:
            
            rm_ls=json.loads(options.rmatom)
            
            from keywords import get_keywords
            
            keyword, space, conf_idx_ls = get_keywords()
            test.priority2gjf(keyword, space, float(options.per), radical = bool(options.radical), rmAtom_ls=rm_ls)
    
        if options.togjf == True:
            
            rm_ls=json.loads(options.rmatom)
        
            from keywords import get_keywords
        
            keyword, space, conf_idx_ls = get_keywords()
            test.get_gjf(keyword, space, conf_idx_ls, radical = bool(options.radical), rmAtom_ls=rm_ls)
        
    

    elif options.pas == True and options.pas_multi == False:
        
        ## part 2: pas test -- more than one molecules 
        
        if bool(options.radical) == True:
            rmat_ls=json.loads(options.rmatom)
            
            test1 = pas(options.path)
            test1.preparation(p_x=float(options.x), p_x_as=float(options.x_as), p_n=int(options.n), p_method = options.m,ra=True, rm_ls=rmat_ls)
            test1.make_prediction( T= options.T)
        
        else:
        
            test1 = pas(options.path)
            test1.preparation(p_x=float(options.x), p_x_as=float(options.x_as), p_n=int(options.n), p_method = options.m)
            test1.make_prediction( T= options.T)
    
    
    elif  options.pas == False and options.pas_multi == True:
        
        #part 2: pas test -- one molecule only 
        
        if bool(options.radical) == True:
            print('Error: the programme is unable to cope with multiple radical molecules at the moment')
            
            pass
        
        sub_folders = [name for name in os.listdir(options.path) if os.path.isdir(os.path.join(options.path, name))]
        sub_folders_path =[options.path+'/'+i for i in sub_folders]
        
        #perform the analyses and compile results together 
        name_ls=[]
        label_ls=[]
        prob_ratio_ls=[]
        confidence_ls=[]
        
        for p in sub_folders_path:
            print()
            print()
            
            try:
                print(p)
                test1 = pas(p)
                test1.preparation(p_x=float(options.x), p_x_as=float(options.x_as), p_n=int(options.n), p_method = options.m)
                test1.make_prediction( T= options.T)
            
                name_ls.append(test1.molname)
                label_ls.append(test1.label)
                prob_ratio_ls.append(test1.prob_ratio)
                confidence_ls.append(test1.confidence)
                
            except:
                print('error with: '+ p)
    
        
        result_df= pd.DataFrame({'name':name_ls, 'label': label_ls, 'prob_ratio': prob_ratio_ls, 'confidence': confidence_ls})
        
        if options.csv == True:
            print()
            result_df.to_csv('pas_'+date_str+'.csv',index=False)
            print('o pas_'+date_str+'.csv created')

    
if __name__ == "__main__":
    main()
    
    

