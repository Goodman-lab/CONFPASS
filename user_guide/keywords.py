#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 19:43:57 2022

@author: chingchinglam
"""

def get_keywords():
    
    keyword = '''%nprocshared=32
%mem=4GB
# opt freq b3lyp/6-31g(d) int=ultrafine empiricaldispersion=gd3 

Title Card Required
    
0 1
'''
    
    space='''
'''
    
    conf_idx_ls =[[ 29, 50, 24, 26, 14, 36, 13, 43, 49]]
    
    return keyword, space, conf_idx_ls