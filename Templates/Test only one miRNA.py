#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' 1. Import libraries '''
from SwitchMiDesigner_main_function import *
import easygui

''' 2. Define parameters '''
file = '/Users/igemuparisbme/miRNA.fasta'
folder = '/Users/igemuparisbme/Output'

# can also be done with easygui
#file = easygui.fileopenbox()
#folder = easygui.diropenbox()

parameters ={"input_seq":file, 
             "output_folder":folder,
             "min_unpaired":4,
             "reporter":'ATGCGTAAA',
             "mol_type":'RNA',
             "len_unpaired":10,
             "len_paired":9,
             "suitable_aa":["L","F","Y","I","N","K"]}

''' 3. Run the function '''

SwitchMiDesigner(parameters)

