from SwitchMiDesigner_main_function import *
#import easygui
import os

directory = r'/Users/igemuparisbme/miRNAs' # path to the folder containing miRNA
folder = '/Users/igemuparisbme/Output' # path to the output folder

# can also be done with easygui
#directory = easygui.diropenbox() 
#folder = easygui.diropenbox() 

parameters ={"input_seq":' ', 
             "output_folder":folder,
             "min_unpaired":4,
             "reporter":'ATGCGTAAA',
             "mol_type":'RNA',
             "len_unpaired":10,
             "len_paired":4,
             "suitable_aa":["L","F","Y","I","N","K"]}

for filename in os.listdir(directory):
    if filename.endswith(".fasta"):
        input_seq=os.path.join(directory, filename)
        try:
            parameters["input_seq"]=input_seq
            SwitchMiDesigner(parameters)
        except:
            continue
    else:
        continue 