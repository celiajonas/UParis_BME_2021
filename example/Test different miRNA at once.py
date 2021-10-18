from SwitchMiDesigner_main_function import *
import os

directory = r'/Users/igemuparisbme/example/miRNA' # path to the folder containing miRNA

folder_12up_8p = r'/Users/igemuparisbme/example/Output/Results - 12 unpaired 8 paired'
folder_10up_9p = r'/Users/igemuparisbme/example/Output/Results - 10 unpaired 9 paired'
folder_10up_4p = r'/Users/igemuparisbme/example/Output/Results - 10 unpaired 4 paired'


parameters ={"input_seq":' ', 
             "output_folder":' ',
             "min_unpaired":4,
             "reporter":'ATGCGTAAA',
             "mol_type":'RNA',
             "len_unpaired":12,
             "len_paired":8,
             "suitable_aa":["L","F","Y","I","N","K"]}

for filename in os.listdir(directory):
    if filename.endswith(".fasta"):
        input_seq=os.path.join(directory, filename)
        
        # Unpaired 12 Paired 8
        try:
            parameters["input_seq"]=input_seq
            parameters["output_folder"]=folder_12up_8p
            parameters["len_unpaired"]=12
            parameters["len_paired"]=8
            SwitchMiDesigner(parameters)
        except:
            continue
        
        # Unpaired 10 Paired 9 
        try:
            parameters["input_seq"]=input_seq
            parameters["output_folder"]=folder_10up_9p
            parameters["len_unpaired"]=10
            parameters["len_paired"]=9
            SwitchMiDesigner(parameters)
        except:
            continue
        
        # Unpaired 10 Paired 4
        try:
            parameters["input_seq"]=input_seq
            parameters["output_folder"]=folder_10up_4p
            parameters["len_unpaired"]=10
            parameters["len_paired"]=4
            SwitchMiDesigner(parameters)
        except:
            continue
    else:
        continue 