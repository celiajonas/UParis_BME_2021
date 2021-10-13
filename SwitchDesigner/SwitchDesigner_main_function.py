from Bio import SeqIO
import os
import numpy as np
import pandas as pd
from SwitchDesigner_helper_functions import *
from nupack import *
    
def SwitchDesigner(parameters):
    ''' Main function to create the toehold switch. Necessite functions from SwitchDesigner_helper_functions.py.
    Create a folder named after the miRNA tested, a .txt file containing the input variables, a .csv file for toehold candidates, a .csv file containing successful toehold and .txt files containing sequence (RNA or DNA) for each toehold
    
    parameters (dict) ; must contain the following keys :
            "input_seq", "len_unpaired", "len_paired", "output_folder", "min_unpaired", "reporter", "mol_type", "suitable_aa"
    
    
    input_seq (str): path to the miRNA sequence (fasta file)
    
    len_unpaired (int): Length of the unpaired region of the recognition sequence 
    
    len_paired (int): Length of the paired region of the recognition sequence
    
    output_folder (str): path to the output folder
    
    min_unpaired (int): Minimum number of unpaired bases in the secondary structure of the target mRNA for a candidate trigger to be considered
    
    reporter (str): DNA or RNA sequence - Reporter gene or tag to be added to the end of the toehold
    
    mol_type (str): 'DNA' or 'RNA'
    
    suitable_aa (list): Capital letters of amino acids that are wanted for the beginning of the protein sequence
    
    return : dataframe of selected toehold switch with all informations collected.
    
    '''
    ### 1. Prepare parameters
    # Retrieve parameters
    input_seq=parameters["input_seq"]
    length_unpaired=parameters["len_unpaired"]
    length_paired=parameters["len_paired"]
    output_folder=parameters["output_folder"]
    min_unpaired=parameters["min_unpaired"]
    reporter = parameters["reporter"]
    mol_type = parameters["mol_type"]
    suitable_aa=parameters["suitable_aa"]
    
    # Retrieve sequence name from the input_variable.py file
    name_sequence=input_seq.split('/')[-1][0:-6] 
    
    # Create a folder "output" containing each sequence result. Create a second folder for each sequence
    seq_folder=output_folder+"/"+name_sequence
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(seq_folder):
        os.makedirs(seq_folder)
        
    # Make a copy of the input parameter file in the new directory
    f= open(seq_folder+"/input_variable.txt","w+")
    variable_str="Sequence : "+str(input_seq)
    variable_str=variable_str+"\nLength unpaired : "+str(length_unpaired)
    variable_str=variable_str+"\nLength paired : "+str(length_paired)
    variable_str=variable_str+"\nOutput_folder : "+str(output_folder)
    variable_str=variable_str+"\nReporter : "+str(reporter)
    variable_str=variable_str+"\nmol_type : "+str(mol_type)
    variable_str=variable_str+"\nMinimum unpaired : "+str(min_unpaired)
    f.write(variable_str)
    f.close() 

    ### 2. Prepare the mRNA sequence
    # Parse the input sequence
    input_seq_record = next(SeqIO.parse(input_seq, 'fasta'))

    # If input sequence is DNA, parse and transcribe
    if mol_type == 'DNA':   
        full_input_seq = DNAtoRNA(str(input_seq_record.seq))
        
    # If it is RNA, make sure it is in upper letters.
    elif mol_type == 'RNA':
        full_input_seq = str(input_seq_record.seq).upper()

    # Calculate the secondary structure with NUPACK
    my_model = Model(material='rna95',celsius=37)
    mfe_structures = mfe(strands=full_input_seq, model=my_model)

    # Store values
    mRNA_dict={}
    mRNA_dict['sequence']=full_input_seq 
    mRNA_dict['structure']=mfe_structures[0].structure 
    mRNA_dict['paired']=mfe_structures[0].structure.pairlist()
    mRNA_dict['energy']=mfe_structures[0].energy

    selected = []
    results = []
    
    ### 3. Find candidates : Loop through the miRNA sequence to find a suitable hairpin start position
    for hairpin_start_pos in range(length_unpaired, len(mRNA_dict['sequence']) - length_paired):

        # Store the sequence of possible hairpin's base
        hairpin_base = mRNA_dict['sequence'][hairpin_start_pos:hairpin_start_pos + 3]
        
        # Discard cases for which the base of the hairpin does not have two weak pairs and a strong one
        if hairpin_base.count('A') + hairpin_base.count('U') != 2:
            continue
        
        # store values
        end_pos = hairpin_start_pos + length_paired
        start_pos = hairpin_start_pos - length_unpaired
        sub_seq = mRNA_dict['sequence'][start_pos:end_pos]
        mRNA_structure=str(mRNA_dict['structure'])
        sub_struc = mRNA_structure[start_pos:end_pos-1]
        
        # Count the number of unpaired bases in the miRNA subsequence 
        unpaired = sub_struc.count('.')
        
        
        # Discard cases for which the number of unpaired bases in the subsequence doesn't match the minimum number of unpaired bases.
        if unpaired < min_unpaired:
            continue 
        
        # Store values
        hairpin_data=[unpaired, sub_struc, sub_seq, start_pos + 1, end_pos, length_unpaired, length_paired]
        selected.append(hairpin_data)
        index=len(selected)-1
        
        ### 4. Generate toehold 
        toehold_dict={}
        toehold_dict=generate_toehold(parameters, sub_seq, seq_folder, index)
        
        # Store values
        toehold_dict["start_pos"]=start_pos + 1
        toehold_dict["end_pos"]=end_pos
            
        ### 5. Check how well the toeholds bind to the target in the mRNA
        binding_data=toehold_binding(mRNA_dict, toehold_dict)
            
        ### 6. Verify the protein produced by the toehold switch    
        protein=translate(toehold_dict["sequence_RNA"])
        binding_data.append(protein)
         
        # Discard cases for which the protein generated by the reporter 1. do not begin by M and 2. the 2nd and 3rd amino acids are not one of the suitable amino acids.  
        if protein[0]!="M":
            continue
        if protein[1] not in suitable_aa:
            continue
        if protein[2] not in suitable_aa:
            continue
        
        # Add this row to the list of results
        results_list=[]
        for i in hairpin_data: 
            results_list.append(i)
        for i in binding_data:
            results_list.append(i)
        results.append(results_list)

    
    ### 7. Create csv files for results
    tmp = pd.DataFrame(np.array(selected), columns=['Non paired count', 'Structure', 'Sequence', 'Start', 'End', 'Length unpaired trigger', 'Length paired trigger'])
    df = tmp.sort_values(by = 'Non paired count', ascending = False)
    df.to_csv(path_or_buf=os.path.join(seq_folder, 'toehold_candidates.csv'), sep = '\t', index = False)
    
    results_df = pd.DataFrame(results, columns=['Non paired count', 'Structure', 'Sequence', 'Start', 'End', 'Length unpaired trigger', 'Length paired trigger', 'Index', 'Binding_energy_toehold_mRNA', 'Binding_energy_toehold', 'Binding_energy_mRNA', 'GC content','Protein_sequence'])
    
    sorted_results = results_df.sort_values(['Binding_energy_toehold_mRNA'], ascending = False)
    sorted_results.to_csv(path_or_buf=os.path.join(seq_folder, 'selected_toeholds_results.csv'), sep = '\t', index = False)
    
    return sorted_results
