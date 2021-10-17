from nupack import *
from Bio.Seq import Seq
import os
import csv

def DNAtoRNA(dnaseq):
    ''' Translate a DNA sequence to an RNA sequence by replacing T with U. '''
    rnaseq=dnaseq.upper().replace('T', 'U')
    return rnaseq

def RNAtoDNA(rnaseq):
    ''' Translate an RNA sequence to a DNA sequence by replacing U with T. '''
    dnaseq=rnaseq.upper().replace('U', 'T')
    return dnaseq

def translate(rnaseq):
    ''' Translate an RNA sequence to an amino-acid sequence if there is a start codon. 
    
    
    rnaseq (str): an RNA sequence 
    
    
    return :  amino-acid sequence (str) 
    '''
    
    # RNA codon table
    table = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
        "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
        "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
       "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
       "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
        "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
        }
    
    # Loop through the RNA sequence to find an AUG codon and store its position
    protein =""
    startpos="NaN"
    for i in range(0, len(rnaseq), 1):
        codon = rnaseq[i:i + 3]
        if codon=="AUG":
            startpos=i
            
        if startpos != "NaN":
            break
    
    # If a start codon is found, translate the RNA sequence to a protein using the AUG codon position and the RNA codon table 
    if startpos != "NaN":
        for i in range(startpos, len(rnaseq), 3):
            codon = rnaseq[i:i + 3]
            protein+= table[codon]
        
    return protein    
    
def generate_toehold(parameters, trigger, seq_folder, index):
    ''' Generate a toehold sequence.
    
    
    parameters : dict containing len_unpaired (int) and reporter (str)
    
    trigger : string (target trigger RNA sequence)
    
    seq_folder : string (path to the sequence file)
    
    index : int (toehold id)
    
    
    return : dictionnary containing informations about the toehold (DNA sequence, RNA sequence, structure, pair list, index and energy) if one is generated 
    '''
    
    len_unpaired=parameters["len_unpaired"]
    reporter = parameters["reporter"]
    
    # Get the recognition sequence for our toehold
    recognition_sequence = get_switch_recognition_seq(trigger)
        
    # Design the starting sequence for the switch
    toehold = get_full_switch(recognition_sequence, reporter, len_unpaired) 
    
   
    if toehold != 'Stop':
        toehold=str(toehold)
        toehold=DNAtoRNA(toehold)
        toehold_DNA=RNAtoDNA(toehold)
        
        # Use nupack to get the secondary structure for our switch
        my_model = Model(material='rna95',celsius=37)
        mfe_structures = mfe(strands=toehold, model=my_model)
        
        write_toehold(seq_folder,toehold,'RNA',index)
        write_toehold(seq_folder,toehold_DNA,'DNA',index)
        
        # Read the secondary structure for the toehold we just generated
        Toehold_Dict={}
        Toehold_Dict['sequence_DNA']=toehold_DNA
        Toehold_Dict['sequence_RNA']=toehold
        Toehold_Dict['structure']=mfe_structures[0].structure 
        Toehold_Dict['paired']=mfe_structures[0].structure.pairlist()
        Toehold_Dict['index']=index
        Toehold_Dict['energy']=mfe_structures[0].energy
        
    return Toehold_Dict    

def toehold_binding(mRNA_dict, toehold_dict):
    ''' This function check how well the toehold binds to the target in the mRNA.
    
    
    mRNA_dict : dict; contains the mRNA sequence "sequence" (str) and energy binding "energy" ()
    
    toehold_dict : dict; contains the toehold RNA sequence "sequence_RNA" (str), "index" (), "start_pos" () and "end_pos" ()
        
    
    return : list containing data on mRNA-Toehold binding (index, toehold-mRNA binding energy, toehold binding energy, mRNA binding energy, gc content) 
    '''
    
    toehold_seq=toehold_dict['sequence_RNA']
    mRNA_seq=mRNA_dict['sequence']
    mRNA_binding_energy=mRNA_dict['energy']
    
    binding_data=[]
    binding_data.append(toehold_dict['index'])
    
    # Run NUPACK
    toehold_mRNA = [toehold_seq, mRNA_seq]
    my_model = Model(material='rna95',celsius=37)
    mfe_structures = mfe(strands=toehold_mRNA, model=my_model)
    
    # Read the secondary structure 
    toehold_mRNA_dict={}
    toehold_mRNA_dict['sequence']=toehold_mRNA
    toehold_mRNA_dict['structure']=mfe_structures[0].structure 
    toehold_mRNA_dict['paired']=mfe_structures[0].structure.pairlist()
    toehold_mRNA_dict['energy']=mfe_structures[0].energy
    
    binding_data.append(toehold_mRNA_dict['energy'])
    
    # Get the binding energy for the toehold on its own
    binding_data.append(toehold_dict['energy'])
            
    # Add the binding energy of the mRNA on its own
    binding_data.append(mRNA_binding_energy)

    # Look at the GC content
    gc_content = round(float(toehold_seq.count('G') + toehold_seq.count('C'))*100/len(toehold_seq), 2)
    binding_data.append(gc_content)
 
    
    return binding_data

def get_switch_recognition_seq(trigger):
    ''' This function receives a target trigger sequence and the type of molecule
    and obtains the RNA trigger for it.
    '''
    trigger_seq = Seq(trigger)
   
    return(trigger_seq.back_transcribe().reverse_complement().transcribe())


def write_toehold(seq_folder,toehold_seq,toehold_type,index):
    ''' Create a .txt file containing the toehold sequence.
    
    
    seq_folder (str): path where the sequence is contained
    
    toehold_seq (str): toehold sequence
    
    toehold_type (str): molecule type (DNA or RNA)
    
    index (int): toehold id 
   
    '''
        
    output_file=os.path.join(seq_folder, 'toehold_{!s}_{!s}.txt'.format(index,toehold_type))
    handle = open(output_file, 'w')
    handle.write('1\n')
    handle.write(str(toehold_seq) + '\n')
    handle.write('1')
    handle.close()


def get_full_switch(recognition_sequence, reporter, length_unpaired):
    ''' This sequence receives the recognition sequence and the reporter gene to design the 
    full toehold switch. 
    
    recognition_sequence : str
    
    reporter : str
    
    length_unpaired : int 
    
    
    return : str
    
    '''
    
    # Define the rest of the parts of the sequences
    part1 = 'GG'
    part2 = recognition_sequence
    
    # Get the reverse complement of the trigger to close the loop.
    rev_comp = str(part2[length_unpaired:].reverse_complement())
    
    # Parts 3 and 5 should be complementary
    # Adjust their length so that the hairpin has 19 total paired bases with 3 weak pairs at the top (AUA - UAU pairs)
    # Get the number of needed bases
    needed = 16 - (len(recognition_sequence) - length_unpaired)
    
    part3 = 'AUA' 
    part4 = 'CAGAAACAGAGGAGA' 
    part5 = 'UAU' 
    
    if needed <= 3:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'U'
        part6 = needed*'A'
        
        # Check how many bases of the complement of part2 should go before the AUG
        comp_pos = 3 - needed
        part6 = part6 + rev_comp[0:comp_pos] + 'AUG' + rev_comp[comp_pos+3:]
    elif needed <= 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'U'
        part6 = needed*'A'
        
        # Check where the complement of part2 should start
        comp_pos = 6 - needed
        part6 = part6 + 'AUG' + rev_comp[comp_pos:]
    elif needed > 6:
        # Add some weak base pairs and positions from the complement
        part2 = part2 + needed*'U'
        part6 = 'AAU' + 'AUG' + (needed-6)*'A' + rev_comp
    

    part7_comp = Seq(reporter[1:5]).back_transcribe().reverse_complement().transcribe()
    part7 = 'ACCUGGCGGCAG' + str(part7_comp)+ 'AAAG'
    
    
    # Test for stop codons and start codon
    stop_codons = ['UGA', 'UAA', 'UAG','AUG'] 
    test_region = part6[6:] + part7
    
    # Split the test region into codons to check one at a time
    test_region_codons = [test_region[i:i+3] for i in range(0, len(test_region), 3)]
                          
    # Loop through the list and check if there are any stop codons in this reading frame
    for codon in test_region_codons:
        # assert not codon in stop_codons, 'The generated switch contains a stop codon' 
        if codon in stop_codons:
            return('Stop')
                          
    part8 = reporter
    
    full_toehold = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8
    
    return(full_toehold)

def testprotein(protein,nb,suitable_aa):
    suitable=True
    if protein[0]!="M":
        #print(f'Amino Acid 1 = {protein[1]} - Protein do not begin by M\nSequence discarded \n ')
        suitable=False
    elif nb > 0:
        for i in range(1,nb+1):
            if suitable==True:
                if protein[i] not in suitable_aa:
                    #print(f'Amino Acid {i+1} = {protein[i]} - AA not favorable \nSequence discarded \n ')
                    suitable=False
    return suitable
        