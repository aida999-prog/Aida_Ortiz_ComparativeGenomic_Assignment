import pandas as pd
from re import sub
from collections import Counter
import numpy as np
from itertools import zip_longest

# This script uses the Multiple Sequence Alignment file which contains the aligned homologs of the UvrABC system protein A of the 25 species

with open("mult_align_25sp", 'r') as f:
        lines = f.readlines()
        char = "*"
        alignm = lines[3:len(lines)] # the first lines are not useful
        for line in alignm: 
            parts = line.split()
            if len(parts)>1 and char not in parts: # because there still exists some blank lines within the file
                aligned_seqs = parts[1]
                if char not in aligned_seqs: # there were still remaining "*" characters
                    print(aligned_seqs, file=open("alignment.txt","a")) # printing directly to file just the aligned sequences

# ------------------

with open("alignment.txt", 'r') as f:
    lines = f.readlines()
    for line in lines:
        line_mod = " ".join(line) # separating each character with a black space to be able to read it as a df
        print(line_mod, file=open("alignment2.txt","a")) 

alignm_in_table = pd.read_fwf("alignment2.txt", sep = " ", header=None)

nb_rows = len(alignm_in_table.index)
i = 0
while i<nb_rows:
    alignm = alignm_in_table.iloc[i:i+25,:] # because of the format of the alignment file each 25 rows the respective continuation of each sequence happen
    data_new1 = alignm.replace(r'^s*$', float('NaN'), regex = True) # converting blank spaces into NaN values
    #print(data_new1)
    df = data_new1.dropna(axis=1) # drops COLUMNS with NaN values
    #print(df)
    cols = df.nunique() # nb of unique characters of each column
    print(cols, file=open("k.txt","a"))
    
    i = i+25

# --------------------

with open("k.txt", 'r+') as f:
    string = "int64"   
    lines = f.readlines()
    for line in lines: 
        parts = line.split()
        if string not in parts: # removing unuseful information
            vals = parts[1]
            print(vals, file = open("k_values.txt","a"))

df = pd.read_fwf("alignment2.txt", sep = " ", header=None)
nb_rows = len(df.index)

i = 0

while i<nb_rows:
    alignm = df.iloc[i:i+25,:]

    for column in alignm.columns:
        col = alignm.iloc[:,column]

        frequency_dict = {} # creating dictionary with the different characters present on each column as key and its frequency as value
        for character in col:
            if character in frequency_dict:
                frequency_dict[character] += 1
            else:
                frequency_dict[character] = 1
        
        foo = Counter(frequency_dict)
        max_val = foo.most_common()[0] # retrieveng the most common one of each column with its value
        
        print(max_val, file = open("n.txt","a"))
        
    i = i+25

# ---------------------------

with open("n.txt", 'r+') as f:
    string = "nan"   
    lines = f.readlines()
    
    for line in lines:
        if string not in line: # removing na values
            output = re.sub(r"[\([{})\]]", "", line)
            parts = output.split()
            n_value = parts[1] # just selecting its respective frequencies
            
            if string not in parts[0]: # removing remaining na values
                print(n_value, file = open("n_values.txt","a"))

# --------------------

with open("k_and_n_values.txt", 'w') as res, open("k_values.txt") as f1, open("n_values.txt") as f2:
    for line1, line2 in zip_longest(f1, f2, fillvalue=""):
        res.write("{} {}\n".format(line1.rstrip(), line2.rstrip())) # joining both k and n values in a file

# each ROW of the FILE correspond to the values for a COLUMN of the ALIGNMENT

# --------------------

df = pd.read_fwf("k_and_n_values.txt", sep = " ", header=None)
df.rename(columns={0:"k_values", 1:"n_values"}, inplace=True) # renaming the columns to clarify

N = 25 # number of sequences in the alignment
df["Variability"] = (N*df["k_values"])/(df["n_values"]) # calculiting variability and adding as a new column to the df
print(df.to_string(), file = open("variability_values.txt","a"))

# conserved regions would be the positions in the alignments with just one type of character in all the sequences at that position (k=1, n=25)
# the first column represents the position in the alignment
variability = df["Variability"]
conserved_regions = df.query('Variability==1.0') # selecting the regions which variability values are bellow the threshold established
print(conserved_regions.to_string(), file = open("conserved_regions.txt","a"))