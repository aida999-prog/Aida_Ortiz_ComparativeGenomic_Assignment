# This function retrieves a txt file with the most significant match for each sequence of the desired species.
# Both input names should be between quotation marks ("__")
# "name_for_the_outputallmatches_file" should be a .txt file name

def best_hits_of_one_species(name_of_the_blast_result_file, name_for_the_outputallmatches_file):  
    target_line="Sequences producing significant alignments:" # target string to find the matches and their E-values
    value = 0
    match = " "
    output_file = open(name_for_the_outputallmatches_file,"a") # new file where 
    with open(name_of_the_blast_result_file, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if target_line in line: # if it contains the string
                line_score=lines[index+2] # selecting the second following row after the match because the first one is a blank line
                parts = line_score.split() # split line into parts
                score = float(parts[-2]) # selecting the E-value and converting it to float
                e_value = float(parts[-1])
                #if (e_value < 0.0): # E-value threhold
                if score > value:
                    value = score # selecting the max
                    match = line_score
                    output_file.write(line_score) # saving the best matches into the new txt file

# This function compares both output files from both species (the ones containing each most significant match for each sequence)
# and looks for the match between them with the higuest E-value. This would correspond to the Best Bidirectional Hit (BBH).
# Both input names should be between quotation marks ("__")

def BBH_of_both_species(one_especies_matches_file, other_species_matches_file):
    with open(one_especies_matches_file, "r") as f1,open(other_species_matches_file, "r") as f2:
        list1 = f1.readlines()
        list2 = f2.readlines()

        for line_1 in list1:
            for line_2 in list2:
                parts1 = line_1.split()
                name1 = parts1[1]
                whole_name1 = parts1[1:-2] # name of the protein
                #print(whole_name1)
                parts2 = line_2.split()
                name2 = parts2[1]
                whole_name2 = parts2[1:-2] # name of the protein
                #print(whole_name2)
                if name2==name1:  # if it is bidirectional
                    name_of_the_BBH = '  '.join(whole_name2)
                    print(name_of_the_BBH, file=open("BBH.txt","a")) # print it to file


best_hits_of_one_species("thermotoga_vs_thermophilus.tsv","matches_thermotoga_thermophilus.txt")
best_hits_of_one_species("thermophilus_vs_thermotoga.tsv", "matches_thermophilus_vs_thermotoga.txt")

BBH_of_both_species("matches_thermotoga_thermophilus.txt", "matches_thermophilus_vs_thermotoga.txt")