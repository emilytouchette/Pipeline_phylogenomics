# python script for my pipeline 
# add code here 

# Import needed modules
import os
import sys
import glob
import subprocess 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

# create string variables that store the input directory and output directory
indir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
outdir = "/scratch/touchete/BB485_Week06/Pipeline_phylogenomics/pipeline_out/"

# get a list of all input files
in_file_names= glob.glob(indir+"*fasta") 
#in_file_names=in_file_names[0:10]

### Run mafft using a “system call”. The following should go within a loop:

for file in in_file_names:
    new_file_path = file.replace(indir, outdir) # Create a new file path pointing to the output directory (this is how we tell mafft what to name the output)
    #print(new_file_path) #Check the file path
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path # Create a command string (this is what get called using the 'system call').
    print(aln_cmd)  # Check the command

    # Run the command below (uncomment) once you've double-checked that it's looking good.
    os.system(aln_cmd) 



### Run iqtree using a “system call” to perform tree inference. The following should go within a loop:
alignments= glob.glob(outdir+"*fasta") 

# Create the command. -nt 2 means two threads. If running this from within a job submission, you could use more threads to make it go faster.
for aln in alignments:
    tree_command = f"iqtree -s {aln} -m TEST -nt 12"
    print(tree_command) # Check the command 
    # Run the command using a 'system call'
    os.system(tree_command) # uncomment once you've check the command



trees= glob.glob(outdir+"*.treefile")

blank_list=[] 

for tree in trees:
    #Read in the tree and store as phylo object 
    temp_tree = Phylo.read(tree, "newick")
    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #Stope the loop once we found the correct tip
            break
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
    
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
        
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp=t 
        elif "Cr_" in t.name:
            Cr_temp=t
        elif "At_" in t.name:
            At_temp=t
        else:
            out_temp=t
            
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
        

    #Use series of if/else statemetns to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    #print(topo_str)

    blank_list.append(topo_str)

#print(blank_list)
n12top_count=blank_list.count("12top")
n13top_count=blank_list.count("13top")
n23top_count=blank_list.count("23top")
print(f"the number of trees with '12top' topology is {n12top_count}")
print(f"the number of trees with '13top' topology is {n13top_count}")
print(f"the number of trees with '23top' topology is {n23top_count}")

