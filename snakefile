import os
from ete3 import Tree


lengths = [200, 500, 800] #sequence lengths 

sys.path.append(os.path.abspath("pythonfiles"))
from initialize_trees import create_branch_specific_newick
from experiments import run_simulation
from experiments import run_mcmc


    
# Rule to simulate trees and alignments without wildcards
#rule all:
    #input:
       # expand("sim/sim{x}.{length}.phy", x=range(0,520), length=lengths)
        

# Rule to run greedy algorithm without wildcards
#rule all:
   # input:
       # expand("results/sim{x}.{length}.json", x=range(50,450), length=lengths)


# Rule to run mcmc without wildcards
rule all:
    input:
        expand("results_mh/sim{x}.{length}_{i}_trace.txt", x=range(0,50), i=range(1,5), length=lengths)




# Rule to generate a tree file
rule generateTree:
    output: "trees/tree{x}.tree"

    shell:
        '''
        iqtree-2.3.2-Windows\\bin\\iqtree2.exe -r 20 {output}
        '''

# Rule to root the tree 
rule rootTree:
    output: "trees/r_tree{x}.tree"
    input: "trees/tree{x}.tree"
    run:
        t = Tree(input[0], quoted_node_names=True, format=1)
        a = list(t.get_children())
        t.set_outgroup(a[0])
        rooted = t.write()
         
        with open(output[0], "w") as f:
            f.write(rooted)

        if int(wildcards[0])> 49:
            #add branch-specific models to .tree file
            create_branch_specific_newick(output[0], "WAG",["WAG","cpREV","mtArt","MtMAM","LG"]) 

       
        
        

# Rule to simulate amino acid sequences with different lengths
rule simAmino:
    output: "sim/sim{x}.{length}.phy"
    input: "trees/r_tree{x}.tree"
    run:
        out_base = os.path.splitext(output[0])[0]
       
        shell(f'''
            iqtree-2.3.2-Windows\\bin\\iqtree2.exe --alisim {out_base} -m WAG -t {input} --length {wildcards.length}
            ''')


# Rule to run greedy algorithm
rule greedySim:
    output: "results/sim{x}.{length}.json"
    input: "trees/r_tree{x}.tree", "sim/sim{x}.{length}.phy"
    run:
        treefile = input[0].replace("trees/","")
        aminofile = input[1].replace("sim/","")
        run_simulation(treefile ,aminofile,testing = False) 


rule mcmc:
    output: "results_mh/sim{x}.{length}_{i}_trace.txt"
    input: "trees/r_tree{x}.tree", "sim/sim{x}.{length}.phy"
    run:
        treefile = input[0].replace("trees/","")
        aminofile = input[1].replace("sim/","")
        run_mcmc(treefile ,aminofile,wildcards.i,testing = False) 


