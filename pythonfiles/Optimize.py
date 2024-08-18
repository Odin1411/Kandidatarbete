# -*- coding: utf-8 -*-
"""
Created on Mon May 27 18:14:09 2024

@author: olivi
"""

from Algorithm import Algorithm
from PhylogeneticTree import BinaryPhyloTree
import time
import numpy as np

def timealg():
    iqtree = "iqtree-2.3.2-Windows\Testing\oldtests\\"
    t = BinaryPhyloTree(iqtree +"tree1.nwk" ,iqtree + "test1.phy", "protein")
    
   # print(t.newicktree)
   # print(t.seqdata.alignment)
    #print(t.root.get_newick(features=["species","name"], format=1))
   # print(t.root)
    time.time()
    alg = Algorithm(["LG","WAG","cpREV"],"protein")
    
   
    #print(t.seqdata.alignment)   
  
    t2 = BinaryPhyloTree(iqtree +"tree2.nwk" ,iqtree + "align_tree2.phy", "nucleotide")
    starttime = time.time()
    alg.FPAinit(t,"WAG")
    print(time.time()-starttime)
    
timealg()
