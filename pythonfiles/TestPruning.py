# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 18:37:10 2024

@author: olivi
"""
from Algorithm import Algorithm
import matplotlib.pyplot as plt
from PhylogeneticTree import BinaryPhyloTree
import PhylogeneticTree as ptree
from ete3 import Tree
import random
from convert_results import convert_tree
from convert_results import convert_truetree
import time

from initialize_trees import check_homogeneity
from initialize_trees import create_node_model_map

class Test:
    def __init__(self,testdata_path = "../testdata/"):
        self.testdata_path = testdata_path

    def test_prob(self): 
        alg = Algorithm({"JC":0.5,"opt1":0.2,"opt2":0.3},"nucleotide")
        x = []
        n = 8000
        for i in range(n):
            x.append(alg.probs.get_model())
        
        print(sum( [1 for h in x if h == "JC"])/n)
        print(sum( [1 for h in x if h == "opt1"])/n)
        print(sum( [1 for h in x if h == "opt2"])/n)
        
        try:
            Algorithm({"JC":0.1,"opt1":0.2,"opt2":0.3},"nucleotide")
        except Exception as e:
            print("exception")
        alg = Algorithm(["WAG"],"protein")
        
        testdata = self.testdata_path
        
        t1 = BinaryPhyloTree(testdata +"r_tree1.tree" ,testdata + "test1_b.phy", "protein")
        
        assignment = alg.assign_branches_model(t1)
        print(assignment)
        alg = Algorithm({"WAG":0,"cpREV":1},"protein")
        assignment = alg.assign_branches_model(t1)
        print(assignment)
    
    
    def test_BPT():
        testdata = "../testdata/"
        
        t1 = BinaryPhyloTree(testdata +"tree1.nwk" ,testdata + "test1.phy", "nucleotide")
       # print(t.newicktree)
       # print(t.seqdata.alignment)
        #print(t.root.get_newick(features=["species","name"], format=1))
       # print(t.root)
        alg = Algorithm({"JC":1},"nucleotide")
        glh = alg.FPAinit(t1)
        t1.change_alignment(testdata+"align_tree1.phy","nucleotide")
        glh = alg.FPAinit(t1) 
        #print(t.seqdata.alignment)   
      
        t2 = BinaryPhyloTree(testdata +"tree2.nwk" ,testdata + "align_tree2.phy", "nucleotide")
        glh = alg.FPAinit(t2)
        
        t2.change_alignment(testdata + "WAG_tree2.phy", "protein")
        alg = Algorithm({"WAG":1},"protein")
        print(alg)
        alg.FPAinit(t2)
        
        alg = Algorithm({"LG":1},"protein")
        alg.FPAinit(t2)
        
        alg = Algorithm(["LG","WAG","cpREV"],"protein")
        
       # print(alg)
       # alg.FPAinit(t2)
    
        
    def test_convert(self):
        testdata = "../testdata/"
        t = BinaryPhyloTree(testdata +"r_tree0.tree" ,testdata + "test0_c.phy", "protein")
        t2 = BinaryPhyloTree(testdata +"r_tree1.tree" ,testdata + "test1_c.phy", "protein")
        #print(convert_truetree(t2.newicktree,"WAG"))
        assert(convert_truetree(t.newicktree,"WAG") == "((T2:0.088146[&&NHX:model=WAG],(T3:0.0278092[&&NHX:model=WAG],T5:0.00345282[&&NHX:model=WAG])1:0.286157[&&NHX:model=WAG])1:0.173125[&&NHX:model=WAG],(T1:0.175735[&&NHX:model=WAG],T4:0.0248948[&&NHX:model=WAG])1:0.173125[&&NHX:model=WAG]);")
        assert(convert_truetree(t2.newicktree,"WAG") == "((T2:0.088146[&&NHX:model=mtMam],(T3:0.0278092[&&NHX:model=cpREV],T5:0.00345282[&&NHX:model=cpREV])1:0.286157[&&NHX:model=cpREV])1:0.173125[&&NHX:model=WAG],(T1:0.175735[&&NHX:model=mtART],T4:0.0248948[&&NHX:model=WAG])1:0.173125[&&NHX:model=WAG]);")
        
        alg = Algorithm(['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'],"protein")
        assignment, rm = alg.assign_branches_model(t,"cpREV")
       
        assert(convert_tree(t.root,assignment,"WAG") == "((T2:0.088146[&&NHX:model=cpREV],(T3:0.0278092[&&NHX:model=cpREV],T5:0.00345282[&&NHX:model=cpREV]):0.286157[&&NHX:model=cpREV]):0.173125[&&NHX:model=cpREV],(T1:0.175735[&&NHX:model=cpREV],T4:0.0248948[&&NHX:model=cpREV]):0.173125[&&NHX:model=cpREV]);")
        assert(convert_tree(t.root,assignment,"WAG","") == "((T2:0.088146[cpREV],(T3:0.0278092[cpREV],T5:0.00345282[cpREV]):0.286157[cpREV]):0.173125[cpREV],(T1:0.175735[cpREV],T4:0.0248948[cpREV]):0.173125[cpREV]);")
        print("convert succeeded")
        
        
    def test_PT(self):
        testdata = "../testdata/"
        t = BinaryPhyloTree(testdata +"r_tree0.tree" ,testdata + "test0_c.phy", "protein")
        print(t.seqdata.alignmentfile)
        print(t.treefile)
   
    def test_branchspecific():
        
        testdata = "../testdata/"
        algH0 = Algorithm({"WAG":1},"protein")
        algH1 = Algorithm({"LG":0.3,"WAG":0.7},"protein")
       # print(alg)
       # alg.FPAinit(t2)
        tbs1 = BinaryPhyloTree(testdata +"aminotree_bs1.nwk" ,testdata + "aminotree_bs1.phy", "protein")
        print(tbs1.newicktree)
        algH0.FPAinit(tbs1)
        algH1.FPAinit(tbs1)
        
    
    def test_experiment():
        testdata = "../testdata/"
        experiment = "test0_c"
        t = BinaryPhyloTree(testdata +"r_tree0.tree" ,testdata + experiment +".phy", "protein")
      
        alg = Algorithm(['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'],"protein")
        
    def test_metropolishasting(self):
        testdata = self.testdata_path
        file = open("tidsinformation.txt","w")
        for i in ["200","500","800"]:
            for k in ["prop1","prop2","hybrid","clade"]:
                experiment = "sim0."+i
                t = BinaryPhyloTree(testdata +"r_tree0.20.tree" ,testdata + experiment +".phy", "protein")
                alg = Algorithm(['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'],"protein")
                starttime = time.time()       
                print(experiment)  
                alg.metropolis_hastings(t,"WAG",experiment,k,N=5000,testing = True)           
        
            print(f'total time: {time.time()-starttime} {experiment} {k}',file=file)
       
    
    def test_simulated_annealing():
        testdata ="../testdata/"
        experiment = "test1_c"
        t = BinaryPhyloTree(testdata +"r_tree1.tree" ,testdata + experiment +".phy", "protein")
      
        alg = Algorithm(['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'],"protein")
        starttime = time.time()       
            
        alg.simulated_annealing(t,"WAG",experiment+"_4",N=1000,testing = True)           
        
        print(f'total time: {time.time()-starttime}')
        
    
    def check_model_at_root(self,t,candidate_assignment,model_at_root):
        root_child = list(t.root.get_children())
        left_rootchild,right_rootchild = root_child[0],root_child[1]
        
        assert(candidate_assignment[left_rootchild] == candidate_assignment[right_rootchild])
        assert(candidate_assignment[left_rootchild] == model_at_root)
        
        
        
    def test_propose_candidate(self):
        testdata  = self.testdata_path
        experiment = "test0_b"
        t = BinaryPhyloTree(testdata +"r_tree0.tree" ,testdata + experiment +".phy", "protein")
      
        alg = Algorithm(['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'],"protein")
        starttime = time.time()       
        for i in range(5000):
            assignment,model_at_root = alg.assign_branches_model(t)
            self.check_model_at_root(t, assignment,model_at_root)
            print(convert_tree(t.root, assignment, "WAG"))
            print(model_at_root)
      
            candidate_root, proposal = alg.propose_candidate(t,assignment,model_at_root) 
            self.check_model_at_root(t,proposal,candidate_root)
            
            print(convert_tree(t.root, proposal, "WAG"))
            print(candidate_root)
            
            candidate_root, proposal = alg.propose_candidate_cladewise(t, assignment, model_at_root)
            self.check_model_at_root(t,proposal,candidate_root)
            
            print(convert_tree(t.root, proposal, "WAG"))
            print(candidate_root)
        
        #print(convert_tree(t.root, proposal2, "WAG"))
      
        
    def test_initalize(self):
        testdata ="../testdata/"
        homo = {"1":"WAG","2":"WAG","3":"WAG","4":"WAG","5":"WAG"}
        hetero1 = {"1":"cpREV","2":"WAG","3":"WAG","4":"WAG","5":"WAG"}
        hetero2 = {"1":"cpREV","2":"WAG","3":"cpRev","4":"WAG","5":"WAG"}
        assert (check_homogeneity(homo) == True)
        assert (check_homogeneity(hetero1) == False)
        assert (check_homogeneity(hetero2) == False)
        tree = Tree(testdata +"r_tree0.tree", quoted_node_names=True, format=1) 
        for i in range(100):
            k = random.choice([3,4,5])
            m = create_node_model_map(tree, "WAG", ["WAG","cpREV","mtArt","MtMAM","LG"],k)
         
            assert check_homogeneity(m) == False
        
        
def test_all():
         test = Test() 
         test.test_initalize()
        # test.test_propose_candidate()
       #  test.test_PT()
        # test.test_metropolishasting()
        # test.test_convert()

test_all()