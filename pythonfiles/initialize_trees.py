# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:21:30 2024

@author: olivi
"""

from ete3 import Tree
from newick import loads
from Algorithm import Probabilities
import random
from convert_results import convert_tree


def root_tree(unroot_treefile, iqtree= "iqtree-2.3.2-Windows\Random\\"):
    """first step in the process, roots a binary tree. 
    Used on a generated trees and adds prefix r_"""

    tree = Tree(iqtree + unroot_treefile, quoted_node_names=True, format=1)
   
    a = list(tree.get_children())
    tree.set_outgroup(a[0])
    rooted = tree.write()
    
    with open(iqtree + "r_" + unroot_treefile, "w") as f:
        f.write(rooted)
    
    print(Tree(iqtree + "r_" + unroot_treefile))
  
    
def check_homogeneity(model_tree_map):
    """check if model assignment is homogeneous"""
    models = list(model_tree_map.values())
    mod = models[0]
    if models.count(mod) == len(models):
        return True
    return False
        
def create_node_model_map(tree, model0,modelprobs, iterations):
     """create a dictionary with node as key and model as value
     """
     nodes = []
     root_child = list(tree.get_children())
     #if we added both root children, the probability of changing the model
     #at the root would be twice as large as any other node 
     left_rootchild,right_rootchild = root_child[0],root_child[1] 
     nodes.append(left_rootchild)
     nodes = list(left_rootchild.iter_descendants()) + list(right_rootchild.iter_descendants())
     nodes.append(left_rootchild)
     
     nodes_model_map = {node:model0 for node in nodes + [right_rootchild]}
     #make sure name format is compatible with iqtree
     translate_modname = {"MTART":"mtART","MTMAM":"mtMam", "WAG":"WAG","LG":"LG","CPREV":"cpREV"}
     ModelProb = Probabilities(modelprobs)
     i = 0 
     while i < iterations: 
      
         node = random.choice(nodes)
         model = ModelProb.get_model()
         model = translate_modname[model.upper()] 
         #if node chosen is a root-child then set node to the root, iterate through roots decendants
         if node == left_rootchild:
             node = tree
             
         for desc in node.iter_descendants():
             nodes_model_map[desc] = model
        
         nodes_model_map[node] = model
         if i == iterations-1:
             if check_homogeneity(nodes_model_map): #make sure that there are more than one model in the tree
                 i = i-1
         i = i+1
     return nodes_model_map
 

def create_branch_specific_newick(treefile, model0,modelprobs, iterations=None):
    """
    adds heterogeneity to the treefile

    Parameters
    ----------
    treefile : .tree file
    model0 : Model used when none specified
    modelprobs : Possible models in the tree 
    iterations : (int) number of iterations
        

    """
    try:
        tree = Tree(treefile, quoted_node_names=True, format=1) 
        if not iterations:
            iterations = random.choice([3,4,5])
        nodes_model_map = create_node_model_map(tree, model0, modelprobs,iterations)
        #convert to a newick format readable by iqtree
        mod_newick = convert_tree(tree,nodes_model_map,model0,modstr="&model=",add_m0 = False)
        print(mod_newick)
        with open(treefile, "w") as f:
           f.write(mod_newick)
    except:
        print("Error: tree already has branch-specific models")
        return
    
    
#create_branch_specific_newick(iqtree + tree, "WAG", ["WAG","cpREV","mtArt","MtMAM","LG"],3)
    



        
     
     
     

