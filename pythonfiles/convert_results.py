# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:16:12 2024

@author: olivi
"""
import matplotlib.pyplot as plt



def add_m0(char,newicktree,model0,modstr = "&&NHX:model="):
    """add nullmodel to the tree"""
    indices = []
    for i in range(len(newicktree)):
        if newicktree[i] == char:
            indices.append(i)
    point = 0
    for i in indices:
        i = i + point
        if newicktree[i-1] !="]":
            newicktree = newicktree[:i] +f'[{modstr}{model0}]' + newicktree[i:]
            point +=len(f'[{modstr}{model0}]')
 
    return newicktree


def convert_truetree(newicktree,model0,modstr ="&&NHX:model="):
    """converts the newicktree from the tree-file to a format readable by ggtree in r"""

    tree = newicktree
    i =  tree.find("[")
    if i!=-1: #if there are models specified in the tree
        #append modstr in brackets
        newicktree = newicktree.replace(tree[i+1:i+8],modstr) 
    newicktree = add_m0(")",newicktree,model0) 
    newicktree = add_m0(",",newicktree,model0) 
    return newicktree


def convert_tree(v,assignment,model0,modstr="&&NHX:model=",add_m0 = True):
    """converts the tree assignment to a newicktree 
    OBS: only works for fully/strictly binary trees
        Parameters:
        v : root node (first iteration)
        assignment : dictionary with nodes as keys and modelnames as values
        model0 : the model used when none specified
        modstr : 
        add_m0 : (bool) if True add model0 in brackets
    """
    
    if not v.is_root():
        model = assignment[v]
  
    if v.is_leaf():  
        s = f'{v.name}:{v.dist}[{modstr}{model}]'
        
        if model == model0 and not add_m0: 
            
            s = f'{v.name}:{v.dist}'
        return s
    
    leftchild,rightchild = list(v.get_children())[0], list(v.get_children())[1]
    
    if v.is_root(): #the root will not have a distance or model
        s=";"
    else:
        s = f'{v.name}:{v.dist}[{modstr}{model}]'
        if model == model0 and not add_m0:
           
            s = f'{v.name}:{v.dist}'
    
    return f"({convert_tree(leftchild,assignment,model0,modstr,add_m0)},{convert_tree(rightchild,assignment,model0,modstr,add_m0)}){s}"
    

