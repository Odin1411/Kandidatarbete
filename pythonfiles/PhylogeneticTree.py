# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:55:45 2024

@author: olivi
"""
import math
import numpy
import os
from ete3 import Tree
from newick import loads
import Substitutionmodels as subs
import random



class Sequencedata:
    """class for sequencedata"""
    def __init__(self):
        "nucleotide,protein"
        self.nucleobases = {"A":0,"C":1,"G":2,"T":3}
        aminochars = "ARNDCQEGHILKMFPSTWYV"
        Aminostates = {}
        for i,a in enumerate(aminochars):
            Aminostates[a] = i
        self.aminoacids = Aminostates
          
    def set_seqtype(self,seqtype):
        """set the sequencetype"""
        assert seqtype == "nucleotide" or seqtype == "protein"
        self.seqtype  = seqtype
        if seqtype == "nucleotide":
            self.n_states = 4
        elif seqtype == "protein":
            self.n_states = 20    
            
        
    def set_n_sites(self,n_sites):
        self.n_sites = n_sites
        
    def set_alignment(self,alignment):
        """alignment is a dictionary with leaf-names as keys and sequences as values"""
        #alignment can look like {"a":"A", "b":"A","c":T} for example
        self.alignment = alignment

    def getStates(self,vname):
        """translate sequence of vname into corresponing state-number"""
        seq_v = self.alignment[vname]
        if self.seqtype == "nucleotide":
            return [self.get_nucleo_position(char) for char in seq_v]
        
        elif self.seqtype == "protein":
            return [self.get_amino_position(char) for char in seq_v]

        
    def get_amino_position(self,aminochar):
        """translate one char into a state-number"""
        return self.aminoacids[aminochar]

    def get_nucleo_position(self,nucleobase):
        """translate one char into a state-number"""
        return self.nucleobases[nucleobase]
    
    def read_alignment(self,alignmentfile,seqtype):
        #taxa has to be one word and each taxa needs a different name
        self.set_seqtype(seqtype)
        self.alignmentfile = alignmentfile.split("\\")[-1]
        
        with open(alignmentfile) as h:
            info = h.readline().split()
            n_taxa = int(info[0])
            self.set_n_sites(int(info[1]))
            alignment = {}
            for i in range(n_taxa): 
                taxa_seq =  h.readline().split()
            #   if len(taxa ==1) :
                   # seq = h.readline().split().strip()
                seq = taxa_seq[-1].strip()
                taxa = taxa_seq[0].strip()
                alignment[taxa] = seq
        self.set_alignment(alignment)
   
      
                


def clean_branchspecific(newicktree):
    """removes information in brackets, makes tree readible by library ete 3"""
    x = newicktree.count("[")
    tree = newicktree
    for i in range(x):
        i = tree.find("[")
        j = tree.find("]")
        tree = tree.replace(tree[i:j+1],"") 
    return tree
        
    

class BinaryPhyloTree:
    """A class for storing binary phylogenetic trees
    
        Attributes
    ----------
    root : ete 3 object, pointing at the root of the tree   
    treefile : file with tree
    seqdata : Sequencedata object. 
    newicktree : the contents of the tree-file 
    
    """
    def __init__(self,treefile,alignmentfile,seqtype):
         
         t = self.readTree(treefile)
         self.root = t 
         self.treefile = treefile.split("\\")[-1]
         self.leaves= [leaf for leaf in t]
         self.seqdata = Sequencedata()
         if alignmentfile:
             self.seqdata.read_alignment(alignmentfile,seqtype)
         self.assert_binary()
        
         
    def assert_binary(self):
        """assert that the tree is strict binary"""
        root = self.root 
        rc = list(root.get_children())
        assert(len(rc) == 2 or len(rc) == 0)
        for v in root.iter_descendants():
             vc = list(v.get_children())
             assert(len(vc) == 2 or len(vc) == 0)
  
    def readTree(self,file):
        """read a tree file and parse into a Tree object"""
        #"treename.nwk or treename.tree"
        with open(file) as h:
            tree = h.readline()
            self.newicktree = tree
            clean_tree = clean_branchspecific(tree)
        try: 
            t = Tree(clean_tree)
            return t
        except:
            return Tree(clean_tree, format= 1)
        
      
    def change_alignment(self,alignmentfile,seqtype):
        self.seqdata.read_alignment(alignmentfile,seqtype)
   
    def printTree(self):
        #tree = getTree(self.root)
        pass
       # print(loads(self.root).ascii_art())
 














