# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:40:17 2024

@author: olivi
"""
from ete3 import Tree
import random
import math
import Substitutionmodels as subs
from pathlib import Path
from convert_results import convert_tree
from convert_results import convert_truetree
import math 
from numpy import random as np_rand
import numpy as np
from scipy.stats import poisson



def get_exponential_temperature(initialT,i,cf =0.995):
    T = initialT* (cf**i)
    if T<=1:
        return 1
    print(f'T: {T}')
    
    return T
    

class Probabilities:
    """a class for probabilities"""
    def __init__(self,model_probs):
        """model_probs is a dictionary with probabilities for respective model.
        If model_probs is a list, models will be chosen uniformly"""
        
        if type(model_probs) == list:
            self.model_probs ="uniform"
            self.models = model_probs
        
        else: 
            self.models = list(model_probs.keys())
            #make sure probabilities sum to 1
            assert sum([p for p in list(model_probs.values())]) == 1 
            self.model_probs = model_probs
     
    
        
    def get_model(self):
        """get a model according to the model probabilities"""
 
        model_probs = self.model_probs
        if model_probs == "uniform":
            return self.get_uniform()
        
        probabilities = sorted(model_probs.items(), key=lambda x: x[1])
    
        t=0
        u = random.uniform(0.0,1.0)
     
        for i in range(len(probabilities)):
            if t<u and u<=probabilities[i][1]+ t:
                return probabilities[i][0]
            t += probabilities[i][1]
            
        return probabilities[0][0]
        

    def get_uniform(self):
        """get model according to uniform distribution"""
        return random.choice(self.models)
 
p = Probabilities({"a":0.4,"b":0.1,"c":0.5})


       
       
class Algorithm:
    """Class for an algorithm
    you need to include all possible models in model_probs"""
    def __init__(self,model_probs,seqtype, tester_results = "../tester_results/"):
        assert seqtype =="nucleotide" or seqtype == "protein"
        self.probs = Probabilities(model_probs)
        self.submodel = subs.Substitutionmodel(seqtype)
        self.tester_results = tester_results
        
        
    def __str__(self):
        return f'{self.probs.model_probs}'
    
 
        
    def assign_branches_model(self,t,model=None):
        """
        Assign branches to models. If model not specified,
        assignment will generated using the modelprobs

        Parameters
        ----------
        t : BinaryPhyloTree
        model : name of a model
        
        Returns
            model assignment and model at the root"""
            
        assigned_models = {}
        if model:
            for node in t.root.traverse():
                assigned_models[node] = model
            model_at_root = model #assign rootchildren to the same model
            return assigned_models,model_at_root
        
        elif isinstance(self.probs.model_probs,dict):
            if 1 in self.probs.model_probs.values(): #if one model have the probability 1
                model = [k[0] for k in self.probs.model_probs.items() if k[1] ==1][0]
                for node in t.root.traverse():
                    assigned_models[node] = model
                model_at_root = model #assign rootchildren to the same model
                return assigned_models,model_at_root
                
        root_child = list(t.root.get_children())
        left_rootchild,right_rootchild = root_child[0],root_child[1]
        model_at_root = self.probs.get_model() #assign rootchildren to the same model
        
        assigned_models[left_rootchild] =  model_at_root
        assigned_models[right_rootchild] =  model_at_root  
        for subtree in [left_rootchild,right_rootchild]:  
            for v in subtree.iter_descendants():
                assigned_models[v] = self.probs.get_model()
        return assigned_models,model_at_root
    
       
    def get_model_at_root(self,t,assignment):
        root_child = list(t.root.get_children())
        left_rootchild,right_rootchild = root_child[0],root_child[1]
        assert(assignment[left_rootchild] == assignment[right_rootchild])
        return assignment[left_rootchild]
    
    def get_stationary_p(self,t,assignment):
        pi = self.get_model_at_root(t,assignment)
        return self.submodel.get_stationary(pi)
    
    def greedy_optimization(self,t):
        """a greedy method for finding the ML assignment"""
        
        prior_assignment,model_at_root = self.assign_branches_model(t)
        pi_at_root = self.submodel.get_stationary(model_at_root)
        Wp = self.FPArec(t,t.root,prior_assignment)
        gLLh_current = self.calculateGlobalLLh(Wp,pi_at_root, t)
        
        print(f'prior assignment: {prior_assignment}, gllh: {gLLh_current}')
        current_assignment = prior_assignment
        candidate_assignment = prior_assignment.copy()
        models = self.probs.models
        root_child = list(t.root.get_children())
        left_rootchild,right_rootchild = root_child[0],root_child[1]

        #check if model at root should be changed
        for candidate_root_model in models: 
            if current_assignment[left_rootchild] != candidate_root_model:
                candidate_assignment[left_rootchild] =  candidate_root_model
                candidate_assignment[right_rootchild] =  candidate_root_model
                candidate_pi = self.submodel.get_stationary(candidate_root_model)
                Wp = self.FPArec(t,t.root,candidate_assignment)
                gLLh_star = self.calculateGlobalLLh(Wp,candidate_pi, t)
              
                candidate_accepted = gLLh_current <  gLLh_star
                if candidate_accepted:
                    #without copy led to some strange behaviour
                   current_assignment = candidate_assignment.copy()
                   pi_at_root = candidate_pi
                   gLLh_current = gLLh_star
       
       
        candidate_assignment = current_assignment.copy()
        #check each branch in the subtree of right child and left child of root
        for subtree in [left_rootchild,right_rootchild]: 
            for v in subtree.iter_descendants():
                for m in models:             
                    if current_assignment[v] != m:
                        candidate_assignment[v] =  m
                       
                        Wp = self.FPArec(t,t.root,candidate_assignment)
                        gLLh_star = self.calculateGlobalLLh(Wp,pi_at_root, t)
                        print("")
                        print("running")
                      
                        candidate_accepted = gLLh_current <  gLLh_star
                        if candidate_accepted:
                           current_assignment = candidate_assignment.copy() 
                           gLLh_current = gLLh_star
                candidate_assignment = current_assignment.copy()
   
        print(current_assignment)
        print(gLLh_current)
        return current_assignment, gLLh_current, prior_assignment
                        
    
       
     
    def FPAinit(self,t,model=None): #we want to save which model is attached to which branch
        """Calculates global log likelihood using FPArec,"""
        assigned_models, model_at_root = self.assign_branches_model(t,model)
       
        pi_at_root = self.submodel.get_stationary(model_at_root)
        Wp = self.FPArec(t,t.root,assigned_models)
        gLLh = self.calculateGlobalLLh(Wp,pi_at_root, t)
      
        
        print(f'tree: {t.newicktree}')
        print(f'number of sites: {t.seqdata.n_sites},  type: {t.seqdata.seqtype}, alignmentfile: {t.seqdata.alignmentfile}')
       
        print(f'global log likelihood: {gLLh}')
        return gLLh
                
    
                
    def FPArec(self,t,v,assigned_models):
        """Performs recursive Felsensteins pruning algorithm
        
        Parameters
        ----------
        t : BinaryPhyloTree
        v : vertex object
        assigned_models : dictionary with vertices as keys and model names as values
        
        Returns
            Wp at root"""

        #m sites, k states, mxk matrix
        n_sites = t.seqdata.n_sites
        n_states = t.seqdata.n_states
        Wn = np.zeros((n_sites, n_states))
        
        if v.is_leaf(): 
            for i,s in enumerate(t.seqdata.getStates(v.name)): 
                Wn[i][s] = 1
            return Wn
        
        leftchild, rightchild = v.get_children()
        Wleft, Wright = self.FPArec(t,leftchild,assigned_models),self.FPArec(t,rightchild,assigned_models)
        bl_l,  bl_r= leftchild.dist, rightchild.dist
        modelleft, modelright = assigned_models[leftchild], assigned_models[rightchild]
        P_left, P_right = self.submodel.calculateP(modelleft, bl_l), self.submodel.calculateP(modelright, bl_r)
        
        for site in range(n_sites):
            Ll = np.dot(P_left, Wleft[site])
            Lr = np.dot(P_right, Wright[site])
            Wn[site] = Ll * Lr
            
        return Wn
              
        
        
    def calculateGlobalLLh(self, Wp,pi,t):
        """calculate the global log likelihood given W matrix at root, and the frequencies of 
        the nucleotides/proteins
        pi is stationary distribution"""
        
        n_sites = t.seqdata.n_sites
        globalLh = 0

        for s in range(n_sites):
            localLh = np.dot(Wp[s], pi)
            globalLh += np.log(localLh)   

        return globalLh
    
    
    def assign_descendants_model(self,candidate_assignment,v,model):
        """assign all descendants of vertex v to model model"""
        for u in v.iter_descendants():
             candidate_assignment[u] = model
        return candidate_assignment
    
    
    def propose_candidate_hybrid(self,t,current_assignment):
        if random.random()<=0.95:
            return self.propose_candidate(t,current_assignment)
        return self.propose_candidate_cladewise(t,current_assignment)
             
    
    def propose_candidate_cladewise(self,t,current_assignment):
        """a clade-wise proposal distribution"""
        nodes = []
        root_child = list(t.root.get_children())
        #if we added both root children, the probability of changing the model
        #at the root would be twice as large as any other node 
        left_rootchild,right_rootchild = root_child[0],root_child[1] 
        nodes.append(left_rootchild)
        for subtree in [left_rootchild,right_rootchild]:  
            for v in subtree.iter_descendants():
                nodes.append(v) 
                
        v = random.choice(nodes)
        models = self.probs.models.copy()
        models.remove(current_assignment[v])
        random_model = random.choice(models)
        candidate_assignment = current_assignment.copy()
       
        if v == left_rootchild:
            #assign the whole tree to a new model
            
            v = t.root
        else:
            candidate_assignment[v] = random_model
        
        self.assign_descendants_model(candidate_assignment,v,random_model)
        
        return candidate_assignment
        
        
        
    
    
    def propose_candidate2(self,t,current_assignment):
        """return a random node, random model"""
        nodes = []
        
        root_child = list(t.root.get_children())
       
        #if we added both root children, the probability of changing the model
        #at the root would be twice as large as any other node 
        left_rootchild,right_rootchild = root_child[0],root_child[1] 
        nodes.append(left_rootchild)
        nodes = list(left_rootchild.iter_descendants()) + list(right_rootchild.iter_descendants())
        nodes.append(left_rootchild)

      #  k_swaps = np_rand.binomial(n=len(nodes)-1, p=1/(3*len(nodes)), size=1)[0] +1 #maybe not depend on how many edges
        k_swaps = np_rand.binomial(n=5, p=1/35, size=1)[0] +1 #maybe not depend on how many edges

        print(f'k swaps: {k_swaps}')
       
        candidate_assignment = current_assignment.copy()
        for i in range(k_swaps):
            v = random.choice(nodes)
            nodes.remove(v)
            
            models = self.probs.models.copy()
            models.remove(current_assignment[v])
            mod = random.choice(models)
            candidate_assignment[v] = mod
            
            if v == left_rootchild: #make sure both root_children are the same model
                candidate_assignment[right_rootchild] = mod
                
  
        return candidate_assignment
    
    
    def propose_candidate(self,t,current_assignment):
        """return a random node, random model"""
        nodes = []
        root_child = list(t.root.get_children())
        #if we added both root children, the probability of changing the model
        #at the root would be twice as large as any other node 
        left_rootchild,right_rootchild = root_child[0],root_child[1] 
        nodes.append(left_rootchild)
        for subtree in [left_rootchild,right_rootchild]:  
            for v in subtree.iter_descendants():
                nodes.append(v) 
                
        v = random.choice(nodes)
        
        models = self.probs.models.copy()
        models.remove(current_assignment[v])
        random_model = random.choice(models)
        
        candidate_assignment = current_assignment.copy()
        candidate_assignment[v] = random_model
        
        if v == left_rootchild: #make sure both root_children are the same model
            candidate_assignment[right_rootchild] = random_model
        
        
        return candidate_assignment
       
    
    def simulated_annealing(self,t,m0,resultfile,initialT=50,N=50,testing=False):
        """
        performs a version of simulated annealing

        Parameters
        ----------
        t : BinaryPhyloTree
             object containing tree to be investigated and sequence data
        m0 : string
            null model (model when none specified)
        resultfile : string
            name of trace-file
        initialT : int
            initial temperature. The default is 50.
        N : int, optional
            number of iterations. The default is 50.
        testing : bool, optional
            if run should be executed as test (True) or an actual simulation (False). The default is False.

    
        """
        algname = "anne_"
        rf = resultfile
        if testing: 
           
            resultfile = "Tester_mh\\" +algname + rf+"_trace.txt" 
            pfile =  "Tester_mh\\" +algname + rf +"_proposals.txt" 
        else:   
            resultfile = "Results_mh\\" +algname + rf +"_trace.txt"
            file_path = Path(rf)
            if file_path.is_file():
                print(f'File "{resultfile}" already exists!')
                return
            
       
        f = open(resultfile,"w") 
        pf = open(pfile,"w") 
        
        prior_assignment = self.assign_branches_model(t)
        pi_at_root = self.get_stationary_p(t, prior_assignment)

        current_assignment = prior_assignment
        Wp = self.FPArec(t,t.root,current_assignment)
        current_gllh = self.calculateGlobalLLh(Wp,pi_at_root, t)
        print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0,"")}')
        print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0)}', file=f) 
        
        for i in range(0,N):
            candidate_assignment = self.propose_candidate_cladewise(t,current_assignment)
            candidate_pi = self.get_stationary_p(t, candidate_assignment)
            Wp = self.FPArec(t,t.root,candidate_assignment)
            candidate_gllh = self.calculateGlobalLLh(Wp,candidate_pi, t)
            print(f'{i} {current_gllh} {candidate_gllh} {convert_tree(t.root,candidate_assignment,m0,"")}', file=pf)
            print(f'candidate {candidate_gllh}')
            print(f'current_gllh {current_gllh}')
           
            log_accept_prob = min(0,candidate_gllh - current_gllh) 
            
            if log_accept_prob == 0:
                print("accepted!")
                current_assignment = candidate_assignment.copy()
                current_gllh = candidate_gllh
                pi_at_root = candidate_pi
            
            else:   
                u = random.random()
                T = get_exponential_temperature(initialT,i)
                
                print(f'u {math.log(u)}')
                print(f'accept prob { log_accept_prob/T}')
                print(f'candidate tree {convert_tree(t.root,candidate_assignment,m0,"")}')
                

                if math.log(u) < (log_accept_prob/T) :
                    print("accepted!")
                    current_assignment = candidate_assignment.copy()
                    current_gllh = candidate_gllh
                    pi_at_root = candidate_pi
                 
            print(f'{i} {current_gllh} {convert_tree(t.root,current_assignment,m0,"")}')
            print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0)}', file=f) 
    
        f.close()
        
    def prop(self,t,current_assignment,prop):
        if prop == "prop1":
            return self.propose_candidate(t,current_assignment)
        elif prop == "prop2":
            return self.propose_candidate2(t,current_assignment)
        elif prop == "hybrid":
            return self.propose_candidate_hybrid(t,current_assignment)
        elif prop == "clade":
            return self.propose_candidate_cladewise(t,current_assignment)
        
    def metropolis_hastings(self,t,m0,experiment_name,N=50,testing=True,local=False):
        """
        Perform metropolis hastings algorithm 

        Parameters
        ----------
        t : BinaryPhyloTree
        m0 : assumed model when none specified
        N : (int) number of iterations"""
      
       
        if testing: 
            resultfile = self.tester_results + experiment_name +"_trace.txt" 
            acceptfile =  self.tester_results + experiment_name +"_acceptratio.txt" 
        
        elif local:
            resultfile = "../results_mh/" + experiment_name +"_trace.txt"
            file_path = Path(resultfile)
            acceptfile = "../results_mh/" + experiment_name +"_acceptratio.txt"
            if file_path.is_file():
                print(f'File "{resultfile}" already exists!')
                return
            
        else:
            resultfile = "results_mh\\" + experiment_name +"_trace.txt"
            file_path = Path(resultfile)
            acceptfile = "results_mh\\" + experiment_name +"_acceptratio.txt"
            if file_path.is_file():
                print(f'File "{resultfile}" already exists!')
                return
       
        f = open(resultfile,"w")   
        af = open(acceptfile,"w")
        
        prior_assignment,model_at_root = self.assign_branches_model(t)
        pi_at_root = self.get_stationary_p(t, prior_assignment)
        current_assignment = prior_assignment
        Wp = self.FPArec(t,t.root,current_assignment)
        current_gllh = self.calculateGlobalLLh(Wp,pi_at_root, t)
      #  print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0,"")}')
        print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0)}', file=f) 
        
        for i in range(1,N+1):
            candidate_assignment = self.propose_candidate(t,current_assignment)
            candidate_pi = self.get_stationary_p(t,candidate_assignment)
            Wp = self.FPArec(t,t.root,candidate_assignment)
            candidate_gllh = self.calculateGlobalLLh(Wp,candidate_pi, t)
            log_accept_prob = min(0,candidate_gllh - current_gllh) 
            
            if log_accept_prob == 0:
                #candidate accepted
                print(1,file=af)
                current_assignment = candidate_assignment.copy()
                pi_at_root = candidate_pi
                current_gllh = candidate_gllh
            
            else:   
                u = random.random()
             
                if math.log(u) < log_accept_prob:
                    #candidate accepted
                    print(1,file=af)
                    current_assignment = candidate_assignment.copy()
                    pi_at_root = candidate_pi
                    current_gllh = candidate_gllh
                else:
                    #candidate rejected
                    print(0,file=af)
            print(" ")
            print(f'iteration {i}')
            print(f'{current_gllh} {convert_tree(t.root,current_assignment,m0)}', file=f) 
    
        f.close()
        af.close()
       
        
 

