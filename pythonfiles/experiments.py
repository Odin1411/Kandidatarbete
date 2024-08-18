# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 14:25:28 2024

@author: olivi
"""

from Algorithm import Algorithm
import matplotlib.pyplot as plt
from PhylogeneticTree import BinaryPhyloTree
import time 
import json
from pathlib import Path
from convert_results import convert_tree
from convert_results import convert_truetree
import os

    
def run_h0(t,models,json_save, filename):
  
    with open(filename,"w") as f:
        f.write("\n")
        f.write(f'tree: {t.newicktree}\n')
        f.write(f'number of sites: {t.seqdata.n_sites} \ntype: {t.seqdata.seqtype} \nalignmentfile: {t.seqdata.alignmentfile}\n')
        json_save["h0"] = {"model":[], "gllh":[]}
        for m in models:
            alg0 = Algorithm({m:1},"protein")
            gllh0 = alg0.FPAinit(t)
            json_save["h0"]["model"].append(m)
            json_save["h0"]["gllh"].append(gllh0)
            f.write(f'{alg0}\n')
            f.write(f'global log likelihood: {gllh0}\n')



def save_experiment(t,alg_star, ml_assign,gllh,filename,model0,json_save,dt = None):

    with open(filename,"a") as f:
    
        f.write("\n \n")
        f.write(f'Experiment: {alg_star}\n')
       # f.write(f'ML assignment: {ml_assign}\n')
        converted_tree = convert_tree(t.root,ml_assign,model0)
        f.write(f'ML assignment: {converted_tree}\n')
        f.write(f'glh: {gllh}\n')
        if dt:
            f.write(f'time to execute: {dt}')
            
        treefile = filename.replace(".txt",".tree")   
        
        json_save["priorprobs"].append(str(alg_star.probs.model_probs))
       
        json_save["mltree"].append(converted_tree)
        json_save["gllh"].append(gllh)
        if dt:
            json_save["time"].append(dt)
        
    

def run_experiment(t,alg_star, model0,json_save,filename=None):
    print("\n")
    print(f'tree: {t.newicktree}\n')
    print(f'number of sites: {t.seqdata.n_sites} \ntype: {t.seqdata.seqtype} \nalignmentfile: {t.seqdata.alignmentfile}\n')
    
    start_time = time.time()
    ml_assign, gllh, prior= alg_star.greedy_optimization(t)
    dt = time.time() - start_time
    
    if filename:
        json_save["prior"].append(convert_tree(t.root,prior,model0))
        save_experiment(t,alg_star,ml_assign,gllh,filename,model0,json_save,dt)
    


def run_different_priors(t,priors,models,seqtype,model0,filename,runH0 = True ):
    start_time = time.time()
    
    json_save = {"prior": [], "priorprobs":[], "mltree": [], "gllh":[] ,"time": []}
    json_save["truetree"] = convert_truetree(t.newicktree, model0)
    json_save["n_sites"] = t.seqdata.n_sites
    json_save["seqtype"] = t.seqdata.seqtype
    json_save["n_leaves"] = len(t.leaves)
    json_save["alignmentfile"]= t.seqdata.alignmentfile
    
    if runH0:
        run_h0(t,models,json_save,filename)
    
    for prior in priors:
        print("\nnew run\n")
        alg = Algorithm(prior, seqtype)
        run_experiment(t, alg,model0,json_save,filename)
   
    dt = time.time() - start_time 
    print(dt)
    f = open(filename,"a")
    f.write("\n \n")
    f.write(f'total time to execute: {dt}')
    f.close()
    jsonfile = filename.replace(".txt",".json") 
    
    json_save["totaltime"] = dt
    print(jsonfile)
    print(json_save)
    
    with open(jsonfile , "w") as f:
        json.dump(json_save, f)
        
    


def run_simulation(treefile,aminofile, testing = True):
    
    if testing:
        run = "_1"
        #path =  os.path.normpath(os.getcwd() + os.sep + os.pardir)
        #iqtree = "..\iqtree-2.3.2-Windows\Testing\\"
        testdata = "testdata"
        resultfile = "\\tester_results\\" + aminofile.replace('.phy', '') + run + ".txt" #result saved to Tester_results, Results_2
        t1 = BinaryPhyloTree(testdata+treefile ,testdata+aminofile, "protein")
    
    else: 
       # path =  os.path.normpath(os.getcwd() + os.sep + os.pardir)
        resultfile ="results\\" + aminofile.replace('.phy', '') +".txt"
        t1 = BinaryPhyloTree("trees\\" +treefile ,"sim\\" +aminofile, "protein")
    
    file_path = Path(resultfile)

    if file_path.is_file():
        print(f'File "{resultfile}" already exists!')
        return
   
        
    model0 = "WAG"
    
    pr1 = {"LG":0.4,"WAG":0.3,"cpREV":0.15,"mtArt":0.15,"MtMAM":0}
    pr2 = {"WAG":0.1,"cpREV":0.2,"mtArt":0.6,"LG":0.1,"MtMAM":0}
    pr3 = {"cpREV":0,"mtArt":0,"MtMAM":0,"LG":0,"WAG":1}
    pr4 = {"cpREV":1,"mtArt":0,"MtMAM":0,"LG":0,"WAG":0}
    pr5 = {"cpREV":0,"mtArt":1,"MtMAM":0,"LG":0,"WAG":0}
    pr6 =  {"cpREV":0,"mtArt":0,"MtMAM":1,"LG":0,"WAG":0}
    pr7 = ["WAG","cpREV","mtArt","MtMAM","LG"]
    pr8 = ["WAG","cpREV","mtArt","MtMAM","LG"]
    pr9 = ["WAG","cpREV","mtArt","MtMAM","LG"]
    pr10 = ["WAG","cpREV","mtArt","MtMAM","LG"]
    priors = [pr1,pr3,pr2,pr4,pr5,pr6,pr7,pr8,pr9,pr10]
    
    run_different_priors(t1,priors,pr7,"protein", model0,resultfile,runH0 = True)

    #t2 = BinaryPhyloTree(iqtree +"aminotree_mtART.nwk" ,iqtree + "mtART_tree.phy", "protein")
    #glh0 = alg0.FPAinit(t2)
    
    #alg = Algorithm({"WAG":0.5,"mtArt":0.5},"protein")
    
    #glh0 = alg0.FPAinit(t2) 
    #ml_assign, gllh = alg.greedy_optimization(t2)
    #run_experiment(t2, alg0, alg, "Results\experiment2.txt")

   
def test_run():
    iqtree = "iqtree-2.3.2-Windows\Simulations\\"
    t1 = BinaryPhyloTree(iqtree +"aminotree_LG.nwk" ,iqtree + "LG_tree.phy", "protein")
    pr1 = {"LG":1,"WAG":0}
    alg = Algorithm(pr1, "protein")
    
    run_experiment(t1, alg)
    

def run_mcmc(treefile,aminofile,chain_i,testing):
     start_time = time.time()
     json_save = {"truetree": [], "totaltime":[]}
     experiment_name = aminofile.replace('.phy', '')  +"_" + chain_i
     t1 = BinaryPhyloTree("trees\\" +treefile ,"sim\\" +aminofile, "protein")
     alg = Algorithm(["WAG","cpREV","mtArt","MtMAM","LG"],"protein")
     alg.metropolis_hastings(t1,"WAG",experiment_name,N = 2000,testing = testing)
     
     dt = time.time() - start_time 
     jsonfile = experiment_name + ".json"
     
     json_save["truetree"] = convert_truetree(t1.newicktree, model0 = "WAG")
     json_save["totaltime"] = dt
    
     with open(jsonfile , "w") as f:
         json.dump(json_save, f)   
   
    
def create_json():
    respath = os.path.normpath(os.getcwd() + os.sep + os.pardir) + "\\results_mh"
    back = os.path.normpath(os.getcwd() + os.sep + os.pardir) 
   
    for L in [200,500]:
        
        for i in range(0,50):
            for chain in range(1,5):
                file_name = respath + f'\\sim{i}.{L}_{chain}_trace.txt'
                file_path = Path(file_name)
               
                if file_path.is_file():
                    print(f'File "{file_path}" already exists!')
                   
                    jsonfile = respath + f'\\sim{i}.{L}_{chain}.json'
                    
                    treefile =  f'r_tree{i}.tree'
                    aminofile = f'sim{i}.{L}.phy'
                    t = BinaryPhyloTree(back + "\\trees\\" +treefile , back + "\\sim\\" +aminofile, "protein")
                    print(treefile)
                    print(aminofile)
                   
                    json_save = {"truetree": [], "totaltime":[]}
                    json_save["truetree"] = convert_truetree(t.newicktree, "WAG")
                    json_save["totaltime"] = None
                    
                    
                    with open(jsonfile , "w") as f:
                        json.dump(json_save, f)

def create_json_local():
    respath = os.path.normpath(os.getcwd() + os.sep + os.pardir) + "/results_mh"
    back = os.path.normpath(os.getcwd() + os.sep + os.pardir) 
   
    for L in [200,500]:
        
        for i in range(0,50):
            for chain in range(1,5):
                file_name = respath + f'/sim{i}.{L}_{chain}_trace.txt'
                file_path = Path(file_name)
               
                if file_path.is_file():
                    print(f'File "{file_path}" already exists!')
                   
                    jsonfile = respath + f'/sim{i}.{L}_{chain}.json'
                    
                    treefile =  f'r_tree{i}.tree'
                    aminofile = f'sim{i}.{L}.phy'
                    t = BinaryPhyloTree(back + "/trees/" +treefile , back + "/sim/" +aminofile, "protein")
                    print(treefile)
                    print(aminofile)
                   
                    json_save = {"truetree": [], "totaltime":[]}
                    json_save["truetree"] = convert_truetree(t.newicktree, "WAG")
                    json_save["totaltime"] = None
                    
                    
                    with open(jsonfile , "w") as f:
                        json.dump(json_save, f)
                        
                        
create_json_local()


def run_mcmc_local(treefile,aminofile,chain_i):
     start_time = time.time()
     json_save = {"truetree": [], "totaltime":[]}
     experiment_name = aminofile.replace('.phy', '')  +"_" + chain_i
     t1 = BinaryPhyloTree("../trees/" +treefile ,"../sim/" +aminofile, "protein")
     alg = Algorithm(["WAG","cpREV","mtArt","MtMAM","LG"],"protein")
     alg.metropolis_hastings(t1,"WAG",experiment_name,N = 2000,testing=False, local =True)
     
     dt = time.time() - start_time 
     jsonfile = "../results_mh" + experiment_name + ".json"
     
     json_save["truetree"] = convert_truetree(t1.newicktree, model0 = "WAG")
     json_save["totaltime"] = dt
    
     with open(jsonfile , "w") as f:
         json.dump(json_save, f)   
           
    
if __name__ == "__main__":
   # create_json()
    #test_run()    
    #for i in range(1,5):
     #   run_mcmc_local("r_tree51.tree","sim51.500.phy",str(i))
     #   run_mcmc_local("r_tree54.tree","sim35.500.phy",str(i))
   # for j in range(50,55):
    #    for i in range(1,5):
     #       run_mcmc_local(f"r_tree{j}.tree",f"sim{j}.200.phy",str(i))
          
    for j in range(59,62):
         for i in range(1,5):
              run_mcmc_local(f"r_tree{j}.tree",f"sim{j}.500.phy",str(i))
              
    for j in range(59,62):
         for i in range(1,5):
              run_mcmc_local(f"r_tree{j}.tree",f"sim{j}.200.phy",str(i))
   # run_simulation(  "r_tree1.tree" ,"test1_b.phy" ,testing = True) 
    pass
  #  create_json()

#k = {'prior': [['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG']], 'mltree': ['(T1:0.0105987[&&NHX:model=WAG],(T2:0.349701[&&NHX:model=WAG],(T3:0.0264179[&&NHX:model=WAG],(T4:0.175924[&&NHX:model=cpREV],T5:0.00922531[&&NHX:model=cpREV]):0.021873[&&NHX:model=cpREV]):0.0370144[&&NHX:model=WAG]):0.0105987[&&NHX:model=WAG])'], 'gllh': [-5690.388416438044], 'time': [154.30851888656616], 'tree': '(T1:0.0105987[&model=mtART],(T2:0.349701,(T3:0.0264179,(T4:0.175924[&model=cpREV],T5:0.00922531[&model=cpREV])1:0.021873[&model=cpREV])1:0.0370144)1:0.0105987);', 'n_sites': 1000, 'seqtype': 'protein', 'alignmentfile': 'r_5.phy', 'h0': {'model': ['WAG', 'cpREV', 'mtArt', 'MtMAM', 'LG'], 'gllh': [-5713.530877828688, -5750.605697386066, -6291.933832040931, -6316.922222332191, -5770.735009641532]}} 
    
 
# t1.change_alignment(iqtree+"align_tree1.phy","nucleotide")
  







