# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:35:16 2024

@author: olivi
"""

import os
import os.path
# importing libraries
import seaborn as sns
import matplotlib.pyplot as plt
("../results_mh")


def check_lengths(length):
    """check that the trace files are of correct length"""
    n_files = 0
    path = "../results_mh/"
    files =[]
    
    for j in range(50,100):
        k = 0
        for i in range(1,5):
            tracefile = f'sim{j}.{length}_{i}_trace.txt'
            gllh = []
            x = []
            e = os.path.isfile(path+tracefile)
            if e:
                n_files = n_files +1
                with open(path + tracefile) as h:
                    for i,line in enumerate(h):
                        l = line.strip()
                        
                        try:
                            gllh.append(float(l.split()[0]))
                            x.append(i)
                        except ValueError:
                            print(l.split()[0])
                            pass
                       # print(l.split())
                try:
                     assert((len(gllh))==2001) 
                    
                except:
                
                    print(f'failed for file {tracefile}')
            else:
               print(f'{tracefile} not found')
            
    print(n_files)
    print("should be: 200 files")
          
 


def check_trace_exists(length):
    files =[]
    
    for j in range(50,100):
        k = 0
        for i in range(1,5):
            file = f'sim{j}.{length}_{i}_trace.txt'
            jsonfile = f'sim{j}.{length}_{i}.json'
            e = os.path.isfile("../results_mh/"+file)
            if e: #and os.path.isfile("../results_mh/"+jsonfile):
                k = k+1
             #   print(jsonfile)
       # print(k)
        files.append(k)
    print(files)
    print(files.count(4))
    
    


        

def plot_density(tracefile,ylim_a=0,ylim_b=None): 
    path = "../results_mh/"
   
    tot_gllh = []
    for i in range(1,5):
        f = f'{path}{tracefile}_{i}_trace.txt'
        gllh = []
        x = []
        with open(f) as h:
            for i,line in enumerate(h):
                l = line.strip()
            
                try:
                    gllh.append(float(l.split()[0]))
                    x.append(i)
                except ValueError:
                    print(l.split()[0])
                    pass
       
        print(len(gllh[-1000:]))
        clean_gllh = gllh[-1000:]
        tot_gllh = tot_gllh + clean_gllh
       
        
           # print(l.split())
    
    sns.distplot(a=tot_gllh, hist=False)
    #plt.plot(x, gllh)
   #plt.ylim(-6527,-6510)
  
#plot_density("sim212.800")     

    
def read_trace(tracefile,ylim_a=0,ylim_b=None): 
    """creates a trace plot"""
    path = "../results_mh/"
    gllh = []
    x = []
    with open(path + tracefile) as h:
        for i,line in enumerate(h):
            l = line.strip()
            
            try:
                gllh.append(float(l.split()[0]))
                x.append(i)
            except ValueError:
                print(l.split()[0])
                pass
           # print(l.split())
    
    
    plt.plot(x, gllh)
   #plt.ylim(-6527,-6510)
    if ylim_b:
       plt.ylim(ylim_a,ylim_b)
    plt.axvline(x=1000, color = "red")
    plt.xlabel("Iteration")
  #  plt.ylabel("Log likelihood")
   # plt.savefig(tracefile+'.png')
    plt.show()


            
def get_acceptratio(acceptfile):
    """calculate acceptratio of a file (with burn-in)"""
    path = "../tester_results/"
    with open(path + acceptfile,"r") as h:
        x = h.readlines()
        x = [int(i.strip()) for i in x]
       # print(x)
        print(len(x))
       
        print(x.count(1)/len(x))
        

def analyse_run(filename,ylim_a=0,ylim_b=None):
    read_trace( filename +"_trace.txt",ylim_a,ylim_b)
    get_acceptratio(filename + "_acceptratio.txt")



def read_true_mod(filename):
    d = {"value":[],"iteration":[]}
    print(d)
    with open(f'../Results_plots/{filename}', "r")  as h: 
        for i,line in enumerate(h):
            l= line
            l = l.strip()
            if len(l)>1:
                value, iteration = l.split(",")[0],l.split(",")[1]
                d["value"].append(value)
                d["iteration"].append(iteration)
    print(d)
    
    print(len(d["value"]))
    print(len(d["iteration"]))
            
            



read_true_mod("true_mod_200.txt")   

#analyse_run("sim214.200_1",-3600,-3500)
#analyse_run("sim214.500_1",-8870,-7000)
#analyse_run("sim214.800_1",-14200,-14100)

#analyse_run("sim215.200_1")
#analyse_run("sim215.500_1")
#analyse_run("sim215.800_1")

#analyse_run("sim212.200_1",-3800,-3750)
#analyse_run("sim212.500_1",-8750,-8700)
#analyse_run("sim212.800_1",-14540,-14480)

#analyse_run("sim212.200_2",-3800,-3750)
#analyse_run("sim212.500_2",-8750,-8700)
#analyse_run("sim212.800_2",-14540,-14480)
#
#analyse_run("sim212.200_3",-3800,-3750)
#analyse_run("sim212.500_3",-8750,-8700)
#analyse_run("sim212.800_3",-14540,-14480)

#file1 ="sim212.800_1"
#file2 = "sim212.800_2"
#file3 = "sim212.800_3"
#file4 = "sim212.800_4"

#analyse_run(file1,-14540,-14480) OBS DOES NOT WORK IF THIS IS NOT UNCOMMENTED
 
#analyse_run(file2,-14540,-14480)

#analyse_run(file3,-14540,-14480)

#analyse_run(file4,-14540,-14480)





if __name__ == "__main__":
    #analyse_run( "test0_cprop1.2000",-6550,-6500)
   # analyse_run( "test0_cclade.2000",-6550,-6500)
    #analyse_run( "test0_c_hybrid.2000",-6550,-6500)
    #analyse_run( "test0_c_prop2.2000",-6550,-6500)
   # analyse_run("sim0.200" + "prop1" +"_oliv",-3920,-3850)
   #analyse_run("sim0.200" + "prop2" +"_oliv",-3920,-3850)
   # analyse_run("sim0.200" + "hybrid" +"_oliv",-3920,-3850)
    #analyse_run("sim0.200" + "clade" +"_oliv",-3920,-3850)
    print("")
   # analyse_run("sim0.500" + "prop1" +"_oliv",-9600,-9500)
   # analyse_run("sim0.500" + "prop2" +"_oliv",-9600,-9500)
  #  analyse_run("sim0.500" + "hybrid" +"_oliv",-9600,-9500)
  #  analyse_run("sim0.500" + "clade" +"_oliv",-9600,-9500)
    print("")
   # analyse_run("sim0.800" + "prop1" +"_oliv",-15600,-15500)
  #  analyse_run("sim0.800" + "prop2" +"_oliv",-15600,-15500)
  #  analyse_run("sim0.800" + "hybrid" +"_oliv",-15600,-15500)
  #  analyse_run("sim0.800" + "clade" +"_oliv",-15600,-15500)
    
                
 
   # read_trace(path + "test0_b_trace.txt",-1310,-1299)
   # read_trace("test0_b_2_trace.txt",-1325,-1310)
   # read_trace("test0_b_3_trace.txt",-1310,-1299)
   # read_trace("test0_b_4_trace.txt",-1325,-1310)
    #read_trace("test0_b_3_trace.txt",-1314,-1298)
   #read_trace("test0_b_prop2_an_2_trace.txt",-1310,-1297)
  