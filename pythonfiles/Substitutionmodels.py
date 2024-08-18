# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:05:19 2024

@author: olivi
"""
import numpy
from modelmatcher import RateMatrix
import numpy as np
import random


np.set_printoptions(precision=3, suppress=True)

model = RateMatrix.instantiate('WAG')
Pt = model.get_replacement_probs(0.1)
#print(Pt)

def calcP(Q,t,n=10000):
    
    #obs this is not a great way of calculating P 
    '''
    Transition matrix for a branch of length t
    '''
    I = numpy.identity(len(Q)) 
    
    M0 = (numpy.add(I,(Q*t)/n )) 
    M = (numpy.add(I,(Q*t)/n )) 
    
    for i in range(n-1):
        
        M = numpy.dot(M0,M)
    M= numpy.round(M, decimals=3)

    return M

class Substitutionmodel:
    """class for substitutionmodels"""
    def __init__(self,seqtype):
        "supports nucleotide and protein models"
        assert seqtype == "nucleotide" or seqtype == "protein"
        self.seqtype = seqtype
        self.nucleobases = {"A":0,"C":1,"G":2,"T":3}
        aminochars = "ARNDCQEGHILKMFPSTWYV"
        Aminostates = {}
        for i,a in enumerate(aminochars):
            Aminostates[a] = i
        self.aminoacids = Aminostates
    
    def JC_Q(self):
        a = numpy.round(1/3,6)

        JC = numpy.array( [ [ -1.0000, a ,a, a] ,
             [a, -1.0000, a, a ] ,
             [a,a, -1.0000, a ] ,
             [a,a,a, -1.0000 ] ]) 
        return JC
    
    def get_stationary(self,name): 
        if self.seqtype == "nucleotide":
            return (self.calculateP(name,1000))[0]
        elif self.seqtype == "protein":
            model = RateMatrix.instantiate(name)
            return model.freq
    
    def calculateP(self,name,t): #only works for aminoacid model
        if self.seqtype == "protein":
            model = RateMatrix.instantiate(name)
            Pt = model.get_replacement_probs(t)
            return Pt
       
        elif self.seqtype == "nucleotide":
            if name == "JC":
                return calcP(self.JC_Q(), t)


 
subm = Substitutionmodel("nucleotide")  
#print(subm.get_stationary("JC")  )

#Subm = Substitutionmodel("protein")
#model = RateMatrix.instantiate("WAG")
#a  = model.freq #model.freq to get frequency
#print(a)

#np.set_printoptions(precision=3, suppress=True)

#models = RateMatrix.get_all_models()
#for m in models:
#   print(m.get_name())

#model = RateMatrix.instantiate('WAG')
#Pt = model.get_replacement_probs(0.1)
#print(Pt)

def JukesCantorQ():
    a = numpy.round(1/3,6)

    JC = numpy.array( [ [ -1.0000, a ,a, a] ,
         [a, -1.0000, a, a ] ,
         [a,a, -1.0000, a ] ,
         [a,a,a, -1.0000 ] ]) 
    return JC
    



def TransitionProbability(Q,i,j,t):
    P = calcP(Q,t) 
    return P[i][j]
