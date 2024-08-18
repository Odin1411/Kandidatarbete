# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:27:10 2024

@author: olivi
"""
import os

print(os.getcwd())
os.chdir('../')
print(os.getcwd())

print( os.path.normpath(os.getcwd() + os.sep + os.pardir))

print(os.getcwd())
