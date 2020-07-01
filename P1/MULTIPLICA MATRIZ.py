# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 22:10:40 2020

@author: murat
"""
import numpy as np
from decimal import *
s=int(input("Qts algarismos significativos?: "))
getcontext().prec = s

m_A=int(input("Quantas linhas tem A?: "))
n_A=int(input("Quantas colunas tem A?: "))
A=[]
for i in range(m_A):
    A.append([])
    for j in range (n_A):
        uou=1*Decimal(input("Insira o elemento %d da linha %d: "%(j+1,i+1)))
        A[i].append(uou)
        
m_B=int(input("Quantas linhas tem B?: "))
n_B=int(input("Quantas colunas tem B?: "))
B=[]
for i in range(m_B):
    B.append([])
    for j in range (n_B):
        uou=1*Decimal(input("Insira o elemento %d da linha %d: "%(j+1,i+1)))
        B[i].append(uou)

result = []
for i in range(m_A):
    result.append([])
    for j in range (n_B):
        result[i].append(0)

# iterate through rows of X
for i in range(len(A)):
   # iterate through columns of Y
   for j in range(len(B[0])):
       # iterate through rows of Y
       for k in range(len(B)):
           result[i][j] += Decimal(A[i][k]) * Decimal(B[k][j])

print(np.array(result).astype(float))
