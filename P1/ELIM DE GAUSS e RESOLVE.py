# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 19:24:27 2020

@author: murat
"""
import numpy as np
import math 
from decimal import *

s=int(input("Qts algarismos significativos?: "))
getcontext().prec = s
d=int(input("Qual a dimensao do sistema?: "))
A=[]
b = []
for i in range(d):
    A.append([])
    for j in range (d):
        uou=1*Decimal(input("Insira o elemento %d da linha %d: "%(j+1,i+1)))
        A[i].append(uou)
    iei=1*Decimal(input("Insira o elemento %d do lado direito: "%(i+1)))
    b.append(iei)


def linearsolver(A,b):
  n = len(A)
  M = A
  L=np.diag(np.zeros(d))
  i = 0
  for x in M:
   x.append(b[i])
   i += 1

  for k in range(n):
   for i in range(k,n):
     if abs(M[i][k]) > abs(M[k][k]):
        M[k], M[i] = M[i],M[k]
        eae=L[k].copy()
        L[k]=L[i]
        L[i]=eae
     else:
        pass

   for j in range(k+1,n):
       q = ((M[j][k]) / M[k][k])
       L[j][k]=Decimal(q)
       for m in range(k, n+1):
          M[j][m] -= Decimal(q) *Decimal(M[k][m])

  x = [0 for i in range(n)]

  x[n-1] =Decimal((M[n-1][n])/M[n-1][n-1])
  for i in range (n-1,-1,-1):
    z = 0
    for j in range(i+1,n):
        z = Decimal(z  + Decimal((M[i][j])*x[j]))
    x[i] =Decimal( Decimal((M[i][n] - z))/M[i][i])
  
  L=L+np.diag(np.ones(d))  
  print ("Solução:")
  print(np.array(x).astype(float))
  print("")
  print("matriz L:")
  print(L)
  print("")
  print("matriz U:")
  print(np.array(M)[:,0:3].astype(float))

(linearsolver(A,b))