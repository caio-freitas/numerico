# -*- coding: utf-8 -*-
 
'''
@author: gabrieldacunha

'''

def gaussSeidel(A, b, x, N, tol, maxIteracoes):
    xAnt = [0.0 for i in range(N)]
    for i in range(maxIteracoes):
        for j in range(N):
            xAnt[j] = x[j]
        for j in range(N):
            soma = 0.0
            for k in range(N):
                if (k != j):
                    soma = soma + A[j][k] * x[k]
            x[j] = (b[j] - soma) / A[j][j]
        difNorma = 0.0
        normaAnt = 0.0
        for j in range(N):
            difNorma = difNorma + abs(x[j] - xAnt[j])
            normaAnt = normaAnt + abs(xAnt[j])  
        if normaAnt == 0.0:
            normaAnt = 1.0
        norma = difNorma / normaAnt
        if (norma < tol) and i != 0:
            print("Usando Gauss-Seidel, a solução converge para [", end="")
            for j in range(N - 1):
                print(x[j], ",", end="")
            print(x[N - 1], "]. Levou", i + 1, "iterações.")
            return
    print("A solução não convergiu antes de atingir o limite de iterações.")
    print("Usando Gauss-Seidel, até então a solução estava em: [", end="")
    for j in range(N - 1):
        print(x[j], ",", end="")
    print(x[N - 1], "].")

def jacobi(A, b, x, N, tol, maxIteracoes):

    xAnt = [0.0 for i in range(N)]
    for i in range(maxIteracoes):
        for j in range(N):
            xAnt[j] = x[j]
        for j in range(N):
            soma = 0.0
            for k in range(N):
                if (k != j):
                    soma = soma + A[j][k] * xAnt[k]
            x[j] = (b[j] - soma) / A[j][j]
        difNorma = 0.0
        normaAnt = 0.0
        for j in range(N):
            difNorma = difNorma + abs(x[j] - xAnt[j])
            normaAnt = normaAnt + abs(xAnt[j])  
        if normaAnt == 0.0:
            normaAnt = 1.0
        norma = difNorma / normaAnt
        if (norma < tol) and i != 0:
            print("Usando Jacobi, a solução converge para [", end="")
            for j in range(N - 1):
                print(x[j], ",", end="")
            print(x[N - 1], "]. Levou", i + 1, "iterações.")
            return
    print("A solução não convergiu antes de atingir o limite de iterações.")
    print("Usando Jacobi, até então a solução estava em: [", end="")
    for j in range(N - 1):
        print(x[j], ",", end="")
    print(x[N - 1], "].")

d=int(input("Qual é a dimensao do sistema?: "))
A=[]
b = []
x0 = []
for i in range(d):
    A.append([])
    for j in range (d):
        linha=float(input("Insira o elemento %d da linha %d da matriz A: "%(j+1,i+1)))
        A[i].append(linha)

for i in range(d):
    elementoB=float(input("Insira o elemento %d do vetor b: "%(i+1)))
    b.append(elementoB)

for i in range(d):
    elementoX=float(input("Insira o elemento %d da solução inicial: "%(i+1)))
    x0.append(elementoX)

repetir = "S"
while(repetir == "S" or repetir == "s"):
    print("1 - Jacobi")
    print("2 - Gauss-Seidel")
    metodo = input("Escolha o método: ")

    matrizA = []
    vetorB = []
    vetorX = []
    for i in range(d):
        matrizA.append([])
        for j in range (d):
            matrizA[i].append(A[i][j])
        vetorB.append(b[i])
        vetorX.append(x0[i])

    maxIteracoes = int(input("Digite o número máximo de iterações: "))
    tol=float(input("Qual é a precisão?: "))

    if(metodo == "1"):
        jacobi(matrizA, vetorB, vetorX, d, tol, maxIteracoes)

    if (metodo == "2"):
        gaussSeidel(matrizA, vetorB, vetorX, d, tol, maxIteracoes)

    repetir = input("Deseja testar outro método/precisão? [S/N]: ")

