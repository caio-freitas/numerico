import numpy as np

A = np.array([[2, -1, 8],
              [-1, 1/2, 1],
              [0, 2, 1]])

b = np.array([[4],
             [1/2],
             [3/2]])

sig = 3
troca = False
linha = np.zeros((1,len(A)))
linha_tr = np.zeros((1,len(A)))

P = np.zeros((len(A), len(A)))
np.fill_diagonal(P, 1)
#print(P)

#--------------------------------------------------------------------------------------------------
#-----------------------------------Condensação pivotal--------------------------------------------
#--------------------------------------------------------------------------------------------------
for i in range (len(A)):

    atual = A[i][i]

    for j in range (len(A)-i):

        if(A[len(A)-j-1][i] > atual):
            atual = A[len(A)-j-1][i]
            num_linha = len(A)-j-1
            for k in range (len(A)):
                linha[0][k] = A[len(A)-j-1][k]
            troca = True

    if (troca == True):

        for k in range (len(A)):
            linha_tr[0][k] = A[i][k]
        for k in range (len(A)):
            A[num_linha][k] = linha_tr[0][k]
            A[i][k] = linha[0][k]
        P[[i,num_linha]] = P[[num_linha, i]]
    troca = False

print("###############################################################################")
print("################ Condensação pivotal - Mudança de linhas ######################")
print("###############################################################################")
print("Matriz A permutada:\n")
print(A, "\n")
print("Matriz P:\n")
print(P, "\n")
b = P.dot(b)
#--------------------------------------------------------------------------------------------------
#------------------------------Arredondamento pelo significativo-----------------------------------
#--------------------------------------------------------------------------------------------------

def Arredondar (num, sig):
    mult = 0
    arredondado = False
    num_mod = abs(num)
    if (num == 0):
        return 0
    if(num < 1 and num > -1):
        while (arredondado == False):
            num_mod = num_mod*10
            if (int(num_mod) != 0):
                arredondado = True
            else:
                mult+=1
        num = np.round(num, sig+mult)
    else:
        while (arredondado == False):
            
            if (int(num_mod) == 0):
                arredondado = True
            else:
                mult+=1
            num_mod = num_mod/10
        if (sig-mult < 0 ):
            num = np.round(num, 0)
        else:
            num = np.round(num, sig-mult)
    return num

#--------------------------------------------------------------------------------------------------
#-----------------------------------Resolvendo o sistema-------------------------------------------
#--------------------------------------------------------------------------------------------------

LU = np.zeros((len(A), len(A)))
mudou_b = False
for i in range (len(A)-1):

    for j in range (len(A)-i-1):

        fator = Arredondar(A[len(A)-j-1][i]/A[i][i], sig)      

        for k in range (i, len(A)):

            if (mudou_b == False):

                temp_b = -1*(Arredondar(b[i][0]*fator, sig)) + b[len(A)-j-1][0]
                b[len(A)-j-1][0] = Arredondar(temp_b, sig)
                mudou_b = True

            temp = -1*(Arredondar(A[i][k]*fator, sig)) + A[len(A)-j-1][k]
            A[len(A)-j-1][k] = Arredondar(temp, sig)

        A[len(A)-j-1][i] = fator
        mudou_b = False

print("###############################################################################")
print("############################# Sist resolvido ##################################")
print("###############################################################################")
print("Matriz A resolvida:\n")
print(A, "\n")
print("Matriz b:\n")
print(b, "\n")

#--------------------------------------------------------------------------------------------------
#-----------------------------------Calculo de L e U-----------------------------------------------
#--------------------------------------------------------------------------------------------------

L = np.zeros((len(A), len(A)))
np.fill_diagonal(L, 1)

U = np.zeros((len(A), len(A)))

for i in range (len(A)-1):
    for k in range (i + 1, len(A)):
        L[k][i] = A[k][i]

for i in range (len(A)):
    for k in range (1+i):
        U[k][i] = A[k][i]

print("###############################################################################")
print("############################# Matrizes L e U ##################################")
print("###############################################################################")
print("Matriz L:\n")
print(L, "\n")
print("Matriz U:\n")
print(U, "\n")

#--------------------------------------------------------------------------------------------------
#-----------------------------------Calculo dos X's------------------------------------------------
#--------------------------------------------------------------------------------------------------

x = np.zeros((len(A), 1))

for i in range (len(A)-1, -1, -1):

    soma = 0

    for k in range (len(A)):
        if(k != i):
            soma = soma + Arredondar(A[i][k]*x[k][0], sig)
    x[i][0] = Arredondar(((b[i][0]-soma)/A[i][i]), sig)

print("###############################################################################")
print("############################# Solução x #######################################")
print("###############################################################################")
print("X:\n")
print(x, "\n")