import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("Agg")

#############################################################################################################
########################################## Funcao da fonte ##################################################
#############################################################################################################

def funcao_fonte(x_pt, t): # funcao de entrada
    global item
    #return 10*(x_pt**2)*(x_pt-1) - 60*x_pt*t + 20*t
    if item == 'a':
        return 10*(np.cos(10*t))*(x_pt**2)*((1-x_pt)**2) - (1+(np.sin(10*t)))*(12*(x_pt**2)-(12*x_pt)+2) # item a
    elif item == 'b':
        return (np.exp(t-x_pt)*(-np.sin(5*t*x_pt)*5*x_pt - 10*t*np.sin(5*t*x_pt)+np.cos(5*t*x_pt)*25*t*t))   # item b
    elif item == 'c':
        if x_pt == int(N*0.25):
            return (10000*(1 - 2*t*t))/delta_x
        else:
            return 0

#############################################################################################################
####################################### Condicao inicial (t=0) ##############################################
#############################################################################################################

def cond_ini (x):
    #return ((x**2)*((1-x)**2))
    if item =='a':
        return ((x**2)*((1-x)**2))
    elif item =='b':
        return np.exp(-x)
    elif item == 'c':
        return funcao_fonte(x, 0)

#############################################################################################################
######################################## Condicoes de contorno ##############################################
#############################################################################################################

def g1(t): # condição de contorno, x=0

    if item == 'a' or item == 'c':
        return 0
    elif item =='b':
        return (np.exp(t))

def g2(t): # condição de contorno, x=1

    if item =='a' or item == 'c':

        return 0
    elif item =='b':

        return (np.exp(t-1)*np.cos(5*t))

#############################################################################################################
########################################## Decomposicao LDLt ################################################
#############################################################################################################

def decomporLDL(ł):
    # L e D Vetores
    D = np.zeros(N-1)
    L = np.zeros(N-2)

    # Condição inicial para prossegir com os cálculos
    D[0] = 1 + 2*ł
    # Cálculo dos valores dos vetores L e D
    for i in range (1, len(D)):
        D[i] = 1+2*ł - (((-ł/D[i-1])**2)*D[i-1])

    for i in range (0, len(L)):
        L[i] = -ł/D[i]

    # Criação das matrizes L e D a partir dos vetores
    D_matriz = np.zeros((N-1, N-1))
    L_matriz = np.zeros((N-1, N-1))

    for i in range (0, N-1):
        D_matriz[i][i] = D[i]
        L_matriz[i][i] = 1      # A matriz L deve ter 1's na diagonal principal

    for i in range (0, N-2):
        L_matriz[i+1][i] = L[i] #Colocando a subdiagonal de L

    LT_matriz = np.transpose(L_matriz)

    return (L_matriz, D_matriz, LT_matriz)

#############################################################################################################
####################################### Resolucao LDLt*x = b ################################################
#############################################################################################################

def solve(A, b, m_type):
    if A.shape[0] != A.shape[1]:
        print("Erro, matriz nao quadrada!")
        return -1
    solution = np.zeros((A.shape[0]))
    if m_type == "lower":
        solution[0] = b[0]/A[0][0]
        for i in range(1, A.shape[0]):
            solution[i] = b[i] - A[i][i-1]*solution[i-1]
        return solution

    elif m_type == "upper":
        solution[A.shape[0] - 1] = b[A.shape[0]-1]/A[A.shape[0] - 1][A.shape[0] - 1]
        for i in range(2, A.shape[0]-1):
            j = A.shape[0] - i
            solution[j] = b[j] - A[j][j+1]*solution[j+1]
        return solution

    elif m_type == "diagonal":
        for i in range(A.shape[0]):
            solution[i] = b[i]/A[i][i]
        return solution


#############################################################################################################
########################################### Resolucao Crank #################################################
#############################################################################################################

#---------------------------------------- Parametros iniciais -----------------------------------------------

N = 640

M = N
lamb = N
item = 'a'
delta_t = 1/N
delta_x = 1/N

#----------------------------------------- Matrizes iniciais ------------------------------------------------

matriz_final = np.zeros((N+1,N+1))

# Problema do tipo Ax = b...

#A = np.zeros ((N-1, N-1))
x = np.zeros ((N-1, 1))
b = np.zeros ((N-1))
x_linha = np.zeros ((N-1))

L_matriz, D_matriz, LT_matriz = decomporLDL(lamb/2)

#----------------------------------------- Condicoes iniciais -----------------------------------------------

for i in range (N+1):
    matriz_final[i][0] = cond_ini(i*delta_x)

for k in range (1, N+1):
    matriz_final[0][k] = g1(k*delta_t)

for k in range (1, N+1):
    matriz_final[N][k] = g2(k*delta_t)

#----------------------------------------- Inicio do calculo ------------------------------------------------

for k in range (0, N):

    for i in range (1, N):
        
        if (i==1):
             b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k)+funcao_fonte(delta_x*i, delta_t*(k+1)))+(lamb/2)*(g1(delta_t*(k+1))+g1(delta_t*k)+matriz_final[2][k]) + (1-lamb)*matriz_final[1][k]
        
        elif (i==N-1):
             b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k)+funcao_fonte(delta_x*i, delta_t*(k+1)))+(lamb/2)*(g2(delta_t*(k+1))+g2(delta_t*k)+matriz_final[N-2][k]) + (1-lamb)*matriz_final[i][k]

        else:
            b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k)+funcao_fonte(delta_x*i, delta_t*(k+1)))+(lamb/2)*(matriz_final[i-1][k] + matriz_final[i+1][k]) + (1-lamb)*matriz_final[i][k]

    # y = solve(L_matriz, b, "lower")
    # z = solve(D_matriz, y, "diagonal")
    # x_linha = solve(LT_matriz, z, "upper")

    x_linha = np.linalg.inv(L_matriz.dot(D_matriz).dot(LT_matriz)).dot(b)

    for j in range (0, N-1):
        x[j][0] = x_linha[j]

    for i in range (1, N):
        matriz_final[i][k+1] = x[i-1][0]
    

#############################################################################################################
##################################### Plotando matriz final #################################################
#############################################################################################################

plt.figure()
#plt.xlabel("Tempo (t)")
#plt.ylabel("Posicao (x)")
plt.matshow(matriz_final,cmap = 'hot')
plt.savefig("Resultado-Crank-AA.png", dpi=400)
plt.show()
plt.close()