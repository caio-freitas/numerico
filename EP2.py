import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
import random

#############################################################################################################
########################################## Funcao da fonte ##################################################
#############################################################################################################

def funcao_fonte(x_pt, t, ponto):

    global item

    if (x_pt >= (ponto - h/2) and x_pt <= (ponto + h/2)):

        r_t = 10*(1 + np.cos(5*t))
        return r_t/delta_x

    else:

        return 0


#############################################################################################################
####################################### Condicao inicial (t=0) ##############################################
#############################################################################################################

def cond_ini (x):

    return 0

#############################################################################################################
######################################## Condicoes de contorno ##############################################
#############################################################################################################

def g1(t): # condição de contorno, x=0


    return 0


def g2(t): # condição de contorno, x=1


    return 0
    

'''
Função que gera resultado da decomposição LDLt para a matriz
tridiagonal presente no método de Crank-Nicolson 
'''
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

'''
Resolução de um sistema Ax = b, onde A pode ser:
- "lower": triangular inferior
- "upper": triangular superior
- "diagonal": diagonal
'''
def solve(A, b, m_type):
    if A.shape[0] != A.shape[1]:
        print("Erro, matriz nao quadrada!")
        return -1
    solution = np.zeros((A.shape[0]))
    if m_type == "lower":
        solution[0] = b[0]/A[0][0]
        
        for i in range(1, A.shape[0]):
            summ = 0
            for j in range(1, i + 1):
                summ += A[i][i-j]*solution[i-j]
            solution[i] = (b[i] - summ)
            
        return solution

    elif m_type == "upper":
        solution[A.shape[0] - 1] = b[A.shape[0]-1]/A[A.shape[0] - 1][A.shape[0] - 1]
        for i in range(2, A.shape[0]+1):
            j = A.shape[0] - i
            summ = 0
            for k in range(1, i):
                summ += A[j][j+k]*solution[j+k]
            solution[j] = (b[j] - summ)
            
        return solution

    elif m_type == "diagonal":
        for i in range(A.shape[0]):
            solution[i] = b[i]/A[i][i]
        return solution

'''
Decomposição LDLt de uma matriz A genérica
'''

def decomposeLDLt(A):
    L = np.diag(np.ones(len(A)))
    D = np.diag(np.ones(len(A)))
    v = np.ones(len(A))
    for i in range(len(A)):
        for j in range(i):
            v[j] = L[i][j]*D[j][j]
        summ = 0
        for k in range(i):
            summ += L[i][k]*v[k]
        D[i][i] = A[i][i] - summ
        
        
        for j in range(i+1, len(A)):
            summ = 0
            for k in range(i):
                summ += L[j][k]*v[k]
            L[j][i] = (A[j][i] - summ)/D[i][i]
    return L, D


#############################################################################################################
######################################### Resolucao Crank-Nicolson ##########################################
#############################################################################################################

#------------------------------ Funcao que resolve por Crank-Nicholson --------------------------------------

def Crank (pontos): # Retorna a matriz com os UTk nas colunas

    pixels_tot = (N-1)*(N-1)*len(pontos) + len(pontos)*(N-1)
    porc = int(pixels_tot/100)
    cont_porc = 0
    cont_pixels = 0
    Temps_T = np.zeros((N-1, len(pontos)))
    contou = False

    for fonte in range (len(pontos)):

        #----------------------------------------- Matrizes iniciais ------------------------------------------------

        matriz_final = np.zeros((N+1,N+1))

        # Problema do tipo Ax = b...

        x = np.zeros ((N-1, 1))
        b = np.zeros ((N-1))
        x_linha = np.zeros ((N-1))

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
                    b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k, pontos[fonte])+funcao_fonte(delta_x*i, delta_t*(k+1), pontos[fonte]))+(lamb/2)*(g1(delta_t*(k+1))+g1(delta_t*k)+matriz_final[2][k]) + (1-lamb)*matriz_final[1][k]
                
                elif (i==N-1):
                    b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k, pontos[fonte])+funcao_fonte(delta_x*i, delta_t*(k+1), pontos[fonte]))+(lamb/2)*(g2(delta_t*(k+1))+g2(delta_t*k)+matriz_final[N-2][k]) + (1-lamb)*matriz_final[i][k]

                else:
                    b[i-1] = (delta_t/2)*(funcao_fonte(delta_x*i, delta_t*k, pontos[fonte])+funcao_fonte(delta_x*i, delta_t*(k+1), pontos[fonte]))+(lamb/2)*(matriz_final[i-1][k] + matriz_final[i+1][k]) + (1-lamb)*matriz_final[i][k]

                cont_pixels += 1
                
                if (cont_pixels%porc == 0):
                    if (contou == False):
                            C = (datetime.datetime.now() - A)*100
                            print("Estimativa de tempo = {}\n".format(C))
                            contou = True
                    cont_porc += 1
                    print(cont_porc, "% concluido")

            y = solve(L_matriz, b, "lower")
            z = solve(D_matriz, y, "diagonal")
            x_linha = solve(LT_matriz, z, "upper")

            for j in range (0, N-1):
                x[j][0] = x_linha[j]

            for i in range (1, N):
                matriz_final[i][k+1] = x[i-1][0]

        #---------------------------------------- Preenchendo com uT(xi) --------------------------------------------

        for i in range (1, N):
            Temps_T[i-1][fonte] = matriz_final[i][N]

    return Temps_T


#############################################################################################################
####################################### Montando o problema #################################################
#############################################################################################################


def Montar (Temps_T, UTs): # Retorna N e b do problema Nx = b, se Uts é um vetor coluna

    Normal = np.transpose(Temps_T).dot(Temps_T)
    b = np.transpose(Temps_T).dot(UTs)

    return (Normal, b)

#############################################################################################################
###################################### Leitura do teste.txt #################################################
#############################################################################################################

def Ler (nome): # 'points' contem os valores dos pontos de aplicacao da fonte e 'content' os Uts

    file = open(nome, 'r+')
    content = []
    points = []
    cont = 0
    for line in file:
        if (cont == 0):
            points = line.split()
        else:
            content.append(line.strip())
        cont+=1

    file.close()

    for i in range (len(content)):
        content[i] = float(content[i])

    for i in range (len(points)):
        points[i] = float(points[i])
    
    return (content, points)

#############################################################################################################
######################################### Erro quadratico ###################################################
#############################################################################################################

def Erro_2 (aks, UTs, Temps_T, nf):

    soma_N = 0
    soma_K = 0

    for i in range (0, N-1):
        
        for k in range (nf):

            soma_K = soma_K + aks[k]*Temps_T[i][k]

        soma_N = soma_N + (UTs[i] - soma_K)**2
        soma_K = 0
    Erro = np.sqrt(delta_x*soma_N)

    return (Erro)

#############################################################################################################
########################################### Cria graficos ###################################################
#############################################################################################################

def Cria_graficos (aks, UTs, Temps_T, nf, Erro_2):

    soma_K = 0
    label_x = np.zeros((N+1))
    resultado = np.zeros((N+1))
    Erro = np.zeros((N+1))
    Erro_2_graf = np.full((N+1), Erro_2)
    
    for i in range (N+1):
        label_x[i] = delta_x*i

    for i in range (1, N):
        
        for k in range (nf):

            soma_K = soma_K + aks[k]*Temps_T[i-1][k]
        
        resultado[i] = soma_K
        soma_K = 0

    plt.plot(label_x, resultado)
    plt.title("Resultado UT(x), N=" + str(N) + ", Item " + item)
    plt.xlabel("Posição")
    plt.ylabel("Temperatura")
    plt.grid(True)
    plt.savefig("Result_item-{}_N-{}.png".format(item, str(N)), dpi=400)
    #plt.show()
    plt.close()

    for i in range (1, N):

        Erro[i] = abs(resultado[i]-UTs[i-1])

    plt.plot(label_x, Erro, '-b', label = 'Módulo do erro instântaneo')
    plt.plot(label_x, Erro_2_graf, '-r', label = 'Erro quadrático médio')
    plt.title("Erro, N="+str(N) + ", Item " + item)
    plt.xlabel("Posição")
    plt.ylabel("Temperatura")
    plt.grid(True)
    plt.legend()
    plt.savefig("Erro_item-{}_N-{}.png".format(item, str(N)), dpi=400)
    #plt.show()
    plt.close()



#############################################################################################################
########################################### Itens pedidos ###################################################
#############################################################################################################

item = input("Digite o item a ser resolvido (a,b,c,d): ")

if (item == 'a'):

    #---------------------------------------- Parametros iniciais -----------------------------------------------

    N = 128
    M = N
    lamb = N
    delta_t = 1/N
    delta_x = 1/N
    L_matriz, D_matriz, LT_matriz = decomporLDL(lamb/2)
    h = delta_x
    pontos = np.array([0.35])
    UTs = np.zeros((N-1, 1))
    A = datetime.datetime.now()

    #---------------------------------------- Inicio da resolucao -----------------------------------------------

    Temps_T = Crank(pontos)

    for i in range (1, N):
        UTs[i-1][0] = 7*Temps_T[i-1][0]

    Normal, b = Montar(Temps_T, UTs)
    NL, ND = decomposeLDLt(Normal)
    y = solve(NL, b, "lower")
    z = solve(ND, y, "diagonal")
    aks = solve(NL.T, z, "upper")
    Erro = Erro_2(aks, UTs, Temps_T, len(pontos))

    print("a1 é igual a: ", aks[0])
    print("O erro quadratico eh: ", Erro)

    Cria_graficos(aks, UTs, Temps_T, len(pontos), Erro)
    print("O processo demorou ", datetime.datetime.now()-A)
    

elif (item == 'b'):

    #---------------------------------------- Parametros iniciais -----------------------------------------------

    N = 128
    M = N
    lamb = N
    delta_t = 1/N
    delta_x = 1/N
    L_matriz, D_matriz, LT_matriz = decomporLDL(lamb/2)
    h = delta_x
    pontos = np.array([0.15, 0.3, 0.7, 0.8])
    UTs = np.zeros((N-1, 1))
    A = datetime.datetime.now()

    #---------------------------------------- Inicio da resolucao -----------------------------------------------

    Temps_T = Crank(pontos)

    for i in range (1, N):
        UTs[i-1][0] = 2.3*Temps_T[i-1][0] + 3.7*Temps_T[i-1][1] + 0.3*Temps_T[i-1][2] + 4.2*Temps_T[i-1][3]

    Normal, B = Montar(Temps_T, UTs)


    NL, ND = decomposeLDLt(Normal)
    y = solve(NL, B, "lower")
    z = solve(ND, y, "diagonal")
    aks = solve(NL.T, z, "upper")
    Erro = Erro_2(aks, UTs, Temps_T, len(pontos))

    print("Os valores 'ak' são: ")
    print(aks)
    print("Com decomposição LDLt: {}".format(aks))
    print("O erro quadratico eh: ", Erro)

    Cria_graficos(aks, UTs, Temps_T, len(pontos), Erro)
    print("O processo demorou ", datetime.datetime.now()-A)
    


elif (item == 'c'):

    #---------------------------------------- Parametros iniciais -----------------------------------------------

    N = int(input("Digite o valor de N: "))
    M = N
    lamb = N
    delta_t = 1/N
    delta_x = 1/N
    L_matriz, D_matriz, LT_matriz = decomporLDL(lamb/2)
    h = delta_x
    todos_UTs, pontos = Ler("teste.txt")
    UTs = []
    A = datetime.datetime.now()

    #---------------------------------------- Inicio da resolucao -----------------------------------------------

    fator = 2048/N
    cont = 1

    for i in range (1, 2048):
        if (i%fator == 0):
            UTs.append(todos_UTs[i])
        cont+=1

    Temps_T = Crank(pontos)
    Normal, B = Montar(Temps_T, UTs)

    NL, ND = decomposeLDLt(Normal)
    y = solve(NL, B, "lower")
    z = solve(ND, y, "diagonal")
    aks = solve(NL.T, z, "upper")
    #aks = np.linalg.inv(Normal).dot(B)
    Erro = Erro_2(aks, UTs, Temps_T, len(pontos))

    print("Os valores 'ak' são: ")
    print(aks)
    print("O erro quadratico eh: ", Erro)

    Cria_graficos(aks, UTs, Temps_T, len(pontos), Erro)
    print("O processo demorou ", datetime.datetime.now()-A)

elif (item == 'd'):

    #---------------------------------------- Parametros iniciais -----------------------------------------------

    N = int(input("Digite o valor de N: "))
    M = N
    lamb = N
    delta_t = 1/N
    delta_x = 1/N
    L_matriz, D_matriz, LT_matriz = decomporLDL(lamb/2)
    h = delta_x
    todos_UTs, pontos = Ler("teste.txt")
    UTs = []
    A = datetime.datetime.now()

    #---------------------------------------- Inicio da resolucao -----------------------------------------------

    fator = 2048/N
    cont = 1

    for i in range (1, 2048):
        if (i%fator == 0):
            UTs.append(todos_UTs[i]*(1 + 2*(random.random()-0.5)*0.01))
        cont+=1

    Temps_T = Crank(pontos)
    Normal, B = Montar(Temps_T, UTs)

    NL, ND = decomposeLDLt(Normal)
    y = solve(NL, B, "lower")
    z = solve(ND, y, "diagonal")
    aks = solve(NL.T, z, "upper")
    Erro = Erro_2(aks, UTs, Temps_T, len(pontos))

    print("Os valores 'ak' são: ")
    print(aks)
    print("O erro quadratico eh: ", Erro)

    Cria_graficos(aks, UTs, Temps_T, len(pontos), Erro)
    print("O processo demorou ", datetime.datetime.now()-A)


else:
    print("Voce digitou um item não válido")