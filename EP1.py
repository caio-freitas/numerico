import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import datetime
import time
from numba import jit

t_total = 1     # Tempo total de "simulacao"
x_total = 1     # Tamanho da "barra"

show=False

item = 'a'

@jit
def funcao_fonte(x_pt, t, N): # funcao de entrada
    global item
    #return 10*(x_pt**2)*(x_pt-1) - 60*x_pt*t + 20*t
    if item == 'a':
        return 10*(np.cos(10*t))*(x_pt**2)*((1-x_pt)**2) - (1+(np.sin(10*t)))*(12*(x_pt**2)-12*x_pt+2) # item a
    elif item == 'b':
        return (np.ex_ptp(t-x_pt)*(-np.sin(5*t*x_pt)*5*x_pt - 10*t*np.sin(5*t*x_pt)+np.cos(5*t*x_pt)*25*t*t))   # item b
    elif item == 'c':
        if x_pt == int(N*0.25):
            return (10000*(1 - 2*t*t))/delta_x
        else:
            return 0

@jit
def funcao_exata(x, t):
    #return (10*t*x*x*(x-1))
    if item =='a':
        return ((1 + np.sin(10*t))*(x**2)*((1-x)**2))
    elif item == 'b':
        return np.exp(t-x)*np.cos(5*t*x)

@jit
def gera_matriz_funcao_exata(M, N, delta_t, delta_x):
    f_exata = np.zeros((int(M)+1, int(N)+1), dtype='float32')
    for i in range(1, int(M)+1):
        t = (i-1)*delta_t
        for j in range(1, int(N)+1):
            x = j*delta_x
            f_exata[i][j] = funcao_exata(x, t)
    for i in range (0,N): # para todos os x's analisados
        f_exata[0][i] = cond_ini(i*delta_x, N)

    for j in range(0, int(M)+1): # para todos os t's analisados
        f_exata[j][0] = g1(j*delta_t)
        f_exata[j][N] = g2(j*delta_t)
    return f_exata

@jit
def cond_ini (x, N):
    #return ((x**2)*((1-x)**2))
    if item =='a':
        return ((x**2)*((1-x)**2))
    elif item =='b':
        return np.exp(-x)       
    elif item == 'c':
        return funcao_fonte(x, 0, N)

@jit
def g1(t): # condição de contorno, x=0
    
    if item == 'a' or item == 'c':
        return 0
    elif item =='b': 
        return (np.exp(t))

@jit
def g2(t): # condição de contorno, x=1
    
    if item =='a' or item == 'c':
        
        return 0
    elif item =='b':

        return (np.exp(t-1)*np.cos(5*t))

'''
A decomposição LDLt recebe

- a: Vetor com a diagonal principal da matriz
- b: Vetor com a subdiagonal da matriz

E retorna
- L: matriz bidiagonal inferior
- D: Matriz diagonal
- Lt: A transposta de L
'''
@jit
def decomporLDL(ł, N):

    # L e D Vetores
    D = np.zeros(N)
    L = np.zeros(N-1)


    # Condição inicial para prossegir com os cálculos
    D[0] = 1 + 2*ł

    # Cálculo dos valores dos vetores L e D
    for i in range (1, len(D)):
        D[i] = 1+2*ł - (((-ł/D[i-1])**2)*D[i-1])
    
    for i in range (0, len(L)):
        L[i] = -ł/D[i]

    # Criação das matrizes L e D a partir dos vetores
    D_matriz = np.zeros((N, N))
    L_matriz = np.zeros((N, N))

    for i in range (0, N):
        D_matriz[i][i] = D[i]
        L_matriz[i][i] = 1      # A matriz L deve ter 1's na diagonal principal
    
    for i in range (0, N-1):
        L_matriz[i+1][i] = L[i] #Colocando a subdiagonal de L

    LT_matriz = np.transpose(L_matriz)

    return (L_matriz, D_matriz, LT_matriz)


def invert_diagonal(M):
    inverse = np.diag(np.ones((M.shape[0])))
    for i in range(M.shape[0]):
        try:
            inverse[i][i] = 1/M[i][i]
        except Exception as e:
            print(e)
    return inverse

@jit
def invert_bidiagonal(M, type):
    inverse = np.diag(np.ones((M.shape[0])))
    if type == "lower":
        M = np.transpose(M)
        inverse = invert_bidiagonal(M, "upper")
        return np.transpose(inverse)

    elif type == "upper":
        for i in range(M.shape[0] - 1):
            inverse[i][i+1] = -M[i][i+1]
        for k in range(2, M.shape[0]):  
            for i in range(M.shape[0] - k):
                inverse[i][k+i] = inverse[i][k+i-1]*inverse[i+1][k+i]/inverse[i+1][k+i-1]

        return inverse


Ns = [10, 20, 40, 80, 160] # numero de pontos analisados
#Ns = [320]
lambdas = [0.25, 0.5]                 # Constante da exponencial

f = open("saidas.txt", "w+")
implicito = True
for lamb in lambdas:
    for N in Ns:
        print("Calculo para N={}, lambda={}".format(N, lamb))
        f.write("Calculo para N={}, lambda={}\n".format(N, lamb))
        M = (N**2)/lamb         # numero de instantes de tempo analisados
        delta_x = x_total/N     # resolucao espacial do analise
        delta_t = t_total/M     # resolucao temporal da analise

        # Variaveis para acompanhar o andamento do processo

        um_porcento = (int((M-1)*(N-2)))/100 #Calculo da porcentagem
        print("Numero de operacoes: {}\n".format(um_porcento*100))
        f.write("Numero de operacoes: {}\n".format(um_porcento*100))
        cont_por = 0
        ind = '%'
        por_ant = 0

        c = datetime.datetime.now()
        
        contou = False
        matriz = np.zeros((int(M)+1, int(N)+1), dtype='float32') #inicializacao da matriz
        for i in range (0,N): # para todos os x's analisados
            matriz[0][i] = cond_ini(i*delta_x, N)

        for j in range(0, int(M)+1): # para todos os t's analisados
            matriz[j][0] = g1(j*delta_t)
            matriz[j][N] = g2(j*delta_t)
            

        a = datetime.datetime.now()
        if not implicito:
            for i in range (1,int(M)+1): # para cada intervalo de tempo
                t = (i-1)*delta_t
                for j in range (1,N):# para cada intervalo de x
                    x = j*delta_x
                    if item == 'c':
                        fx = funcao_fonte(j, t, N)
                    else:
                        fx = funcao_fonte(x, t, N)
                    matriz[i][j] = matriz[i-1][j] + (delta_t*(((matriz[i-1][j-1]-(2*matriz[i-1][j])+matriz[i-1][j+1])/((delta_x)**2)) + fx))
                    
                    cont_por += (1.0/um_porcento)
                    
                    if (int(cont_por) != int(por_ant)): # Calculo da estimativa de tempo e porcentagem
                        
                        if (contou == False):
                            a = (datetime.datetime.now() - a)*100
                            print("Estimativa de tempo = {}\n".format(a))
                            contou = True
                        
                        print(int(cont_por))
                    
                    por_ant = cont_por
        ############## metodo implicito ###########################
        if implicito:
            for i in range (1,int(M)): # para cada intervalo de tempo
                t = (i-1)*delta_t
                x = j*delta_x

                v = np.zeros((N))
                v[0] = matriz[i-1][0] + delta_t*funcao_fonte(0, i, N) + lamb*g1((i)*delta_t)
                v[N-1] = matriz[i-1][N] + delta_t*funcao_fonte(1, i, N) + lamb*g2((i)*delta_t)
                L, D, Lt = decomporLDL(lamb, N)

                for j in range(1, N-1):
                    v[j] = matriz[i-1][j] + delta_t*funcao_fonte(j, i, N)
                aux = np.dot(invert_bidiagonal(L, "lower"), v)
                aux = np.dot(invert_diagonal(D), aux)
                matriz[i][:N] = np.dot(invert_bidiagonal(Lt, "upper"), aux)


        ###############################################################
        d = datetime.datetime.now()
        
        print("\nO calculo demorou {} segundos".format(d-c))
        f.write("\nO calculo demorou {} segundos".format(d-c))
        matriz = np.transpose(matriz)
        
        ################ plotar matriz ########################
        plt.figure()
        plt.xlabel("Tempo (t)")
        plt.ylabel("Posicao (x)")
        plt.matshow(matriz, fignum = 0, interpolation = 'none', cmap = 'hot', origin = 'lower', aspect="auto")
        plt.savefig("./fig/%s/Resultado-N_%d_Lambda_%f.png" %(item, N,lamb), bbox_inches='tight', dpi=400)
        if show:
            plt.show()
        plt.close()

        ############ gerar matriz com função exata e plotar erro entre matrizes#########
        matriz_ideal = gera_matriz_funcao_exata(M, N, delta_t, delta_x)
        matriz_ideal = np.transpose(matriz_ideal)
        # plt.xlabel("Tempo (t)")
        # plt.ylabel("Posicao (x)")
        # plt.matshow(matriz - matriz_ideal, fignum = 0, interpolation = 'none', cmap = 'Reds', origin = 'lower', aspect="auto")
        # plt.savefig("./fig/%s/IDEAL_%d_Lambda_%f.png" %(item,N,lamb), bbox_inches='tight', dpi=400)
        # if show:
        #     plt.show()
        # plt.close()

        ################################# calcula todos os erros e plota erro em T=1 #################
        plt.plot((matriz-matriz_ideal)[:, int(M)])
        plt.xlabel("Posicao (x)")
        plt.ylabel("Erro em T=1")
        plt.grid()
        plt.savefig("./fig/%s/ErroEm1-N_%d_Lambda_%f.png" %(item,N,lamb), bbox_inches='tight', dpi=400)
        if show:
            plt.show()
        plt.close()

        ############################### plotar solução obtida de 0.1 em 0.1 #########################
        plt.figure(1)
        for t in range(11):
            plt.plot(matriz[:, t*int(M/10)], label="t=%s" %str(t/10))
            if show:
                plt.show()
        plt.grid()
        plt.legend(loc='upper right')
        plt.xlabel("Posicao (x)")
        plt.ylabel("Temperatura")
        plt.savefig("./fig/%s/GraficosSolucao-N_%d_Lambda_%f.png" %(item, N,lamb), bbox_inches='tight', dpi=100)
        plt.close()

f.close()
print("Terminouuu")
