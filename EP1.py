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

item = 'b'

@jit
def funcao_fonte(x, t): # funcao de  entrada
    global item
    #return 10*(x**2)*(x-1) - 60*x*t + 20*t
    if item == 'a':
        return 10*(np.cos(10*t))*(x**2)*((1-x)**2) - (1+(np.sin(10*t)))*(12*(x**2)-12*x+2) # item a
    if item == 'b':
        return (np.exp(t-x)*(np.cos(5*t*x)*(11*t-(25*t*t)+5*x-1)))                         # item b

@jit
def funcao_exata(x, t):
    #return (10*t*x*x*(x-1))
    return ((1 + np.sin(10*t))*(x**2)*((1-x)**2))

@jit
def gera_matriz_funcao_exata(M, N, delta_t, delta_x):
    f_exata = np.zeros((int(M), int(N)), dtype='float64')
    for i in range(1, int(M)):
        t = (i-1)*delta_t
        for j in range(1, int(N-1)):
            x = j*delta_x
            f_exata[i][j] = funcao_exata(x, t)
    return f_exata

@jit
def cond_ini (x):
    #return ((x**2)*((1-x)**2))
    return np.exp(-x)       

@jit
def g1(t): # condição de contorno, x=0
    return (np.exp(t))

@jit
def g2(t): # condição de contorno, x=1
    return (np.exp(t-1)*np.cos(5*t))


N = 80         # numero de pontos analisados
lamb = 0.5     # Constante da exponencial

Ns = [10, 20, 40, 80, 160, 320]
lambdas = [0.25, 0.5]

# Ns = [10]
# lambdas = [0.5]
#Ns = [160]
#lambdas = [0.5]
f = open("saidas.txt", "w+")

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
        matriz = np.zeros((int(M), int(N)), dtype='float64') #inicializacao da matriz
        for i in range (0,N):
            matriz[0][i] = cond_ini(i*delta_x)

        for j in range(0, int(M)):
            matriz[j][0] = g1(j*delta_t)
            matriz[j][N-1] = g2(j*delta_t)
            

        a = datetime.datetime.now()
        
        for i in range (1,int(M)): # para cada intervalo de tempo
            t = (i-1)*delta_t
            for j in range (1,N-1):# para cada intervalo de x
                x = j*delta_x
                fx = funcao_fonte(x, t)
                matriz[i][j] = matriz[i-1][j] + (delta_t*(((matriz[i-1][j-1]-(2*matriz[i-1][j])+matriz[i-1][j+1])/((delta_x)**2))+fx))
                
                cont_por += (1.0/um_porcento)
                
                if (int(cont_por) != int(por_ant)):
                    print(int(cont_por))
                    if (contou == False):
                        a = (datetime.datetime.now() - a)*100
                        print("Estimativa de tempo = {}\n".format(a))
                        contou = True
                
                por_ant = cont_por
                
        d = datetime.datetime.now()
        print("\nO calculo demorou {} segundos".format(d-c))
        f.write("\nO calculo demorou {} segundos".format(d-c))
        #print(matriz[int(M-1)][int(N/2)]) ###########################################################################
        matriz = np.transpose(matriz)
        plt.figure()
        plt.xlabel("Tempo (t)")
        plt.ylabel("Posicao (x)")
        plt.matshow(matriz, fignum = 0, interpolation = 'none', cmap = 'hot', origin = 'lower', aspect="auto")
        plt.savefig("./fig/Resultado-N_%d_Lambda_%f.jpg" %(N,lamb), bbox_inches='tight', dpi=400)
        if show:
            plt.show()
        plt.close()


        matriz_ideal = gera_matriz_funcao_exata(M, N, delta_t, delta_x)
        matriz_ideal = np.transpose(matriz_ideal)
        plt.xlabel("Tempo (t)")
        plt.ylabel("Posicao (x)")
        plt.matshow(matriz-matriz_ideal, fignum = 0, interpolation = 'none', cmap = 'Reds', origin = 'lower', aspect="auto")
        #plt.matshow(matriz_ideal, fignum = 0, interpolation = 'none', cmap = 'hot', origin = 'lower', aspect="auto")
        plt.savefig("./fig/IDEAL_%d_Lambda_%f.jpg" %(N,lamb), bbox_inches='tight', dpi=400)
        if show:
            plt.show()
        plt.close()


        plt.plot((matriz-matriz_ideal)[:, int(M)-1])
        plt.xlabel("Posicao (x)")
        plt.ylabel("Erro em T=1")
        plt.grid()
        plt.savefig("./fig/ErroEm1-N_%d_Lambda_%f.jpg" %(N,lamb), bbox_inches='tight', dpi=400)
        if show:
            plt.show()
        plt.close()

f.close()
print("Terminouuu")
