import automatocelular as ac
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

#tamanho da matriz, n^2 é a população simulada
n = 1000
#fator de mudança da população simulada para o Brasil
constante = float(200000000/(n**2))
#tempo de simulação
t = 366
tempo = np.arange(0,t,1)

#casos: numero de casos na simulação, brasil: número de casos pela mudança de população
brasil = np.zeros([t], dtype=np.float32)
casos = np.zeros([t], dtype=np.int32)
#mortos e mortos brasil é a mesma coisa
mortos = np.zeros([n,n], dtype=np.int32)
mortosbrasil = np.zeros([t], dtype=np.int32)
#variação pura
variacaomortos = np.zeros([t], dtype=np.int32)
#média móvel e tempo da média móvel (plotagem)
mediamovel = np.zeros([t-6], dtype=np.int32)
tmediamovel = np.arange(0,t-6,1)

#geração das matrizes de população e timer (para a ativação do ciclo de infecção)
matriz = ac.automatocelular.geracao_inicial(n)
matriz_timer = ac.automatocelular.geracao_timer(matriz,n)

for i in range(t):
    #calculos da infecção, atualização do timer e mortes ao longo do tempo
    matriz = ac.automatocelular.calculo_infeccao(matriz,matriz_timer,n)
    matriz_timer = ac.automatocelular.atualizacao_timer(matriz,matriz_timer,n)

    casos[i] = ac.automatocelular.soma_infectados(matriz,n)
    brasil[i] = casos[i]*constante

    mortos = ac.automatocelular.atualizacao_mortos(matriz,matriz_timer,mortos,n)
    mortosbrasil[i] = (ac.automatocelular.soma_mortos(mortos,n))*constante

    #calculos de variação de mortes e média móvel
    if i == 0:
        variacaomortos[i] = mortosbrasil[i]
    else:
        variacaomortos[i] = mortosbrasil[i]-mortosbrasil[i-1]
    if i >=6:
        mediamovel[i-6] = (variacaomortos[i]+variacaomortos[i-1]+variacaomortos[i-2]+variacaomortos[i-3]+variacaomortos[i-4]+variacaomortos[i-5]+variacaomortos[i-6])/7

#plots que ajudam a entender melhor e verificar se ta certo
plt.plot(tmediamovel,mediamovel)
plt.xlim(0,t-6)
plt.show()
#plt.plot(tempo,mortosbrasil)
#plt.xlim(20,t)
plt.show()
plt.imshow(mortos, interpolation='none')
plt.show()