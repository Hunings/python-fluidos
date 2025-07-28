import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
'''
Este código realiza uma simulação 2D de um fluido em um duto utilizando diferenças finitas centradas,
e um solver para a pressão utilizando o método iterativo de Jacobi.

Abaixo são definidas os parâmetros constantes do problema padrão, mas um código auxiliar deve ser utilizado para
rodar uma simulação específica (um exemplo de código está junto a esse no repositório do GitHub duto_ressalto.py).

A simulação principal consiste basicamente na função simulacao(), que inicia os arrays e o loop para a evolução temporal, além do cálculo da pressão,
junto com condições de contorno para velocidade e pressão, definidas como as funções condicoes_contorno_pressao(p) e condicoes_contorno_velocidades(u, v).
'''
comprimento = 0
altura = 0
nx = 0
ny = 0
Re = 0
dx = 0
dy = 0
tol = 0
u_max = 0
v_max = 0
tau = 0
dt = 0
passos_tempo = 0
t_final = 0
parede = 0
fator = 0
it_pressao = 0
plotar_a_cada = 0
pos_x = 0
pos_y = 0
sx = 0
sy = 0

def parametros():
    print("Parâmetros da simulação a seguir")
    print("Altura: ", altura)
    print("Comprimento: ", comprimento)
    print("Número de pontos x: ", nx)
    print("Número de pontos y: ", ny)
    print("Número de Reynolds: ", Re)
    print("Espaçamento x: ", dx)
    print("Espaçamento y: ", dy)
    print("Passo de tempo: ", dt)
    print("Tolerância na diferença de pressão: ", tol)
    print("Velocidade máxima em x: ", u_max)
    print("Velocidade máxima em y: ", v_max)
    print("Total de iterações no tempo: ", passos_tempo)
    print("Tamanho da entrada bifurcação em pontos (se simulando bifurcação): ", parede)
    print("Tempo: ", t_final)
    input("Pressione qualquer tecla para continuar")
    return 
def condicoes_contorno_pressao_bfs(p):
    #Paredes Neumann homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    #Paredes Duto Neumann homogênea
    p[:sx+1, :sy+1] = np.pi #Qualquer valor não deve interferir no resto da simulação
    p[sx, :sy+1] = p[sx+1, :sy+1] #Parede leste do ressalto
    p[:sx+1, sy] = p[:sx+1, sy+1] #Parede norte do ressalto
    return p
def condicoes_contorno_velocidades_bfs(u, v):
    #Paredes No-Slip
    u[:, 0] = 0 # Sul
    v[:, 0] = 0
    u[:, -1] = 0 # Norte
    v[:, -1] = 0
    #Entrada Dirichlet
    u[0, sy+1:-1] = 1. # Oeste
    v[0, sy+1:-1] = 0.
    #Saída Neumann homogênea
    u[-1, :] = u[-2, :] # Leste
    v[-1, :] = v[-2, :]
    #Velocidade no interior da borda
    u[:sx+1, :sy+1] = 0
    v[:sx+1, :sy+1] = 0
    return u, v
def condicoes_contorno_velocidades_obs(u, v):
    #Paredes No-Slip
    u[:, -1] = 0 #Norte
    v[:, -1] = 0
    u[:, 0] = 0 #Sul
    v[:, 0] = 0
    # Entrada Dirichlet
    u[0, 1:-1] = 1.0 #Oeste
    v[0, 1:-1] = 0.0
    # Saída Neumann homogênea
    u[-1, 1:-1] = u[-2, 1:-1] #Leste
    v[-1, 1:-1] = v[-2, 1:-1]
    #Obstáculo No-Slip
    u[pos_x:pos_x+8, pos_y:pos_y+8] = 0
    v[pos_x:pos_x+8, pos_y:pos_y+8] = 0
    return u, v
def condicoes_contorno_pressao_obs(p):
    #Paredes Neumann homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    #Obstáculo Neumann
    p[pos_x:pos_x+8, pos_y+8] = p[pos_x:pos_x+8, pos_y+9] #Norte
    p[pos_x:pos_x+8, pos_y] = p[pos_x:pos_x+8, pos_y-1] #Sul
    p[pos_x, pos_y:pos_y+8] = p[pos_x-1, pos_y:pos_y+8] #Oeste
    p[pos_x+8, pos_y:pos_y+8] = p[pos_x+9, pos_y:pos_y+8] #Leste
    #Interior do obstáculo
    p[pos_x:pos_x+8, pos_y:pos_y+8] = np.pi
    return p
def condicoes_contorno_velocidades_bif(u, v, fator):
    #Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0
    u[:, -1] = 0 #Norte
    v[:, -1] = 0
    u[0, :] = 0 #Oeste
    v[0, :] = 0
    #Entrada Dirichlet
    u[0, (fator-1)*parede:fator*parede+1] = 1
    #Saída Neumann homogênea
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    return u, v
def condicoes_contorno_pressao_bif(p, fator):
    #Paredes Neumman homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    #Entrada Neumann homogênea
    p[0, (fator-1)*parede:fator*parede+1] = p[1, (fator-1)*parede:fator*parede+1]
    return p
def condicoes_contorno_pressao_duto(p):
    #Paredes Neumann homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    return p
def condicoes_contorno_velocidades_duto(u, v):
    #Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0 
    u[:, -1] = 0 #Norte
    v[:, -1] = 0
    #Entrada Dirichlet
    u[0, :] = 1. #Oeste
    v[0, :] = 0.
    #Saída Neumann homogênea
    u[-1, :] = u[-2, :] #Leste
    v[-1, :] = v[-2, :]
    return u, v
def condicoes_contorno_velocidades_cav(u, v):
    # Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0
    u[0, :] = 0 #Oeste
    v[0, :] = 0
    u[-1, :] = 0 #Leste
    v[-1, :] = 0
    #Entrada Dirichlet
    u[1:-1, -1] = 1. #Norte
    v[:, -1] = 0.
    return u, v
def pressao(fonte, p, p_novo, dx, dy):
    j = 0
    deltap = normal2 = 1
    while j < it_pressao and normal2 > tol:
            #Atualiza pressão iterativamente
            p_novo[1:-1, 1:-1] = (
              (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
              +
              (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, 1:-1]))/(2 * (dx**2 + dy**2))
            #Condições de Contorno Pressão
            p_novo = condicoes_contorno_pressao_duto(p_novo)

            #Calcula medidas do erro
            deltap = np.max(abs(p-p_novo))
            normal2 = np.linalg.norm(p-p_novo, ord=np.inf)/np.linalg.norm(p_novo)

            #Atualiza pressão e avança na iteração
            p[:] = p_novo
            j+=1
    return p, deltap, normal2
def malha(comprimento, altura, nx, ny):
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)
    return X, Y
def plotar_contorno(X, Y, V, Re, t_final, titulo):
    plt.figure(figsize=(15, 2))
    plt.contourf(X, Y, V, levels=200, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"{titulo}, Re = {Re}, t = {t_final}")
    plt.colorbar(orientation='horizontal')
    plt.show()
    return
def plotar_streamlines(X, Y, u, v, V, Re, t_final):
    plt.figure(figsize=(15, 2))
    plt.streamplot(X.T, Y.T, u.T, v.T, color=V.T, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"Re = {Re}, t = {t_final}")
    plt.colorbar(orientation='horizontal')
    plt.show()
    return
def plotar_vetores(X, Y, u, v, V, Re, t_final):
    plt.figure(figsize=(15, 2))
    plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2], V[::2, ::2], scale=60, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"Re = {Re}, t = {t_final}")
    plt.colorbar(orientation='horizontal')
    plt.show()  
    return
def informacoes(i, passos_tempo, tempo_transcorrido, t_final, V_max, deltap, normal2):
    print('It:', i, '/', passos_tempo, f"t = {tempo_transcorrido:.3} / {t_final}", '||V|| =', V_max, '|∆p| =', deltap, 'Norma-L2 =', normal2)
    return
def evolucao(X, Y, V):
    plt.contourf(X, Y, V, levels=10, cmap='jet')
    plt.colorbar()
    plt.pause(0.0001)
    plt.clf()
    return
def simulacao(u0, v0, p0):
    #Conta tempo de simulação
    inicio = perf_counter()

    # Malha
    X, Y = malha(comprimento, altura, nx, ny)

    # Condições Iniciais
    u = u0*np.ones((nx, ny))
    v = v0*np.ones((nx, ny))
    p = p0*np.ones((nx, ny))

    # Arrays 
    u_, v_ = np.zeros((nx, ny)), np.zeros((nx, ny))
    u_novo, v_novo = np.zeros((nx, ny)), np.zeros((nx, ny))
    difusao_x, difusao_y = np.zeros((nx, ny)), np.zeros((nx, ny))
    conveccao_x, conveccao_y = np.zeros((nx, ny)), np.zeros((nx, ny))
    dpdx, dpdy = np.zeros((nx, ny)), np.zeros((nx, ny))
    fonte = np.zeros((nx, ny))

    u, v = condicoes_contorno_velocidades_duto(u, v)
    p = condicoes_contorno_pressao_duto(p)

    parametros()
    tempo_transcorrido = 0
    plotar_evolucao = bool(input('Plotar evolução temporal? [0/1]'))
    if plotar_evolucao:
        plt.figure(figsize=(15,2))        
    # Iteração 
    for i in range(passos_tempo):
        difusao_x[1:-1, 1:-1] = 1/Re * ((u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx**2 +
                                        (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2]) / dy**2)
    
        conveccao_x[1:-1, 1:-1] = (v[1:-1, 2:] * u[1:-1, 2:] - v[1:-1, :-2] * u[1:-1, :-2])/(2*dy) + (
            u[2:, 1:-1]**2 - u[:-2, 1:-1]**2
        )/(2*dx)
    
        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1]) / dx**2 +
                                        (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v[2:, 1:-1]*u[2:, 1:-1] - v[:-2, 1:-1]*u[:-2, 1:-1])/(2*dx) + (
            v[1:-1, 2:]**2 - v[1:-1, :-2]**2
        )/(2*dy)
        
        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

        u_, v_ = condicoes_contorno_velocidades_duto(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt

        # Resolve a Pressão iterativamente
        p_novo = np.copy(p)

        p, deltap, normal2 = pressao(fonte, p, p_novo, dx, dy)

        dpdx[1:-1, 1:-1] = (p_novo[2:, 1:-1] - p_novo[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p_novo[1:-1, 2:] - p_novo[1:-1, :-2])/(2*dy)
       
        u_novo[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v_novo[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_duto(u, v)

        u, v, p = u_novo, v_novo, p_novo

        tempo_transcorrido += dt

        V = (u**2 + v**2)**(0.5)
        V_max = np.max(V)

        informacoes(i, passos_tempo, tempo_transcorrido, t_final, V_max, deltap, normal2)

        if V_max > 50 or np.isnan(V_max):
            break
        if plotar_evolucao and i % plotar_a_cada == 0:
            evolucao(X, Y, V)
    fim = perf_counter()
    tempo = fim - inicio
    return X, Y, u, v, p, V, tempo