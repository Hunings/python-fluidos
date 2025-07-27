import numpy as np
import matplotlib.pyplot as plt
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

def simulacao(u0, v0, p0):
    # Malha
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)

    # Condições Iniciais
    u_ant = u0*np.ones((nx, ny))
    v_ant = v0*np.ones((nx, ny))
    p_ant = p0*np.ones((nx, ny))

    # Arrays 
    u_, v_ = np.zeros_like(u_ant), np.zeros_like(v_ant)
    u, v = np.zeros_like(u_ant), np.zeros_like(v_ant)
    difusao_x, difusao_y = np.zeros_like(u_ant), np.zeros_like(v_ant)
    conveccao_x, conveccao_y = np.zeros_like(u_ant), np.zeros_like(v_ant)
    dpdx, dpdy = np.zeros_like(p_ant), np.zeros_like(p_ant)
    fonte = np.zeros_like(p_ant)

    u_ant, v_ant = condicoes_contorno_velocidades_duto(u_ant, v_ant)
    p_ant = condicoes_contorno_pressao_duto(p_ant)

    # Iteração 
    parametros()
    plotar_evolucao = bool(input('Plotar evolução temporal? [0/1]'))
    tt = 0
    for i in range(passos_tempo):
        difusao_x[1:-1, 1:-1] = 1/Re * ((u_ant[2:, 1:-1] - 2*u_ant[1:-1, 1:-1] + u_ant[:-2, 1:-1]) / dx**2 +
                                        (u_ant[1:-1, 2:] - 2*u_ant[1:-1, 1:-1] + u_ant[1:-1, :-2]) / dy**2)
    
        conveccao_x[1:-1, 1:-1] = (v_ant[1:-1, 2:] * u_ant[1:-1, 2:] - v_ant[1:-1, :-2] * u_ant[1:-1, :-2])/(2*dy) + (
            u_ant[2:, 1:-1]**2 - u_ant[:-2, 1:-1]**2
        )/(2*dx)
    
        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u_ant[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v_ant[2:, 1:-1] - 2*v_ant[1:-1, 1:-1] + v_ant[:-2, 1:-1]) / dx**2 +
                                        (v_ant[1:-1, 2:] - 2*v_ant[1:-1, 1:-1] + v_ant[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v_ant[2:, 1:-1]*u_ant[2:, 1:-1] - v_ant[:-2, 1:-1]*u_ant[:-2, 1:-1])/(2*dx) + (
            v_ant[1:-1, 2:]**2 - v_ant[1:-1, :-2]**2
        )/(2*dy)
        
        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v_ant[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

        u_, v_ = condicoes_contorno_velocidades_duto(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt

        # Resolve a Pressão iterativamente
        p = np.copy(p_ant)
        p_novo = np.copy(p)
        j = 0
        deltap = 1
        while j < it_pressao and deltap > tol:
            p_novo[1:-1, 1:-1] = (
              (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
              +
              (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, 1:-1]))/(2 * (dx**2 + dy**2))
            # Condições de Contorno Pressão
            p_novo = condicoes_contorno_pressao_duto(p_novo)
            deltap = np.max(np.abs(p-p_novo))
            p[:] = p_novo
            j+=1

        dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(2*dy)
       
        u[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_duto(u, v)

        u_ant, v_ant, p_ant = u, v, p

        V = (u**2 + v**2)**(0.5)
        tt += dt
        V_max = np.max(V)
        print('It:', i, '/', passos_tempo, f"t = {tt:.3} / {t_final}", '||u|| =', V_max, '|∆p| =', deltap)
        if V_max > 50 or np.isnan(V_max):
            break
        if plotar_evolucao and i % plotar_a_cada == 0:
            plt.contourf(X, Y, V, levels=10, cmap='jet')
            plt.colorbar()
            plt.draw()
            plt.pause(0.005)
            plt.clf()
    plt.show()
    return X, Y, u, v, p, V