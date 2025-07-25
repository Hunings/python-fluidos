import numpy as np
import matplotlib.pyplot as plt
'''
Este código realiza uma simulação 2D de um fluido em um duto com ressalto utilizando diferenças finitas centradas,
e um solver para a pressão utilizando o método iterativo de Jacobi.

Abaixo são definidas os parâmetros constantes do problema padrão, mas um código auxiliar deve ser utilizado para
rodar uma simulação específica (um exemplo de código está junto a esse no repositório do GitHub duto_ressalto.py).

A simulação principal consiste basicamente na função simulacao(), que inicia os arrays e o loop para a evolução temporal, além do cálculo da pressão,
junto com condições de contorno para velocidade e pressão, definidas como as funções condicoes_contorno_pressao_bfs(p) e condicoes_contorno_velocidades_bfs(u, v).
'''
comprimento = 15
altura = 1
nx = 150
ny = 100
Re = 100
dx = comprimento / (nx-1)
dy = altura / (ny-1)

u_max = 5
v_max = 5
tol = 1e-2
dt = 1e-5
t_final = 1
passos_tempo = 0

it_pressao = 200
plotar_a_cada = 1   

sy = int(nx/7)
sx = int(nx/2)

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
    print(f"Tamanho do degrau {sx}x{sy}")
    input("Pressione qualquer tecla para continuar")
    return 
def condicoes_contorno_pressao_bfs(p):
    # Paredes   
    p[-1, :] = p[-2, :] # leste
    p[:, -1] = p[:, -2] # norte

    # Paredes Duto
    p[sx:, 0] = p[sx:, 1] #dpdy = 0 sul
    p[:sx+1, :sy+1] = 0# p = 0 no ressalto
    p[sx, :sy+1] = p[sx+1, :sy+1] #dpdx = 0,  parede leste ressalto
    p[:sx, sy] = p[:sx, sy+1] #dpdy = 0, parede norte do ressalto

    p[0, sy+1:-1] = p[1, sy+1:-1]

    return p
def condicoes_contorno_velocidades_bfs(u, v):
    u[:, 0] = 0 # sul
    v[:, 0] = 0
    
    u[:, -1] = 0 # norte
    v[:, -1] = 0
    
    # Entrada
    u[0, sy+1:-1] = 1
    v[0, sy+1:-1] = 0
    # Saída
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]

    # Velocidade no interior da borda
    u[:sx+1, :sy+1] = 0
    v[:sx+1, :sy+1] = 0
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

    u_ant, v_ant = condicoes_contorno_velocidades_bfs(u_ant, v_ant)
    p_ant = condicoes_contorno_pressao_bfs(p_ant)

    # Iteração 
    parametros()
    plotar_evolucao = bool(input('Plotar evolução temporal? [0/1]'))
    tt = 0
    for i in range(passos_tempo):
        difusao_x[1:-1, 1:-1] = 1/Re * ((u_ant[2:, 1:-1] - 2*u_ant[1:-1, 1:-1] + u_ant[:-2, 1:-1]) / dx**2 +
                                        (u_ant[1:-1, 2:] - 2*u_ant[1:-1, 1:-1] + u_ant[1:-1, :-2]) / dy**2)
    
        conveccao_x[1:-1, 1:-1] = ((v_ant[1:-1, 2:] * u_ant[1:-1, 2:] - v_ant[1:-1, :-2] * u_ant[1:-1, :-2])/(2*dy) + 
            (u_ant[2:, 1:-1]**2 - u_ant[:-2, 1:-1]**2)/(2*dx))

        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u_ant[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v_ant[2:, 1:-1] - 2*v_ant[1:-1, 1:-1] + v_ant[:-2, 1:-1]) / dx**2 +
                                        (v_ant[1:-1, 2:] - 2*v_ant[1:-1, 1:-1] + v_ant[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v_ant[2:, 1:-1]*u_ant[2:, 1:-1] - v_ant[:-2, 1:-1]*u_ant[:-2, 1:-1])/(2*dx) + (
            v_ant[1:-1, 2:]**2 - v_ant[1:-1, :-2]**2
        )/(2*dy) 

        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v_ant[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

        u_, v_ = condicoes_contorno_velocidades_bfs(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt

        # Resolve a Pressão iterativamente
        p = np.copy(p_ant)
        p_novo = np.copy(p)
        j = 0
        deltap=1
        while j < it_pressao and deltap > tol:
            p_novo[1:-1, 1:-1] = (
              (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
              +
              (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, 1:-1]))/(2 * (dx**2 + dy**2))
            # Condições de Contorno Pressão
            p_novo = condicoes_contorno_pressao_bfs(p_novo)
            deltap = np.max(np.abs(p-p_novo))
            p[:] = p_novo
            j+=1

        dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(2*dy)

        u[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_bfs(u, v)

        u_ant[:], v_ant[:], p_ant[:] = u, v, p
        velocidade_modulo = (u**2 + v**2)**(0.5)
        tt += dt
        V_max = np.max(velocidade_modulo)
        u_maximo = np.max(u)
        v_maximo = np.max(v)
        print('It:', i, '/', passos_tempo, f"t = {tt:.3} / {t_final}", '||u|| =', u_maximo, '||v|| =', v_maximo, '|∆p| =', deltap)
        if V_max > 100 or np.isnan(V_max):
            break
        if plotar_evolucao and i % plotar_a_cada == 0:
            plt.contourf(X, Y, velocidade_modulo, levels=100, cmap='jet')
            plt.colorbar()
            plt.plot(X[:sx+1, sy], Y[:sx+1, sy], c='black')
            plt.plot(X[sx, :sy+1], Y[sx, :sy+1], c='black')
            plt.draw()
            plt.pause(0.005)
            plt.clf()
    plt.show()
    print(u)
    return X, Y, u, v, p, velocidade_modulo
