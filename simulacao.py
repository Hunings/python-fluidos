import numpy as np
import matplotlib.pyplot as plt

#A pressão ainda está interferindo

comprimento = 10
altura = 10
nx = 10
ny = 10
Re = 10
dx = comprimento / (nx-1)
dy = altura / (ny-1)
tol = 1e-2

u_max = 5
v_max = 5
tau = 0.3
dt = tau*min(Re/2*(1/dx**2 + 1/dy**2), dx/u_max, dy/u_max)
passos_tempo = 10000
t_final = 1

parede = int(ny/3)
fator = 2

it_pressao = 100
plotar_a_cada = 1
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
def bifurcacao_velocidades(u, v, fator):
    # Paredes
    u[:, 0] = 0
    u[:, -1] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    # Entrada
    u[0, :] = 0
    v[0, :] = 0
    u[0, (fator-1)*parede:fator*parede+1] = 1
    
    # Saída
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    return u, v
def bifurcacao_pressao(p, fator):
    #Saída
    p[-1, :] = p[-2, :]
    #Topo
    p[:, -1] = p[:, -2]
    #Base
    p[:, 0] = p[:, 1]
    #Entrada
    p[0, :] = p[1, :]
    p[0, (fator-1)*parede:fator*parede+1] = p[1, (fator-1)*parede:fator*parede+1]
    return p
def condicoes_contorno_pressao_duto(p):
    p[-1, :] = p[-2, :]
    p[:, -1] = p[:, -2]
    p[:, 0] = p[:, 1]
    p[0, :] = p[1, :]
    return p
def condicoes_contorno_pressao_cavidade(p):
    p[-1, :] = p[-2, :]
    p[:, -1] = p[:, -2]
    p[:, 0] = p[:, 1]
    p[0, :] = p[1, :] 
    return p

def condicoes_contorno_velocidades_duto(u, v):
   # Paredes
    u[:, 0] = 0
    u[:, -1] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    # Entrada
    u[0, :] = 1
    v[0, :] = 0
    
    # Saída
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    return u, v
def condicoes_contorno_velocidades_cavidade(u, v):
    # Paredes
    u[:, 0] = 0
    u[0, :] = 0
    u[-1, :] = 0
    v[0, :] = 0
    v[-1, :] = 0
    v[:, 0] = 0
    # Parede de entrada de fluido
    u[1:-1, -1] = 1
    v[:, -1] = 0
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
    
        
        dpdx[1:-1, 1:-1] = (p_ant[2:, 1:-1] - p_ant[:-2, 1:-1])/(2*dx) 

        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u_ant[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v_ant[2:, 1:-1] - 2*v_ant[1:-1, 1:-1] + v_ant[:-2, 1:-1]) / dx**2 +
                                        (v_ant[1:-1, 2:] - 2*v_ant[1:-1, 1:-1] + v_ant[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v_ant[2:, 1:-1]*u_ant[2:, 1:-1] - v_ant[:-2, 1:-1]*u_ant[:-2, 1:-1])/(2*dy) + (
            v_ant[1:-1, 2:]**2 - v_ant[1:-1, :-2]**2
        )/(2*dy)
        
        dpdy[1:-1, 1:-1] = (p_ant[1:-1, 2:] - p_ant[1:-1, :-2])/(2*dy)

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

        velocidade_modulo = (u**2 + v**2)**(0.5)
        tt += dt
        V_max = np.max(velocidade_modulo)
        print('It:', i, '/', passos_tempo, f"t = {tt:.3} / {t_final}", '||u|| =', V_max, '|∆p| =', deltap)
        if V_max > 50 or np.isnan(V_max):
            break
        if plotar_evolucao and i % plotar_a_cada == 0:
            plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='viridis')
            plt.colorbar()
            plt.quiver(X, Y, u, v, color='white')
            plt.draw()
            plt.pause(0.005)
            plt.clf()
    plt.show()
    return X, Y, u, v, p, velocidade_modulo
if __name__ == '__main__':
    print(simulacao(1, 0, 1))