import numpy as np
import matplotlib.pyplot as plt

comprimento = 50
altura = 20
nx = 50
ny = 20
Re = 1
dx = comprimento / nx
dy = altura / ny

u_max = 3
v_max = 3
tau = 0.01
dt = tau * min(Re/2*(1/dx**2 + 1/dy**2), dx/u_max, dy/v_max)
passos_tempo = 5000

it_pressao = 100
plotar_a_cada = 1

proporcao_ressalto_duto_comprimento = 0.3
proporcao_ressalto_duto_altura = 0.3
altura_ressalto_pontos = int(proporcao_ressalto_duto_altura * ny)
comprimento_ressalto_pontos = int(proporcao_ressalto_duto_comprimento * nx)

def condicoes_contorno_pressao_duto(p):
    p[-1, :] = p[-2, :]
    p[:, -1] = p[:, -2]
    p[:, 0] = p[:, 1]
    p[0, :] = p[1, :]
    return p
def condicoes_contorno_pressao_cavidade(p):
    p[-1, :] = p[-2, :]
    p[:, -1] = 0
    p[:, 0] = p[:, 1]
    p[0, :] = p[1, :] 
    return p
def condicoes_contorno_pressao_bfs(p):
    p[-1, :] = -p[-2, :]
    p[:, -1] = p[:, -2]
    p[:, 0] = p[:, 1]
    p[0, :] = p[1, :] 
    p[:comprimento_ressalto_pontos, :altura_ressalto_pontos] = 0
    p[:comprimento_ressalto_pontos, altura_ressalto_pontos] = p[:comprimento_ressalto_pontos, altura_ressalto_pontos+1]
    p[comprimento_ressalto_pontos, :altura_ressalto_pontos] = p[comprimento_ressalto_pontos+1, :altura_ressalto_pontos]
    return p

def condicoes_contorno_velocidades_duto(u, v):
   # Paredes
    u[:, 0] = 0
    u[:, -1] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    # Entrada
    u[0, 1:-1] = 1
    u[-1, :] = u[-2, :]

    v[0, :] = 0
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
    u[:, -1] = 1
    v[:, -1] = v[:, -2]
    return u, v
def condicoes_contorno_velocidades_bfs(u, v):
    u[:, 0] = 0
    u[:, -1] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    
    u[0, altura_ressalto_pontos:-1] = 1
    u[-1, :] = u[-2, :]
    
    v[0, altura_ressalto_pontos:-1] = 0
    v[-1, :] = v[-2, :]

    u[:comprimento_ressalto_pontos, :altura_ressalto_pontos] = 0
    v[:comprimento_ressalto_pontos, :altura_ressalto_pontos] = 0

    return u, v
def simulacao(u0, v0, p0):
    # Malha
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)

    # Condições Iniciais
    u_anterior = u0*np.ones((nx, ny))
    v_anterior = v0*np.ones((nx, ny))
    p_anterior = p0*np.ones((nx, ny))

    #u_anterior[:, 0] = -u_anterior[:, 1]
    #u_anterior[:, -1] = -u_anterior[:, -2]

    # Arrays 
    u_, v_, p_ = np.zeros_like(u_anterior), np.zeros_like(v_anterior), np.zeros_like(p_anterior)
    u, v, p = np.zeros_like(u_anterior), np.zeros_like(v_anterior), np.zeros_like(p_anterior)
    difusao_x, difusao_y = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    conveccao_x, conveccao_y = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    dpdx, dpdy = np.zeros_like(p_anterior), np.zeros_like(p_anterior)
    fonte = np.zeros_like(p_anterior)

    # Iteração 
    for i in range(passos_tempo):
        difusao_x[1:-1, 1:-1] = 1/Re * ((u_anterior[2:, 1:-1] - 2*u_anterior[1:-1, 1:-1] + u_anterior[:-2, 1:-1]) / dx**2 +
                                        (u_anterior[1:-1, 2:] - 2*u_anterior[1:-1, 1:-1] + u_anterior[1:-1, :-2]) / dy**2)
        
        conveccao_x[1:-1, 1:-1] = (u_anterior[2:, 1:-1]**2 - u_anterior[:-2, 1:-1]**2)/(2*dx) + (u_anterior[1:-1, 2:]*v_anterior[1:-1, 2:] - u_anterior[1:-1, :-2]*v_anterior[1:-1, :-2])/(2*dy)
        
        dpdx[1:-1, 1:-1] = (p_anterior[2:, 1:-1] - p_anterior[:-2, 1:-1])/(2*dx) # dx ou 2dx?

        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u_anterior[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1] - dpdx[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v_anterior[2:, 1:-1] - 2*v_anterior[1:-1, 1:-1] + v_anterior[:-2, 1:-1]) / dx**2 +
                                        (v_anterior[1:-1, 2:] - 2*v_anterior[1:-1, 1:-1] + v_anterior[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v_anterior[1:-1, 2:]**2 - v_anterior[1:-1, :-2]**2)/(2*dy) + (u_anterior[2:, 1:-1]*v_anterior[2:, 1:-1] - u_anterior[:-2, 1:-1]*v_anterior[:-2, 1:-1])/(2*dx)
        
        dpdy[1:-1, 1:-1] = (p_anterior[1:-1, 2:] - p_anterior[1:-1, :-2])/(2*dy)

        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v_anterior[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1] - dpdy[1:-1, 1:-1])

        # Condições de Contorno velocidades temporárias
        u_, v_ = condicoes_contorno_velocidades_duto(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt # dx ou 2dx?
        
        # Resolve Pressão

        for j in range(it_pressao):
            p_ = np.copy(p)
            p[1:-1, 1:-1] = ((dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]) + dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]) - dx**2 * dy**2 * fonte[1:-1, 1:-1])/ (2 * (dx**2 + dy**2)))

            # Condições de Contorno Pressão
            p = condicoes_contorno_pressao_duto(p)

        p += p_anterior
        dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(2*dy)
        u[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_duto(u, v)

        u_anterior, v_anterior, p_anterior = np.copy(u), np.copy(v), np.copy(p)

        velocidade_modulo = (u**2 + v**2)**(0.5)
        if i % plotar_a_cada == 0:
            plt.contourf(X, Y, velocidade_modulo, levels=10, cmap='inferno')
            plt.colorbar()
            plt.draw()
            plt.pause(0.05)
            plt.clf()
    plt.show()
if __name__ == '__main__':
    print(simulacao(0, 0, 1))