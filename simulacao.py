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

u_max = 5
v_max = 5
tau = 0.3
dt = tau*min(Re/2*(1/dx**2 + 1/dy**2), dx/u_max, dy/u_max)
passos_tempo = 10000

it_pressao = 100
plotar_a_cada = 1

altura_ressalto_pontos = 3
comprimento_ressalto_pontos = 3

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
    u[0, 1:-1] = 1
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
    u[:, -1] = 1
    v[:, -1] = v[:, -2]
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

    # Arrays 
    u_, v_ = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    u, v = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    difusao_x, difusao_y = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    conveccao_x, conveccao_y = np.zeros_like(u_anterior), np.zeros_like(v_anterior)
    dpdx, dpdy = np.zeros_like(p_anterior), np.zeros_like(p_anterior)
    fonte = np.zeros_like(p_anterior)

    u_anterior, v_anterior = condicoes_contorno_velocidades_duto(u_anterior, v_anterior)
    p_anterior = condicoes_contorno_pressao_duto(p_anterior)

    # Iteração 
    for i in range(passos_tempo):
        difusao_x[1:-1, 1:-1] = 1/Re * ((u_anterior[2:, 1:-1] - 2*u_anterior[1:-1, 1:-1] + u_anterior[:-2, 1:-1]) / dx**2 +
                                        (u_anterior[1:-1, 2:] - 2*u_anterior[1:-1, 1:-1] + u_anterior[1:-1, :-2]) / dy**2)
    
        conveccao_x[1:-1, 1:-1] = (v_anterior[1:-1, 2:] * u_anterior[1:-1, 2:] - v_anterior[1:-1, :-2] * u_anterior[1:-1, :-2])/(2*dy) + (
            u_anterior[2:, 1:-1]**2 - u_anterior[:-2, 1:-1]**2
        )/(2*dx)
    
        
        dpdx[1:-1, 1:-1] = (p_anterior[2:, 1:-1] - p_anterior[:-2, 1:-1])/(2*dx) 

        # Velocidade u antes da correção da pressão
        u_[1:-1, 1:-1] = u_anterior[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

        difusao_y[1:-1, 1:-1] = 1/Re * ((v_anterior[2:, 1:-1] - 2*v_anterior[1:-1, 1:-1] + v_anterior[:-2, 1:-1]) / dx**2 +
                                        (v_anterior[1:-1, 2:] - 2*v_anterior[1:-1, 1:-1] + v_anterior[1:-1, :-2]) / dy**2)
        
        conveccao_y[1:-1, 1:-1] = (v_anterior[2:, 1:-1]*u_anterior[2:, 1:-1] - v_anterior[:-2, 1:-1]*u_anterior[:-2, 1:-1])/(2*dy) + (
            v_anterior[1:-1, 2:]**2 - v_anterior[1:-1, :-2]**2
        )/(2*dy)
        
        dpdy[1:-1, 1:-1] = (p_anterior[1:-1, 2:] - p_anterior[1:-1, :-2])/(2*dy)

        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v_anterior[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

        u_, v_ = condicoes_contorno_velocidades_duto(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt

        # Resolve a Pressão iterativamente
        p = np.copy(p_anterior)
        for j in range(it_pressao):
            p_old = np.copy(p)
            p[1:-1, 1:-1] = (
          (
              dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1])
              +
              dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])
              -
              dx**2 * dy**2 * fonte[1:-1, 1:-1]
          )
          /
          (2 * (dx**2 + dy**2))
          )
            # Condições de Contorno Pressão
            p = condicoes_contorno_pressao_duto(p)
            erro = np.linalg.norm((p - p_old)/p, ord=np.inf)

        dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(2*dy)
       
        u[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_duto(u, v)

        u_anterior, v_anterior, p_anterior = u, v, p

        velocidade_modulo = (u**2 + v**2)**(0.5)
        if i % plotar_a_cada == 0:
            print(i)
            plt.pcolormesh(X, Y, velocidade_modulo, levels=60, cmap='viridis')
            plt.colorbar()
            plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2], color='white')
            plt.draw()
            plt.pause(0.005)
            plt.clf()
    plt.show()
    plt.contourf(X, Y, velocidade_modulo, levels=10)
    plt.quiver(X, Y, u, v)
    plt.colorbar()
    plt.show()
    print(X[comprimento_ressalto_pontos, altura_ressalto_pontos])
if __name__ == '__main__':
    print(simulacao(1, 0, 1))