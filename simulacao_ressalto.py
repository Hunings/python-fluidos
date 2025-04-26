import numpy as np
import matplotlib.pyplot as plt

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

alt_bfs = 5 # em número de pontos
comp_bfs = 5

def condicoes_contorno_pressao_bfs(p):
    # Paredes 
    p[-1, :] = 0
    p[:, -1] = p[:, -2]

    # Paredes Duto
    p[comp_bfs:, 0] = p[comp_bfs:, 1]
    p[:comp_bfs+1, :alt_bfs+1] = 0 # p = 0 no ressalto
    p[comp_bfs, :alt_bfs+1] = p[comp_bfs+1, :alt_bfs+1] #dpdx = 0,  +1 no alt_bfs porque :alt_bfs vai até alt_bfs-1
    p[:comp_bfs+1, alt_bfs] = p[:comp_bfs+1, alt_bfs] #dpdy = 0
    p[comp_bfs:, 0] = p[comp_bfs:, 1]

    p[0, alt_bfs+1:-1] = p[1, alt_bfs+1:-1]

    return p
def condicoes_contorno_velocidades_bfs(u, v):
    u[:, 0] = 0
    u[:, -1] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    
    # Entrada
    u[0, alt_bfs+1:-1] = 1
    v[0, alt_bfs+1:] = 0
    # Saída
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]

    # Velocidade no interior da borda
    u[:comp_bfs+1, :alt_bfs+1] = 0
    v[:comp_bfs+1, :alt_bfs+1] = 0

    return u, v
    
def simulacao(u0, v0, p0):
    # Malha
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X_, Y_ = np.meshgrid(x, y)
    X, Y = np.transpose(X_), np.transpose(Y_)

    # Condições Iniciais
    u_ant = u0*np.ones((nx, ny))
    v_ant = v0*np.ones((nx, ny))
    p_ant = p0*np.ones((nx, ny))

    # Arrays 
    u_, v_ = np.zeros_like(u_ant), np.zeros_like(v_ant)
    u, v = np.zeros_like(u_ant), np.zeros_like(v_ant)
    difusao_x, difusao_y = np.zeros_like(v_ant), np.zeros_like(v_ant)
    conveccao_x, conveccao_y = np.zeros_like(u_ant), np.zeros_like(v_ant)
    dpdx, dpdy = np.zeros_like(p_ant), np.zeros_like(p_ant)
    fonte = np.zeros_like(p_ant)
    mascara = np.ones_like(p_ant)
    mascara[:comp_bfs, :alt_bfs] = 0

    u_ant, v_ant = condicoes_contorno_velocidades_bfs(u_ant, v_ant)
    p_ant = condicoes_contorno_pressao_bfs(p_ant)

    # Iteração 
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
        
        conveccao_y[1:-1, 1:-1] = (v_ant[2:, 1:-1]*u_ant[2:, 1:-1] - v_ant[:-2, 1:-1]*u_ant[:-2, 1:-1])/(2*dy) + (
            v_ant[1:-1, 2:]**2 - v_ant[1:-1, :-2]**2
        )/(2*dy) 

        #Velocidade v antes da correção da pressão
        v_[1:-1, 1:-1] = v_ant[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

        u_, v_ = condicoes_contorno_velocidades_bfs(u_, v_)

        fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt*mascara[1:-1, 1:-1]

        # Resolve a Pressão iterativamente
        p = np.copy(p_ant)
        for j in range(it_pressao):
            p_old = np.copy(p)
            p[1:-1, 1:-1] = np.where(mascara[1:-1, 1:-1], ((
              dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1])
              +
              dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])
              -
              dx**2 * dy**2 * fonte[1:-1, 1:-1])
              /
            (2 * (dx**2 + dy**2))
            ), p_old[1:-1, 1:-1])
            # Condições de Contorno Pressão
            p = condicoes_contorno_pressao_bfs(p)

        dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*dx)
        dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(2*dy)

        u[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
        v[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

        # Condições de Contorno Velocidades Finais
        u, v = condicoes_contorno_velocidades_bfs(u, v)

        u_ant, v_ant, p_ant = u, v, p

        velocidade_modulo = (u**2 + v**2)**(0.5)
        if i % plotar_a_cada == 0:
            print(i)
            #plt.plot(X[:comp_bfs+1, alt_bfs], Y[:comp_bfs+1, alt_bfs], c='black')
            #plt.plot(X[comp_bfs, :alt_bfs+1], Y[comp_bfs, :alt_bfs+1], c='black')
            plt.contourf(X, Y, velocidade_modulo, levels=60, cmap='viridis')
            plt.colorbar()
            plt.plot(X[:comp_bfs+1, alt_bfs], Y[:comp_bfs+1, alt_bfs], c='black')
            plt.plot(X[comp_bfs, :alt_bfs+1], Y[comp_bfs, :alt_bfs+1], c='black')
            #plt.quiver(X[::4, ::4], Y[::4, ::4], u[::4, ::4], v[::4, ::4], color='white')
            plt.draw()
            plt.pause(0.005)
            plt.clf()
    plt.show()
    plt.contourf(X, Y, velocidade_modulo, levels=10)
    plt.plot(X[:comp_bfs+1, alt_bfs], Y[:comp_bfs+1, alt_bfs], c='black')
    plt.plot(X[comp_bfs, :alt_bfs+1], Y[comp_bfs, :alt_bfs+1], c='black')
    plt.quiver(X[::4, ::4], Y[::4, ::4], u[::4, ::4], v[::4, ::4], color='white')
    plt.colorbar()
    plt.show()
    plt.streamplot(X.T, Y.T, u.T, v.T, cmap='viridis', density=2)
    plt.colorbar()
    plt.show()
if __name__ == '__main__':
    print(simulacao(1, 0, 0))
