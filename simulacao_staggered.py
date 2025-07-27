import matplotlib.pyplot as plt
import numpy as np

comprimento = 1
altura = 1
nx = 50
ny = 50
Re = 1000
dx = comprimento / (nx-1) # tamanho dividido pelo NÚMERO DE CÉLULAS
dy = altura / (ny-1) # tamanho dividido pelo NÚMERO DE CÉLULAS

u_max = 5
v_max = 5
tol = 1e-2
dt = 1e-2
t_final = 1
passos_tempo = 2000

it_pressao = 100
plotar_a_cada = 1   

sy = int(nx/5)
sx = int(nx/5)

def condicoes_contorno_V(u, v):
    #Sul
    u[:, 0] = -u[:, 1] # No-Slip
    v[:, 0] = 0 # No-Slip

    #Norte
    u[:, -1] = -u[:, -2] # No-slip
    v[:, -1] = 0 # No-slip

    #Leste
    v[-1, :-1] = v[-2, :-1] # Neumann
    u[-1, :-1] = u[-2, :-1] # Neumann

    #Oeste
    v[0, :-1] = -v[1, :-1] # Dirichlet = 0
    u[0, :-1] = 1.0 # Dirichlet

    return u, v
def condicoes_contorno_V_cavidade(u, v):
    #Sul
    u[:, 0] = -u[:, 1] # No-Slip
    v[:, 0] = 0 # No-Slip

    #Norte
    u[:, -1] = 2 - u[:, -2] # Dirichlet = 1
    v[:, -1] = 0 # Dirichlet = 0

    #Leste
    v[-1, :-1] = -v[-2, :-1] # No-slip
    u[-1, :-1] = 0 # No-slip

    #Oeste
    v[0, :-1] = -v[1, :-1] # No-slip
    u[0, :-1] = 0 # No-slip

    return u, v
def condicoes_contorno_p(p):
    #Sul
    p[:, 0] = p[:, 1]
    #Norte
    p[:, -1] = p[:, -2]
    #Leste
    p[-1, :] = p[-2, :]
    #Oeste
    p[0, :] = p[1, :]

    return p

def simulacao(u0, v0, p0):
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)

    u_ini = u0*np.ones((nx, ny+1)) # 7x6
    v_ini = v0*np.ones((nx+1, ny)) #8x5
    p_ini = p0*np.ones((nx+1, ny+1)) #8x6

    u_ini, v_ini = condicoes_contorno_V(u_ini, v_ini)
    p_ini = condicoes_contorno_p(p_ini)

    u, v, p = np.zeros_like(u_ini), np.zeros_like(v_ini), np.zeros_like(p_ini)
    u_, v_ = np.zeros_like(u), np.zeros_like(v)
    u_novo, v_novo = np.zeros_like(u), np.zeros_like(v)

    u = np.copy(u_ini)
    v = np.copy(v_ini)
    p = np.copy(p_ini)

    t=0
    plotar_evolucao = bool(input())
    for it in range(passos_tempo):
        difusao_x = 1/Re * ( # 5x4 nesse caso
            (u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/dx**2
            +
            (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2])/dy**2
        )
        difusao_y = 1/Re * ( # 6x3 nesse caso
            (v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1])/dx**2
            +
            (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, :-2])/dy**2
        )
        du2dx = ((u[1:-1, 1:-1] + u[2:, 1:-1])**2 - (u[1:-1, 1:-1] + u[:-2, 1:-1])**2)/(4*dx) # 5x4
        duvdy = ((v[1:-2, 1:] + v[2:-1, 1:])*(u[1:-1, 1:-1] + u[1:-1, 2:]) - (v[1:-2, :-1] + v[2:-1, :-1])*(u[1:-1, :-2] + u[1:-1, 1:-1]))/(4*dy)
        conveccao_x = du2dx + duvdy
        dv2dy2 = ((v[1:-1, 1:-1] + v[1:-1, 2:])**2 - (v[1:-1, 1:-1] + v[1:-1, :-2])**2)/(4*dy) #6x3
        duvdx = ((u[1:, 1:-2] + u[:-1, 2:-1])*(v[1:-1, 1:-1] + v[2:, 1:-1]) - (u[:-1, 1:-2] + u[:-1, 2:-1])*(v[:-2, 1:-1] + v[1:-1, 1:-1]))/(4*dx)
        conveccao_y = dv2dy2 + duvdx

        u_[1:-1, 1:-1] = u[1:-1, 1:-1] + dt*(difusao_x - conveccao_x) #5x4 considerando os [1:-1, 1:-1]
        v_[1:-1, 1:-1] = v[1:-1, 1:-1] + dt*(difusao_y - conveccao_y) #6x3 considerando os [1:-1, 1:-1]

        u_, v_ = condicoes_contorno_V(u_, v_)

        fonte = ((u_[1:, 1:-1] - u_[:-1, 1:-1])/dx + (v_[1:-1, 1:] - v_[1:-1, :-1])/dy)/dt #6x4  TALVEZ O PROBLEMA ESTEJA AQUI! na verdade com certeza está daqui pra baixo
        
        p_novo = np.copy(p) # pressão é 8x6
        j = 0
        deltap=1
        while j < it_pressao and deltap > tol:
            p_novo[1:-1, 1:-1] = (
                (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
                +
                (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
                -
                (dx**2 * dy**2 * fonte))/(2 * (dx**2 + dy**2))
            # Condições de Contorno Pressão
            deltap = np.max(abs(p-p_novo))
            p_novo = condicoes_contorno_p(p_novo)
            p[:] = p_novo
            j+=1

        dpdx = (p_novo[2:-1, 1:-1] - p_novo[1:-2, 1:-1])/(dx) #7x6
        dpdy = (p_novo[1:-1, 2:-1] - p_novo[1:-1, 1:-2])/(dy) #8x5
        u_novo[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx * dt # 7x6
        v_novo[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy * dt # 8x5
        
        u_novo, v_novo = condicoes_contorno_V(u_novo, v_novo)

        u[:], v[:], p[:] = u_novo, v_novo, p_novo
        
        t += dt

        u_maximo = np.max(u[1:-1, 1:-1])
        v_maximo = np.max(v[1:-1, 1:-1])
        print('It:', it, '/', passos_tempo, f"t = {t:.3} / {t_final}", '||u|| =', u_maximo, '||v|| =', v_maximo, '|∆p| =', deltap)
        um, vm = (u[:, :-1] + u[:, 1:])/2, (v[1:, :] + v[:-1, :])/2
        V = (um**2+vm**2)**(1/2)
        p_med = (p_novo[:-1, 1:] + p_novo[1:, 1:] + p_novo[:-1, :-1] + p_novo[1:, :-1])/4
        if plotar_evolucao:
            plt.contourf(X, Y, V, levels=100, cmap='jet')
            plt.colorbar()
            plt.draw()
            plt.pause(0.1)
            plt.clf()
    plt.show()
    return X, Y, um, vm, V, p_med
X, Y, u, v, V, p = simulacao(1, 0, 0)
plt.figure(figsize=(11, 10))
plt.streamplot(X.T, Y.T, u.T, v.T, color=V.T, cmap='jet')
plt.xlabel('X')
plt.ylabel('Y')  
plt.title(f"Re = {Re} t = {t_final}, staggered")
plt.colorbar()
plt.show()
plt.title('Velocidade horizontal u')
plt.contourf(X, Y, u, levels=200, cmap='jet')
plt.colorbar()
plt.show()
plt.title('Velocidade vertical v')
plt.contourf(X, Y, v, levels=200, cmap='jet')
plt.colorbar()
plt.show()
plt.title('Pressão')
plt.contourf(X, Y, p, levels=200, cmap='jet')
plt.colorbar()
plt.show()