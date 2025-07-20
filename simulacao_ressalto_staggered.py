from numpy import *
from matplotlib.pyplot import *

comprimento = 15
altura = 1
nx = 7
ny = 5
Re = 100
dx = comprimento / (nx-1) # tamanho dividido pelo NÚMERO DE CÉLULAS
dy = altura / (ny-1) # tamanho dividido pelo NÚMERO DE CÉLULAS

u_max = 5
v_max = 5
tol = 1e-2
dt = 1e-5
t_final = 1
passos_tempo = 10

it_pressao = 200
plotar_a_cada = 1   

sy = int(nx/5)
sx = int(nx/5)

def condicoes_contorno_V(u, v):
    #Sul
    u[:, 0] = -u[:, 1] # No-Slip
    v[:, 0] = 0 # No-Slip

    #Norte
    u[:, -1] = -u[:, -2] # No-Slip
    v[:, -1] = 0 # No-Slip

    #Leste
    v[-1, :] = -v[-2, :] # No-Slip
    u[-1, :] = u[-2, :] # Neumann

    #Oeste
    v[1, :] = v[0, :] # Neumann
    u[1, :] = 1.0

    return u, v
def condicoes_contorno_p(p):
    #Sul
    p[:, 1] = p[:, 2]
    #Norte
    p[:, -1] = p[:, -2]
    #Leste
    p[-1, :] = p[-2, :]
    #Oeste
    p[0, :] = p[1, :]

    return p

def simulacao(u0, v0, p0):
    x = linspace(0.0, comprimento, nx)
    y = linspace(0.0, altura, ny)
    X, Y = meshgrid(x, y)
    X, Y = transpose(X), transpose(Y)

    u_ini = u0*ones((nx, ny+1)) # 7x6
    v_ini = v0*ones((nx+1, ny)) #8x5
    p_ini = p0*ones((nx+1, ny+1)) #6x4

    u, v = zeros_like(u_ini), zeros_like(v_ini)
    u_, v_ = zeros_like(u), zeros_like(v)
    fonte = zeros_like(p_ini)
    dpdx = zeros_like(p_ini)
    dpdy = zeros_like(p_ini)

    u = u_ini
    v = v_ini
    p = p_ini
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
        dv2dy2 = ((v[1:-1, 1:-1] + v[1:-1, 2:])**2 - (v[1:-1, 1:-1] + v[1:-1, :-2])**2)/(4*dy)
        duvdx = duvdx = ((u[1:, 1:-2] + u[:-1, 2:-1])*(v[1:-1, 1:-1] + v[2:, 1:-1]) - (u[:-1, 1:-2] + u[:-1, 2:-1])*(v[:-2, 1:-1] + v[1:-1, 1:-1]))/(4*dx)
        conveccao_y = dv2dy2 + duvdx

        u_[1:-1, 1:-1] = u[1:-1, 1:-1] + dt*(difusao_x - conveccao_x)
        v_[1:-1, 1:-1] = v[1:-1, 1:-1] + dt*(difusao_y - conveccao_y)

        print(shape(u_), shape(v_))
        u_, v_ = condicoes_contorno_V(u_, v_)

        fonte = ((u_[1:, 1:-1] - u_[:-1, 1:-1])/(2*dx) + (v_[1:-1, 1:] - v_[1:-1, :-1])/(2*dy))/dt
        print(shape(fonte))
        p_novo = copy(p)
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
            p_novo = condicoes_contorno_p(p_novo)
            deltap = max(abs(p-p_novo))
            p[:] = p_novo
            j+=1
        #dpdx[1:-1, 1:-1] = (p[2:, 1:-1] - p[:-2, 1:-1])/(dx)
        #dpdy[1:-1, 1:-1] = (p[1:-1, 2:] - p[1:-1, :-2])/(dy)

        #u_novo = u_[1:-1, 1:-1] - dpdx * dt
        #v_novo = v_[1:-1, 1:-1] - dpdy * dt

        #u_novo, v_novo = condicoes_contorno_V(u_novo, v_novo)

        #u, v, p = u_novo, v_novo, p_novo

        #t += dt

        #print('It:', it, '/', passos_tempo, f"t = {t:.3} / {t_final}", '||u|| =', max(u), '||v|| =', max(v), '|∆p| =', deltap)
    print(shape(u), shape(v), shape(p), shape(difusao_x), shape(conveccao_x), shape(fonte), shape(dpdx), shape(dpdy))
    u_vert = (u[:, 1:] + u[:, :-1]) / 2
    v_vert = (v[1:, :] + v[:-1, :]) / 2
    p_vert_x = (p_ini[1:, 1:] +p_ini[:-1, 1:]) / 2
    return
simulacao(10, 1, 0)