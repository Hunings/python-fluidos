import matplotlib.pyplot as plt
import simulacao as s
import numpy as np

comprimento = 24
altura = 2
nx = 201
ny = 21
Re = 200
u_max = 1
v_max = 1
tau = 0.2
t_final = 40
it_pressao = 10000
tol = 1e-2
plotar_a_cada = 1
s.bfs_x = int(nx/6)
s.bfs_y = int(ny/2)  
u0 = 0
v0 = p0 = 0

s.condicoes_contorno_pressao_duto = s.condicoes_contorno_pressao_bfs
s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_bfs

X, Y, u, v, p, velocidade_modulo, tempo, it = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

print(f"Tempo de execução: {(tempo):2f} segundos")

dx = comprimento / (nx-1)
dy = altura / (ny-1)

dvdx, dudy = np.zeros_like(u), np.zeros_like(v)
dvdx[1:-1, 1:-1] = (v[2:, 1:-1] - v[:-2, 1:-1] ) / (2*dx)
dvdx[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2] ) / (2*dy)
w = 1/2 * (dvdx - dudy)

s.plotar_contorno(X, Y, velocidade_modulo, Re, t_final, 'Módulo da Velocidade', False)
s.plotar_streamlines(X, Y, u, v, velocidade_modulo, Re, t_final, False)
#s.plotar_vetores(X, Y, u, v, velocidade_modulo, Re, t_final, 100, 2, False)
s.plotar_contorno(X, Y, p, Re, t_final, 'Vorticidade', False)