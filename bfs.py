import matplotlib.pyplot as plt
import simulacao as s
import numpy as np


comprimento = 35
altura = 2
nx = 201
ny = 21
Re = 100
u_max = 5
v_max = 5
tau = 0.5
t_final = 10
it_pressao = 1000
tol = 1e-4
plotar_a_cada = 10
s.sx = int((nx-1)/6)
s.sy = int((ny-1)/2)  
u0 = 1
v0 = p0 = 0

s.condicoes_contorno_pressao_duto = s.condicoes_contorno_pressao_bfs
s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_bfs

X, Y, u, v, p, velocidade_modulo, tempo = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

print(f"Tempo de execução: {(tempo):2f} segundos")

dvdx, dudy = np.zeros_like(u), np.zeros_like(v)
dvdx[1:-1, 1:-1] = (v[2:, 1:-1] - v[:-2, 1:-1] ) / (2*s.dx)
dvdx[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2] ) / (2*s.dy)
w = 1/2 * (dvdx - dudy)

s.plotar_contorno(X, Y, velocidade_modulo, Re, t_final, 'Módulo da Velocidade')
s.plotar_streamlines(X, Y, u, v, velocidade_modulo, Re, t_final)
s.plotar_vetores(X, Y, u, v, velocidade_modulo, Re, t_final)