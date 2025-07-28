import matplotlib.pyplot as plt
import simulacao as s
import numpy as np


s.comprimento = 35
s.altura = 2
s.nx = 201
s.ny = 21
s.Re = 100
s.dx = s.comprimento/(s.nx-1)
s.dy = s.altura/(s.ny-1)
s.u_max = 5
s.v_max = 5
tau = 0.5
s.dt = tau*min(s.Re/2*(1/s.dx**2 + 1/s.dy**2), s.dx/s.u_max, s.dy/s.u_max)

t_final = 1
s.passos_tempo = int(t_final/s.dt)

s.it_pressao = 1000
s.tol = 1e-4
s.plotar_a_cada = 10
s.sx = int((s.nx-1)/6)
s.sy = int((s.ny-1)/2)  

s.condicoes_contorno_pressao_duto = s.condicoes_contorno_pressao_bfs
s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_bfs

X, Y, u, v, p, velocidade_modulo, tempo = s.simulacao(1, 0, 0)

print(f"Tempo de execução: {(tempo):2f} segundos")

dvdx, dudy = np.zeros_like(u), np.zeros_like(v)
dvdx[1:-1, 1:-1] = (v[2:, 1:-1] - v[:-2, 1:-1] ) / (2*s.dx)
dvdx[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2] ) / (2*s.dy)
w = 1/2 * (dvdx - dudy)

s.plotar_contorno(X, Y, velocidade_modulo, s.Re, t_final, 'Módulo da Velocidade')
s.plotar_streamlines(X, Y, u, v, velocidade_modulo, s.Re, t_final)
s.plotar_vetores(X, Y, u, v, velocidade_modulo, s.Re, t_final)