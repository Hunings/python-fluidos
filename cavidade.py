import simulacao as s
import matplotlib.pyplot as plt

comprimento = 1
altura = 1
nx = 30
ny = 30
Re = 100
u_max = 1
v_max = 1
tau = 0.5
t_final = 200

tol = 1e-2
it_pressao = 100
plotar_a_cada = 1

u0 = 0
v0 = p0 = 0

s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_cav

X, Y, u, v, p, velocidade_modulo, tempo = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, 'quadrado')

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

s.plotar_contorno(X, Y, velocidade_modulo, Re, t_final, 'Módulo da Velocidade', 1)
s.plotar_streamlines(X, Y, u, v, velocidade_modulo, Re, t_final, 1)
s.plotar_vetores(X, Y, u, v, velocidade_modulo, Re, t_final, 1, 1)