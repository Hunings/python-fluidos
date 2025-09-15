import simulacao as s
import matplotlib.pyplot as plt

comprimento = 1
altura = 1
nx = 51
ny = 51
Re = 10
u_max = 1
v_max = 1
tau = 0.5
t_final = 100

tol = 1e-4
it_pressao = 300
plotar_a_cada = 10

u0 = 1
v0 = p0 = 0

s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_cav

X, Y, u, v, p, V, tempo, x1r = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, forma='quadrado', caso='duto')

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

#s.plotar_contorno(X, Y, V, Re, t_final, nx, ny, comprimento, altura, x1r, tau, it_pressao,  'Módulo da Velocidade', True)
s.plotar_streamlines(X, Y, u, v, Re, t_final, nx, ny, comprimento, altura, x1r, tau, it_pressao, 'streamlines', True)
#s.plotar_vetores(X, Y, u, v, V, Re, t_final, 60, 1, True)