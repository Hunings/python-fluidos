import simulacao as s
import matplotlib.pyplot as plt

#Define as constantes da simulação
comprimento = 15
altura = 1
nx = 100
ny = 50
Re = 100
u_max = 2
v_max = 2
tau = 0.9
t_final = 100
tol = 1e-4
it_pressao = 100
plotar_a_cada = 10
u0, v0, p0 = 1, 0, 0

#Condições de contorno internas padrão

#Executa a simulação
X, Y, u, v, p, V, tempo = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

#Visualização 

s.plotar_contorno(X, Y, V, Re, t_final, 'Módulo da Velocidade', 0)
s.plotar_streamlines(X, Y, u, v, V, Re, t_final, 0)
s.plotar_vetores(X, Y, u, v, V, Re, t_final, 60, 2, 0)
