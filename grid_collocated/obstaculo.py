import simulacao as s
import matplotlib.pyplot as plt

#Define as constantes da simulação
comprimento = 15
altura = 1
nx = 300
ny = 150
Re = 300
u_max = 3
v_max = 3
tau = 0.5
t_final = 100
it_pressao = 100
plotar_a_cada = 100
tol = 1e-3
u0, v0, p0 = 1, 0, 0
s.pos_obs_x = 10
s.pos_obs_y = 10
s.tam_obs = 10


# Modifica as condições de contorno internas

s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_obs
s.condicoes_contorno_pressao_duto = s.condicoes_contorno_pressao_obs

#Executa a simulação
X, Y, u, v, p, V, tempo = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

#Visualização 

s.plotar_contorno(X, Y, V, Re, t_final, 'Módulo da Velocidade', 0)
s.plotar_streamlines(X, Y, u, v, V, Re, t_final, 0)
s.plotar_vetores(X, Y, u, v, V, Re, t_final, 30, 1, 0)