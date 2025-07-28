import simulacao as s
import matplotlib.pyplot as plt
from functools import partial

#Define as constantes da simulação
comprimento = 10
altura = 1
nx = 200
ny = 21
Re = 100
s.abertura = int(ny/3)
u_max = 10
v_max = 10
tol = 1e-2
tau = 0.5
t_final = 10
it_pressao = 100
plotar_a_cada = 100

u0 = 1
v0 = 0
p0 = 0

# Modifica as condições de contorno internas

s.condicoes_contorno_velocidades_duto = s.condicoes_contorno_velocidades_bif
s.condicoes_contorno_pressao_duto = s.condicoes_contorno_pressao_bif

#Executa a simulação
X, Y, u, v, p, V, tempo = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

# Tempo
print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

s.plotar_contorno(X, Y, V, Re, t_final, 'Módulo da Velocidade', False)
s.plotar_streamlines(X, Y, u, v, V, Re, t_final, False)
s.plotar_vetores(X, Y, u, v, V, Re, t_final, 40, 1, 0)
