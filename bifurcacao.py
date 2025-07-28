import simulacao as sim
import matplotlib.pyplot as plt
from functools import partial

#Define as constantes da simulação
comprimento = 10
altura = 1
nx = 200
ny = 20
Re = 100
sim.parede = int(ny/3)
sim.fator = 5 # n - 3
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

sim.condicoes_contorno_velocidades_duto = sim.condicoes_contorno_velocidades_bif
sim.condicoes_contorno_pressao_duto = sim.condicoes_contorno_pressao_bif

#Executa a simulação
X, Y, u, v, p, velocidade_modulo, tempo = sim.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

# Tempo
print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

sim.plotar_contorno(X, Y, velocidade_modulo, Re, t_final, 'Módulo da Velocidade')
sim.plotar_streamlines(X, Y, u, v, velocidade_modulo, Re, t_final)
