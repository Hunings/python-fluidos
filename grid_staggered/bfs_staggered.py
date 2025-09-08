import simulacao_staggered as sim
import matplotlib.pyplot as plt

#Define as constantes da simulação
comprimento = 24
altura = 2
nx = 200
ny = 60
Re = 100
u_max = 2
v_max = 2
tau = 1
t_final = 500
it_pressao = 100
plotar_a_cada = 1
sim.bfs_x = int(nx/6)
sim.bfs_y = int(ny/2)
u0 = 1
p0 = 0
v0 = 0
tol = 1e-3

sim.condicoes_contorno_velocidades = sim.condicoes_contorno_velocidades_bfs
sim.condicoes_contorno_pressao = sim.condicoes_contorno_pressao_bfs

#Executa a simulação
X, Y, u, v, velocidade_modulo, p, tempo = sim.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, forma='duto')

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

sim.plotar_contorno(X, Y, velocidade_modulo, Re, t_final, nx, ny, comprimento, altura, tau, it_pressao, 'Módulo da Velocidade (staggered)', False)
sim.plotar_contorno(X, Y, p, Re, t_final, nx, ny, comprimento, altura, tau, it_pressao, 'Pressão (staggered)', False)
sim.plotar_vetores(X, Y, u, v, velocidade_modulo, Re, t_final, 50, 2, False)