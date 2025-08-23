import simulacao_staggered as sim
import matplotlib.pyplot as plt

#Define as constantes da simulação
sim.comprimento = 24
sim.altura = 2
sim.nx = 201
sim.ny = 41
sim.Re = 100
u_max = 2
v_max = 2
tau = 0.5
sim.t_final = 10
sim.it_pressao = 200
sim.plotar_a_cada = 1
sim.bfs_x = int(sim.nx/3)
sim.bfs_y = int(sim.ny/3)
tol = 1e-3

sim.condicoes_contorno_velocidades = sim.condicoes_contorno_velocidades_bfs
sim.condicoes_contorno_pressao = sim.condicoes_contorno_pressao_bfs

#Executa a simulação
X, Y, u, v, velocidade_modulo, p, tempo = sim.simulacao(1, 0, 0, tau)

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

plt.figure(figsize=(15, 2))
plt.contourf(X, Y, u, cmap='jet', levels=200)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(orientation='horizontal')
plt.show()
plt.figure(figsize=(15, 2))
plt.contourf(X, Y, v, cmap='jet', levels=200)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(orientation='horizontal')
plt.show()
plt.figure(figsize=(15, 2))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='jet')
plt.colorbar(orientation='horizontal')
plt.xlabel('X')
plt.ylabel('Y')  
plt.show()