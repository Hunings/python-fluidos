import simulacao_staggered as sim
import matplotlib.pyplot as plt

#Define as constantes da simulação
sim.comprimento = 20
sim.altura = 1
sim.nx = 200
sim.ny = 20
sim.Re = 100
sim.dx = sim.comprimento / (sim.nx -1)
sim.dy = sim.altura / (sim.ny -1)

u_max = 2
v_max = 2
tau = 0.9
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 100
sim.passos_tempo = int(t_final / sim.dt)
sim.it_pressao = 1000
sim.plotar_a_cada = 10

#Executa a simulação
X, Y, u, v, velocidade_modulo, p, tempo = sim.simulacao(1, 0, 0)

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

plt.figure(figsize=(15, 2))
plt.contourf(X, Y, velocidade_modulo, cmap='jet', levels=200)
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
plt.figure(figsize=(15, 2))
plt.contourf(X, Y, p, levels=200, cmap='jet')
plt.colorbar(orientation='horizontal')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()