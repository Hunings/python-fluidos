import simulacao_ressalto_staggered as sim
import matplotlib.pyplot as plt

#Define as constantes da simulação
sim.comprimento = 35
sim.altura = 2
sim.nx = 70
sim.ny = 50
sim.Re = 100
sim.dx = sim.comprimento / (sim.nx -1)
sim.dy = sim.altura / (sim.ny -1)

u_max = 3
v_max = 3
tau = 0.5
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 100
sim.passos_tempo = int(t_final / sim.dt)
sim.it_pressao = 1000
sim.plotar_a_cada = 100

# Modifica as condições de contorno internas

#Executa a simulação
X, Y, u, v, velocidade_modulo, p = sim.simulacao(3, 0, 0)

#Visualização 

plt.figure(figsize=(12, 6))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=40)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='jet')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')  
plt.show()
plt.figure(figsize=(12, 6))
plt.contourf(X, Y, p)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()