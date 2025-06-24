import simulacao as sim
import matplotlib.pyplot as plt

#Define as constantes da simulação
sim.comprimento = 10
sim.altura = 1
sim.nx = 30
sim.ny = 60
sim.Re = 1000
sim.dx = sim.comprimento / (sim.nx -1)
sim.dy = sim.altura / (sim.ny -1)

u_max = 3
v_max = 3
tau = 0.5
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 100
sim.passos_tempo = int(t_final / sim.dt)
sim.it_pressao = 100
sim.plotar_a_cada = 100

# Modifica as condições de contorno internas

sim.condicoes_contorno_velocidades_duto = sim.condicoes_contorno_velocidades_duto
sim.condicoes_contorno_pressao_duto = sim.condicoes_contorno_pressao_duto

#Executa a simulação
X, Y, u, v, p, velocidade_modulo = sim.simulacao(1, 0, 1)

#Visualização 

plt.figure(figsize=(12, 6))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=40)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='viridis')
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