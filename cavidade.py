import simulacao as sim
import matplotlib.pyplot as plt
from functools import partial

sim.comprimento = 1
sim.altura = 1
sim.nx = 100
sim.ny = 100
sim.Re = 3000
sim.dx = sim.comprimento / (sim.nx - 1)
sim.dy = sim.altura / (sim.ny - 1)

u_max = 3
v_max = 3
tau = 0.1
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
sim.t_final = 200

sim.tol = 1e-2
sim.passos_tempo = int(sim.t_final/sim.dt)
sim.it_pressao = 100
sim.plotar_a_cada = 100

sim.condicoes_contorno_velocidades_duto = partial(sim.condicoes_contorno_velocidades_cavidade)
sim.condicoes_contorno_pressao_duto = partial(sim.condicoes_contorno_pressao_cavidade)

X, Y, u, v, p, velocidade_modulo = sim.simulacao(0, 0, 0)

#Visualização 

plt.figure(figsize=(11, 10))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=20)
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {sim.t_final}")
plt.colorbar()
plt.show()
plt.figure(figsize=(11, 10))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='jet')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')  
plt.title(f"Re = {sim.Re} t = {sim.t_final}")
plt.show()
plt.figure(figsize=(11, 10))
plt.contourf(X, Y, p)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()