import matplotlib.pyplot as plt
import simulacao_ressalto as sim
import numpy as np

sim.comprimento = 10
sim.altura = 1
sim.nx = 40
sim.ny = 80
sim.Re = 4000
sim.dx = sim.comprimento/(sim.nx-1)
sim.dy = sim.altura/(sim.ny-1)
u_max = 100
v_max = 100
tau = 0.1
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 10
sim.passos_tempo = int(t_final/sim.dt)

sim.it_pressao = 100
sim.tol = 1e-2
sim.plotar_a_cada = 1

sim.sx = int(sim.nx/4)
sim.sy = int(sim.ny/4)

X, Y, u, v, p, velocidade_modulo = sim.simulacao(1, 0, 0)
print(np.shape(velocidade_modulo))
plt.figure(figsize=(50, 5))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=30)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.show()
plt.figure(figsize=(50, 5))
plt.streamplot(np.transpose(X), np.transpose(Y), np.transpose(u), np.transpose(v), color=np.transpose(velocidade_modulo), cmap='viridis')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')  
plt.show()
plt.figure(figsize=(50, 5))
plt.contourf(X, Y, p)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()