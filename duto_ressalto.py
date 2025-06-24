import matplotlib.pyplot as plt
import simulacao_ressalto as sim
import numpy as np


sim.comprimento = 10
sim.altura = 1
sim.nx = 40
sim.ny = 80
sim.Re = 100
sim.dx = sim.comprimento/(sim.nx-1)
sim.dy = sim.altura/(sim.ny-1)
sim.u_max = 5
sim.v_max = 5
tau = 1 
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/sim.u_max, sim.dy/sim.u_max)

t_final = 1000
sim.passos_tempo = int(t_final/sim.dt)

sim.it_pressao = 100
sim.tol = 1e-2
sim.plotar_a_cada = 10
sim.sx = int(sim.nx/4)
sim.sy = int(sim.ny/4)

X, Y, u, v, p, velocidade_modulo = sim.simulacao(1, 0, 0)

dvdx, dudy = np.zeros_like(u), np.zeros_like(v)
dvdx[1:-1, 1:-1] = (v[2:, 1:-1] - v[:-2, 1:-1] ) / (2*sim.dx)
dvdx[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2] ) / (2*sim.dy)
w = 1/2 * (dvdx - dudy)

print(np.shape(velocidade_modulo))
plt.figure(figsize=(50, 5))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=30)
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {t_final}")
plt.colorbar()
plt.show()
plt.figure(figsize=(50, 5))
plt.streamplot(np.transpose(X), np.transpose(Y), np.transpose(u), np.transpose(v), color=np.transpose(velocidade_modulo), cmap='viridis')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {t_final}")
plt.show()