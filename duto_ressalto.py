import matplotlib.pyplot as plt
import simulacao_ressalto as sim
import numpy as np

sim.comprimento = 10
sim.altura = 1
sim.nx = 200
sim.ny = 200
sim.Re = 1000
sim.dx = sim.comprimento/(sim.nx-1)
sim.dy = sim.altura/(sim.ny-1)
u_max = 5
v_max = 5
tau = 0.1
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 10
sim.passos_tempo = int(t_final/sim.dt)

sim.it_pressao = 100
sim.tol = 1e-1
sim.plotar_a_cada = 1
sim.sx = int(sim.nx/2)
sim.sy = int(sim.ny/2)

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
plt.contourf(X, Y, velocidade_modulo)
plt.colorbar()
plt.show()
plt.contourf(X, Y, w)
plt.colorbar()
plt.show()