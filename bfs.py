import matplotlib.pyplot as plt
import simulacao as sim
import numpy as np


sim.comprimento = 35
sim.altura = 2
sim.nx = 201
sim.ny = 21
sim.Re = 100
sim.dx = sim.comprimento/(sim.nx-1)
sim.dy = sim.altura/(sim.ny-1)
sim.u_max = 5
sim.v_max = 5
tau = 0.1
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/sim.u_max, sim.dy/sim.u_max)

t_final = 10
sim.passos_tempo = int(t_final/sim.dt)

sim.it_pressao = 300
sim.tol = 1e-2
sim.plotar_a_cada = 10
sim.sx = int((sim.nx-1)/6)
sim.sy = int((sim.ny-1)/2)  

sim.condicoes_contorno_pressao_duto = sim.condicoes_contorno_pressao_bfs
sim.condicoes_contorno_velocidades_duto = sim.condicoes_contorno_velocidades_bfs

X, Y, u, v, p, velocidade_modulo, tempo = sim.simulacao(0, 0, 0)

print(f"Tempo de execução: {(tempo):2f} segundos")

dvdx, dudy = np.zeros_like(u), np.zeros_like(v)
dvdx[1:-1, 1:-1] = (v[2:, 1:-1] - v[:-2, 1:-1] ) / (2*sim.dx)
dvdx[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2] ) / (2*sim.dy)
w = 1/2 * (dvdx - dudy)

plt.figure(figsize=(50, 5))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=30)
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {t_final}")
plt.colorbar()
plt.show()
plt.figure(figsize=(50, 5))
plt.streamplot(np.transpose(X), np.transpose(Y), np.transpose(u), np.transpose(v), color=np.transpose(velocidade_modulo), cmap='jet')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {t_final}")
plt.show()
plt.contourf(X, Y, u, levels=200, cmap='jet')
plt.title(f"Re = {sim.Re} t = {t_final}")
plt.colorbar(orientation='horizontal')
plt.show()
plt.contourf(X, Y, v, levels=200, cmap='jet')