import simulacao as sim
import matplotlib.pyplot as plt
from functools import partial

#Define as constantes da simulação
sim.comprimento = 10
sim.altura = 1
sim.nx = 40
sim.ny = 80
sim.Re = 4000
sim.dx = sim.comprimento / (sim.nx -1)
sim.dy = sim.altura / (sim.ny -1)
sim.parede = int(sim.ny/3)
sim.fator = 2 # n - 3
sim.u_max = 10
sim.v_max = 10
sim.tol = 1e-2

sim.dt = sim.tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/sim.u_max, sim.dy/sim.u_max)
sim.t_final = 10
sim.passos_tempo = int(sim.t_final / sim.dt)
sim.it_pressao = 100
sim.plotar_a_cada = 100

# Modifica as condições de contorno internas

sim.condicoes_contorno_velocidades_duto = partial(sim.bifurcacao_velocidades, fator=sim.fator)
sim.condicoes_contorno_pressao_duto = partial(sim.bifurcacao_pressao, fator=sim.fator)

#Executa a simulação
X, Y, u, v, p, velocidade_modulo = sim.simulacao(1, 0, 1)

#Visualização 
print()
plt.figure(figsize=(12, 6))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=40)
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Re = {sim.Re} t = {sim.t_final}")
plt.colorbar()
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='viridis')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')  
plt.title(f"Re = {sim.Re} t = {sim.t_final}")
plt.show()
plt.figure(figsize=(12, 6))
plt.contourf(X, Y, p)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()