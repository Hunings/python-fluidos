import simulacao as sim
import matplotlib.pyplot as plt
from functools import partial

#Define as constantes da simulação
sim.comprimento = 10
sim.altura = 1
sim.nx = 80
sim.ny = 120
sim.Re = 2000
sim.dx = sim.comprimento / (sim.nx -1)
sim.dy = sim.altura / (sim.ny -1)
sim.parede = int(sim.ny/7)
sim.fator = 4
sim.u_max = 10
sim.v_max = 10

sim.dt = sim.tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/sim.u_max, sim.dy/sim.u_max)
t_final = 10
sim.passos_tempo = int(t_final / sim.dt)
sim.it_pressao = 100
sim.plotar_a_cada = 100

# Modifica as condições de contorno internas

sim.condicoes_contorno_velocidades_duto = partial(sim.bifurcacao_velocidades, fator=sim.fator)
sim.condicoes_contorno_pressao_duto = partial(sim.bifurcacao_pressao, fator=sim.fator)

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