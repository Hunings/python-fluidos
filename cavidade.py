import simulacao as sim
import matplotlib.pyplot as plt

sim.comprimento = 1
sim.altura = 1
sim.nx = 80
sim.ny = 80
sim.Re = 3000

u_max = 3
v_max = 3
tau = 0.1
sim.dt = tau*min(sim.Re/2*(1/sim.dx**2 + 1/sim.dy**2), sim.dx/u_max, sim.dy/u_max)
t_final = 200

sim.passos_tempo = int(t_final/sim.dt)
sim.it_pressao = 100
sim.plotar_a_cada = 10

sim.condicoes_contorno_velocidades_duto = sim.condicoes_contorno_velocidades_cavidade
sim.condicoes_contorno_pressao_duto = sim.condicoes_contorno_pressao_cavidade

X, Y, u, v, p, velocidade_modulo = sim.simulacao(0, 0, 1)

#Visualização 

plt.figure(figsize=(11, 10))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=20)
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.show()
plt.figure(figsize=(11, 10))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo.T, cmap='viridis')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')  
plt.show()
plt.figure(figsize=(11, 10))
plt.contourf(X, Y, p)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
#CC utilizada foi Neumann nas bordas