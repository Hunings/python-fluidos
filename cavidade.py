import simulacao
import matplotlib.pyplot as plt

simulacao.comprimento = 1
simulacao.altura = 1
simulacao.nx = 30
simulacao.ny = 30
simulacao.Re = 10

simulacao.u_max = 5
simulacao.v_max = 5
simulacao.tau = 0.1
simulacao.passos_tempo = 1000
simulacao.it_pressao = 100
simulacao.plotar_a_cada = 1

simulacao.condicoes_contorno_velocidades_duto = simulacao.condicoes_contorno_velocidades_cavidade
simulacao.condicoes_contorno_pressao_duto = simulacao.condicoes_contorno_pressao_cavidade

X, Y, u, v, p, velocidade_modulo = simulacao.simulacao(-1, 0, 1)

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