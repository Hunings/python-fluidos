import simulacao
import matplotlib.pyplot as plt

#Define as constantes da simulação
simulacao.comprimento = 10
simulacao.altura = 1
simulacao.nx = 30
simulacao.ny = 30
simulacao.Re = 100

simulacao.u_max = 3
simulacao.v_max = 3
simulacao.tau = 0.1
simulacao.passos_tempo = 10000
simulacao.it_pressao = 100
simulacao.plotar_a_cada = 10

# Modifica as condições de contorno internas

simulacao.condicoes_contorno_velocidades_duto = simulacao.condicoes_contorno_velocidades_duto
simulacao.condicoes_contorno_pressao_duto = simulacao.condicoes_contorno_pressao_duto

print(simulacao.dt)

#Executa a simulação
X, Y, u, v, p, velocidade_modulo = simulacao.simulacao(1, 0, 1)

#Visualização 

plt.figure(figsize=(12, 6))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=40)
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo, cmap='viridis')
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