import simulacao
import matplotlib.pyplot as plt

simulacao.comprimento = 10
simulacao.altura = 1
simulacao.nx = 50
simulacao.ny = 50
simulacao.Re = 10

simulacao.u_max = 3
simulacao.v_max = 3
simulacao.tau = 0.1
simulacao.passos_tempo = 100
simulacao.it_pressao = 100
simulacao.plotar_a_cada = 10

simulacao.condicoes_contorno_velocidades_duto = simulacao.condicoes_contorno_velocidades_duto
simulacao.condicoes_contorno_pressao_duto = simulacao.condicoes_contorno_pressao_duto

print(simulacao.dt)
X, Y, u, v, p, velocidade_modulo = simulacao.simulacao(1, 0, 1)

plt.figure(figsize=(12, 6))
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2], color='red')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, color=velocidade_modulo, cmap='viridis')
plt.colorbar()  
plt.show()