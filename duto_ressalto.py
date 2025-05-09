import simulacao_ressalto
import matplotlib.pyplot as plt

simulacao_ressalto.comprimento = 10
simulacao_ressalto.altura = 1
simulacao_ressalto.nx = 29
simulacao_ressalto.ny = 29
simulacao_ressalto.Re = 100

simulacao_ressalto.u_max = 5
simulacao_ressalto.v_max = 5
simulacao_ressalto.tau = 0.1
simulacao_ressalto.passos_tempo = 3000

simulacao_ressalto.it_pressao = 100
simulacao_ressalto.plotar_a_cada = 1

simulacao_ressalto.comp_bfs = 15    
simulacao_ressalto.alt_bfs = 15

simulacao_ressalto.condicoes_contorno_velocidades_duto = simulacao_ressalto.condicoes_contorno_velocidades_bfs
simulacao_ressalto.condicoes_contorno_pressao_duto = simulacao_ressalto.condicoes_contorno_pressao_bfs

X, Y, u, v, p, velocidade_modulo = simulacao_ressalto.simulacao(1, 0, 0)

plt.figure(figsize=(100, 10))
plt.quiver(X, Y, u, v, velocidade_modulo, scale=30)
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