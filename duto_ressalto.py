import simulacao_ressalto
import matplotlib.pyplot as plt

simulacao_ressalto.comprimento = 10
simulacao_ressalto.altura = 1
simulacao_ressalto.nx = 50
simulacao_ressalto.ny = 50
simulacao_ressalto.Re = 10

simulacao_ressalto.u_max = 5
simulacao_ressalto.v_max = 5
simulacao_ressalto.tau = 0.1
simulacao_ressalto.passos_tempo = 100000

simulacao_ressalto.it_pressao = 100
simulacao_ressalto.plotar_a_cada = 1

simulacao_ressalto.comp_bfs = 15    
simulacao_ressalto.alt_bfs = 30

simulacao_ressalto.condicoes_contorno_velocidades_duto = simulacao_ressalto.condicoes_contorno_velocidades_bfs
simulacao_ressalto.condicoes_contorno_pressao_duto = simulacao_ressalto.condicoes_contorno_pressao_bfs

X, Y, u, v, p, velocidade_modulo = simulacao_ressalto.simulacao(1, 0, 0)

plt.figure(figsize=(12, 6))
plt.contourf(X, Y, velocidade_modulo, levels=10)
plt.plot(X[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], Y[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], c='black')
plt.plot(X[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], Y[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], c='black')
plt.quiver(X[::2, ::1], Y[::2, ::1], u[::2, ::1], v[::2, ::1], color='white')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.figure(figsize=(12, 6))
plt.pcolormesh(X, Y, velocidade_modulo, cmap='viridis')
plt.plot(X[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], Y[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], c='black')
plt.plot(X[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], Y[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], c='black')
#plt.quiver(X[::1, ::4], Y[::1, ::4], u[::1, ::4], v[::1, ::4], color='white')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.figure(figsize=(12, 6))
plt.streamplot(X.T, Y.T, u.T, v.T, cmap='viridis', density=2)
plt.plot(X[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], Y[:simulacao_ressalto.comp_bfs+1, simulacao_ressalto.alt_bfs], c='black')
plt.plot(X[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], Y[simulacao_ressalto.comp_bfs, :simulacao_ressalto.alt_bfs+1], c='black')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()