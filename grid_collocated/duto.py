import simulacao as s
import matplotlib.pyplot as plt

#Define as constantes da simulação
comprimento = 20
altura = 1
nx = 201
ny = 31
Re = 10
u_max = 2
v_max = 2
tau = 0.1
t_final = 6
tol = 1e-2
it_pressao = 200
plotar_a_cada = 100
u0, v0, p0 = 1, 0, 0

#Condições de contorno internas padrão

#Executa a simulação
X, Y, u, v, p, V, tempo, it = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0)

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

s.plotar_contorno(X, Y, V, Re, t_final, 'Módulo da Velocidade', 0)
s.plotar_streamlines(X, Y, u, v, V, Re, t_final, 0)
s.plotar_vetores(X, Y, u, v, V, Re, t_final, 60, 2, 0)
s.plotar_contorno(X, Y, v, Re, t_final, 'v', 0)
s.plotar_contorno(X, Y, p, Re, t_final, 'Pressão', 0)

