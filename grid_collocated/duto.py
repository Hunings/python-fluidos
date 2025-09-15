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
t_final = 10
tol = 1e-5
it_pressao = 200
plotar_a_cada = 10
u0, v0, p0 = 1, 0, 0

#Condições de contorno internas padrão

#Executa a simulação
X, Y, u, v, p, V, tempo, x1r = s.simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, caso='duto')

print(f"Tempo de execução: {(tempo):2f} segundos")

#Visualização 

s.plotar_contorno(X, Y, V, Re, t_final,  nx, ny, comprimento, altura, x1r, tau, it_pressao, 'Módulo da Velocidade', False)
s.plotar_streamlines(X, Y, u, v, Re, t_final, nx, ny, comprimento, altura, x1r, tau, it_pressao, 'streamlines', False)

