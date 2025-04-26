import simulacao
import matplotlib.pyplot as plt

simulacao.comprimento = 20
simulacao.altura = 1
simulacao.nx = 30
simulacao.ny = 30
simulacao.Re = 1

simulacao.u_max = 3
simulacao.v_max = 3
simulacao.tau = 0.1
simulacao.passos_tempo = 10
simulacao.it_pressao = 100
simulacao.plotar_a_cada = 10

simulacao.condicoes_contorno_velocidades_duto = simulacao.condicoes_contorno_velocidades_duto
simulacao.condicoes_contorno_pressao_duto = simulacao.condicoes_contorno_pressao_duto

print(simulacao.dt)
simulacao.simulacao(1, 0, 1)

