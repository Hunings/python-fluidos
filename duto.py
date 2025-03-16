import simulação
import matplotlib.pyplot as plt

simulação.comprimento = 10
simulação.altura = 1
simulação.nx = 40
simulação.ny = 40
simulação.Re = 1

simulação.u_max = 1.5
simulação.v_max = 1.5
simulação.tau = 0.1
simulação.passos_tempo = 5000
simulação.it_pressao = 50
simulação.plotar_a_cada = 10
simulação.dt = 0.01

simulação.condicoes_contorno_velocidades_duto = simulação.condicoes_contorno_velocidades_duto
simulação.condicoes_contorno_pressao_duto = simulação.condicoes_contorno_pressao_duto

print(simulação.dt)
simulação.simulacao(1, 0, 1)

