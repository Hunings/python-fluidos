import simulação

simulação.comprimento = 5
simulação.altura = 5
simulação.nx = 100
simulação.ny = 100
simulação.Re = 1

simulação.u_max = 2
simulação.v_max = 2
simulação.tau = 0.1
simulação.passos_tempo = 5000
simulação.it_pressao = 100
simulação.plotar_a_cada = 1

simulação.condicoes_contorno_velocidades_duto = simulação.condicoes_contorno_velocidades_cavidade
simulação.condicoes_contorno_pressao_duto = simulação.condicoes_contorno_pressao_cavidade

simulação.simulacao(0, 0, 0)

#CC utilizada foi Neumann nas bordas