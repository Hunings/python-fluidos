import simulação

simulação.comprimento = 10
simulação.altura = 1
simulação.nx = 100
simulação.ny = 100
simulação.Re = 10

simulação.u_max = 3
simulação.v_max = 3
simulação.tau = 0.1
simulação.passos_tempo = 1000
simulação.it_pressao = 1000
simulação.plotar_a_cada = 1

simulação.proporcao_ressalto_duto_comprimento = 1
simulação.proporcao_ressalto_duto_altura = 1

simulação.dt = 1e-6


simulação.condicoes_contorno_velocidades_duto = simulação.condicoes_contorno_velocidades_bfs
simulação.condicoes_contorno_pressao_duto = simulação.condicoes_contorno_pressao_bfs

simulação.simulacao(1, 0, 0)