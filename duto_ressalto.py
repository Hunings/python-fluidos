import simulacao

simulacao.comprimento = 10
simulacao.altura = 1
simulacao.nx = 150
simulacao.ny = 200
simulacao.Re = 1

simulacao.u_max = 3
simulacao.v_max = 3
simulacao.tau = 0.8
simulacao.passos_tempo = 1000
simulacao.it_pressao = 1000
simulacao.plotar_a_cada = 1

simulacao.comprimento_ressalto_pontos = 30
simulacao.altura_ressalto_pontos = 30

simulacao.condicoes_contorno_velocidades_duto = simulacao.condicoes_contorno_velocidades_bfs
simulacao.condicoes_contorno_pressao_duto = simulacao.condicoes_contorno_pressao_bfs

simulacao.simulacao(1, 0, 0)