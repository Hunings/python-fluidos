import simulacao_ressalto

simulacao_ressalto.comprimento = 10
simulacao_ressalto.altura = 1
simulacao_ressalto.nx = 50
simulacao_ressalto.ny = 75
simulacao_ressalto.Re = 10

simulacao_ressalto.u_max = 5
simulacao_ressalto.v_max = 5
simulacao_ressalto.tau = 1
simulacao_ressalto.passos_tempo = 10000

simulacao_ressalto.it_pressao = 100
simulacao_ressalto.plotar_a_cada = 1

simulacao_ressalto.comp_bfs = 10    
simulacao_ressalto.alt_bfs = 30

simulacao_ressalto.condicoes_contorno_velocidades_duto = simulacao_ressalto.condicoes_contorno_velocidades_bfs
simulacao_ressalto.condicoes_contorno_pressao_duto = simulacao_ressalto.condicoes_contorno_pressao_bfs

simulacao_ressalto.simulacao(1, 0, 0)