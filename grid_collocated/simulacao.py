import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
'''
Este código realiza uma simulação 2D de um fluido utilizando diferenças finitas centradas,
e um solver para a equação de Poisson para a pressão que aplica o método iterativo de Jacobi.

Abaixo estão as funções que aplicam as condições de contorno, criam a malha, resolvem as equações e rodam o loop principal.
Para rodar o código é necessário um código auxiliar que inicializa as variáveis e executa algumas das funções. Exemplos destes
códigos auxiliares estão junto com este no GitHub.

A simulação principal consiste basicamente na função simulacao(), que inicia os arrays e o loop para a evolução temporal, e executa as funções para
aplicar as condições de contorno para velocidade e pressão.
'''

abertura = 0
bfs_x = 0
bfs_y = 0

def salvar(u, v, p, it):
    query = input('Quer salvar? [Digite 1 para sim]')
    if query:
        nome = input('Nome do arquivo: ')
        np.savez(rf"E:\Documentos\IC\Códigos\python-fluidos\grid_collocated\dados\{nome}", u=u, v=v, p=p, iteracao=it)
    return
def carregar(u, v, p):
    nome = input('Nome: ')
    dados = np.load(rf"E:\Documentos\IC\Códigos\python-fluidos\grid_collocated\dados\{nome}.npz")
    u[:] = dados["u"]
    v[:] = dados["v"]
    p[:] = dados["p"]
    iteracao = int(dados["iteracao"])
    return u, v, p, iteracao
def estima_x1r(u, X, bfs_x, dx):
    ivort=bfs_x
    for i in range(bfs_x+5,4*bfs_x-1):
        if u[i,1]<0 and u[i+1,1]>0:
            ivort=i
            break
    return X[ivort,1] - bfs_x*dx
def condicoes_contorno_pressao_bfs(p):
    #Paredes Neumann homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    #Paredes Duto Neumann homogênea
    p[:bfs_x+1, bfs_y] = p[:bfs_x+1, bfs_y+1] #Parede norte do ressalto
    p[bfs_x, :bfs_y+1] = p[bfs_x+1, :bfs_y+1]   #Parede leste do ressalto
    return p
def condicoes_contorno_velocidades_bfs(u, v):
    #Paredes No-Slip
    u[:, 0] = 0. # Sul
    v[:, 0] = 0.
    u[:, -1] = 0. # Norte
    v[:, -1] = 0.
    #Entrada Dirichlet
    Umax = 1.
    y = np.linspace(0, 2, 100-bfs_y) # aqui ainda tenho que colocar altura e ny manualmente
    u[0, bfs_y:] = Umax*(1-((y - 1)**2))
    #Saída Neumann homogênea
    u[-1, 1:-1] = u[-2, 1:-1] # Leste
    v[-1, 1:-1] = v[-2, 1:-1]
    #Velocidade no interior da borda
    u[:bfs_x+1, :bfs_y+1] = 0.
    v[:bfs_x+1, :bfs_y+1] = 0.
    return u, v
def condicoes_contorno_velocidades_bif(u, v):
    #Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0
    u[:, -1] = 0 #Norte
    v[:, -1] = 0
    #Entrada Dirichlet
    u[0, abertura:2*abertura] = 1.
    #Saída Neumann homogênea
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    return u, v
def condicoes_contorno_pressao_bif(p):
    #Paredes Neumman homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    #Entrada Neumann homogênea
    p[0, abertura:2*abertura] = p[1, abertura:2*abertura]
    return p
def condicoes_contorno_pressao_duto(p):
    #Paredes Neumann homogênea
    p[-1, :] = p[-2, :] #Leste
    p[:, -1] = p[:, -2] #Norte
    p[:, 0] = p[:, 1] #Sul
    p[0, :] = p[1, :] #Oeste
    return p
def condicoes_contorno_velocidades_duto(u, v):
    #Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0 
    u[:, -1] = 0 #Norte
    v[:, -1] = 0
    #Entrada Dirichlet
    Umax = 1.
    y = np.linspace(0, 2, 31) # aqui ainda tenho que colocar altura e ny manualmente
    u[0, bfs_y:] = Umax*(1-((y - 1)**2))
    u[0, :] = 1.
    #Saída Neumann homogênea
    u[-1, :] = u[-2, :] #Leste
    v[-1, :] = v[-2, :]
    return u, v
def condicoes_contorno_velocidades_cav(u, v):
    # Paredes No-Slip
    u[:, 0] = 0 #Sul
    v[:, 0] = 0
    u[0, :] = 0 #Oeste
    v[0, :] = 0
    u[-1, :] = 0 #Leste
    v[-1, :] = 0
    #Entrada Dirichlet
    u[1:-1, -1] = 1. #Norte
    v[:, -1] = 0.
    return u, v
def pressao(fonte, p, p_novo, dx, dy, it_pressao, tol, residuo):
    j = 0
    normal2 = 1
    denominador = (2 * (dx**2 + dy**2))
    while j < it_pressao and normal2 > tol:
            #Atualiza pressão iterativamente
            p_novo[1:-1, 1:-1] = (
              (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
              +
              (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, 1:-1]))/denominador
            #Condições de Contorno Pressão
            p_novo[:] = condicoes_contorno_pressao_duto(p_novo)
            
            #Calcula medidas do erro
            normal2 = np.linalg.norm(p[1:-1, 1:-1]-p_novo[1:-1, 1:-1], ord=2)
            residuo[1:-1, 1:-1] = ((p_novo[2:,1:-1] - 2*p_novo[1:-1,1:-1] + p_novo[:-2,1:-1]) / dx**2 
                                   + (p_novo[1:-1,2:] - 2*p_novo[1:-1,1:-1] + p_novo[1:-1,:-2]) / dy**2) - fonte[1:-1, 1:-1]
            #Atualiza pressão e avança na iteração
            p[:] = p_novo
            j+=1
    return p, normal2, residuo
def pressao_bfs(fonte, p, p_novo, dx, dy, it_pressao, tol, residuo):
    j = 0
    normal2 = 1
    denominador = (2 * (dx**2 + dy**2))
    while j < it_pressao and normal2 > tol:
            #Atualiza pressão iterativamente
            p_novo[1:-1, bfs_y+1:-1] = (
              (dy**2 * (p[2:, bfs_y+1:-1] + p[:-2, bfs_y+1:-1]))
              +
              (dx**2 * (p[1:-1, bfs_y+2:] + p[1:-1, bfs_y:-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, bfs_y+1:-1]))/denominador
            p_novo[bfs_x+1:-1, 1:bfs_y+1] = (
              (dy**2 * (p[bfs_x+2:, 1:bfs_y+1] + p[bfs_x:-2, 1:bfs_y+1]))
              +
              (dx**2 * (p[bfs_x+1:-1, 2:bfs_y+2] + p[bfs_x+1:-1, :bfs_y]))
              -
              (dx**2 * dy**2 * fonte[bfs_x+1:-1, 1:bfs_y+1]))/denominador
             
            #Condições de Contorno Pressão
            p_novo[:] = condicoes_contorno_pressao_bfs(p_novo)
            
            #Calcula medidas do erro
            normal2 = (np.linalg.norm(p[1:-1, bfs_y+1:-1] - p_novo[1:-1, bfs_y+1:-1], ord=2) 
            + 
            np.linalg.norm(p[bfs_x+1:-1, 1:bfs_y+1] - p_novo[bfs_x+1:-1, 1:bfs_y+1], ord=2))
            residuo[1:-1, bfs_y+1:-1] = ((p_novo[2:, bfs_y+1:-1] - 2*p_novo[1:-1, bfs_y+1:-1] + p_novo[:-2, bfs_y+1:-1]) / dx**2 
                                   + (p_novo[1:-1, bfs_y+2:] - 2*p_novo[1:-1, bfs_y+1:-1] + p_novo[1:-1, bfs_y:-2]) / dy**2) - fonte[1:-1, bfs_y+1:-1]
            residuo[bfs_x+1:-1, 1:bfs_y+1] = ((p_novo[bfs_x+2:, 1:bfs_y+1] - 2*p_novo[bfs_x+1:-1, 1:bfs_y+1] + p_novo[bfs_x:-2, 1:bfs_y+1]) / dx**2 
                                   + (p_novo[bfs_x+1:-1, 2:bfs_y+2] - 2*p_novo[bfs_x+1:-1, 1:bfs_y+1] + p_novo[bfs_x+1:-1, :bfs_y]) / dy**2) - fonte[bfs_x+1:-1, 1:bfs_y+1]
            #Atualiza pressão e avança na iteração
            p[:] = p_novo
            j+=1
    return p, normal2, residuo
def malha(comprimento, altura, nx, ny):
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)
    return X, Y
def plotar_contorno(X, Y, V, Re, t_final, nx, ny, comprimento, altura, x1r, tau, it_pressao, titulo, quadrado):
    if not quadrado:
        c, a = 15, 2
    else:
        c, a = 8, 12
    plt.figure(figsize=(c, a))
    plt.contourf(X, Y, V, levels=200, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(rf"{titulo}, Re = {Re}, t = {t_final}, ({comprimento} x {altura}), ({nx} x {ny}) pontos, X/r = {x1r:.3f}, $\tau$ = {tau}, itp = {it_pressao}")
    plt.colorbar(orientation='horizontal')
    plt.show()
    return
def plotar_streamlines(X, Y, u, v, Re, t_final, nx, ny, comprimento, altura, x1r, tau, it_pressao, titulo, quadrado):
    if not quadrado:
        c, a = 15, 2
    else:
        c, a = 8, 12
    plt.figure(figsize=(c, a))
    plt.streamplot(X.T, Y.T, u.T, v.T, density=0.5, linewidth=0.5, color='k', broken_streamlines=False, arrowstyle='-')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(rf"{titulo}, Re = {Re}, t = {t_final}, ({comprimento} x {altura}), ({nx} x {ny}) pontos, X/r = {x1r:.3f}, $\tau$ = {tau}, itp = {it_pressao}")
    plt.show()
    return
def informacoes(i, passos_tempo, tempo_transcorrido, t_final, V_max, normal2, residuo):
    print('It=', i, '/', passos_tempo, f"t={tempo_transcorrido:.3}/{t_final}", '||V||=', V_max, 'NormaL²=', normal2, '∇·u ', np.linalg.norm(residuo), end=' ')
    return
def evolucao(X, Y, V):
    plt.contourf(X, Y, V, levels=100, cmap='jet')
    plt.colorbar()
    plt.pause(1e-10)
    plt.clf()
    return

def simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, forma='duto', caso='bfs'):
    #Constantes calculadas
    dx = comprimento/(nx-1)
    dy = altura/(ny-1)
    dt = tau*min(Re/2/(1/dx**2 + 1/dy**2), dx/u_max, dy/v_max)
    passos_tempo = int(t_final/dt)
    print(bfs_x, bfs_y)
    print(f"dx = {dx:.3f}, dy = {dy:.3f}, dt = {dt:.3f}, bfs_x = {bfs_x*dx}, bfs_y = {bfs_y*dy}")

    #Conta tempo de simulação
    inicio = perf_counter()

    # Malha
    X, Y = malha(comprimento, altura, nx, ny)

    # Condições Iniciais
    continuar = input("Quer carregar dados? [Digite 1 para sim]")
    if continuar:
        u, v, p = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))
        u, v, p, iteracao_anterior = carregar(u, v, p)
    else:
        if caso == 'bfs':
            u = u0*np.zeros((nx, ny))
            Umax = 1.
            y = np.linspace(0, altura, ny-bfs_y)
            u[0, bfs_y:] = Umax*(1-((y - 1)**2))
            v = v0*np.ones((nx, ny))
            p = p0*np.ones((nx, ny))
        else:
            u = u0*np.ones((nx, ny))
            v = v0*np.ones((nx, ny))
            p = p0*np.ones((nx, ny))
    # Arrays 
    u_, v_ = np.zeros((nx, ny)), np.zeros((nx, ny))
    u_novo, v_novo = np.zeros((nx, ny)), np.zeros((nx, ny))
    difusao_x, difusao_y = np.zeros((nx, ny)), np.zeros((nx, ny))
    conveccao_x, conveccao_y = np.zeros((nx, ny)), np.zeros((nx, ny))
    dpdx, dpdy = np.zeros((nx, ny)), np.zeros((nx, ny))
    fonte = np.zeros((nx, ny))
    residuo = np.zeros((nx, ny))

    u, v = condicoes_contorno_velocidades_duto(u, v)
    p = condicoes_contorno_pressao_duto(p)

    tempo_transcorrido = 0
    plotar_evolucao = bool(input('Plotar evolução temporal? [Digite 1 para sim]'))
    if plotar_evolucao:
        if forma == 'quadrado':
            plt.figure(figsize=(10, 8))
        else:
            plt.figure(figsize=(15,2))
    # Iteração 
    try:
        for it in range(passos_tempo):
            difusao_x[1:-1, 1:-1] = 1/Re * ((u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx**2 +
                                            (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2]) / dy**2)
            
            conveccao_x[1:-1, 1:-1] = (v[1:-1, 2:] * u[1:-1, 2:] - v[1:-1, :-2] * u[1:-1, :-2])/(2*dy) + (
                u[2:, 1:-1]**2 - u[:-2, 1:-1]**2
            )/(2*dx)

            # Velocidade u antes da correção da pressão
            u_[1:-1, 1:-1] = u[1:-1, 1:-1] + dt*(difusao_x[1:-1, 1:-1] - conveccao_x[1:-1, 1:-1])

            difusao_y[1:-1, 1:-1] = 1/Re * ((v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1]) / dx**2 +
                                            (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, :-2]) / dy**2)
            
            conveccao_y[1:-1, 1:-1] = (v[2:, 1:-1]*u[2:, 1:-1] - v[:-2, 1:-1]*u[:-2, 1:-1])/(2*dx) + (
                v[1:-1, 2:]**2 - v[1:-1, :-2]**2
            )/(2*dy)
            
            #Velocidade v antes da correção da pressão
            v_[1:-1, 1:-1] = v[1:-1, 1:-1] + dt*(difusao_y[1:-1, 1:-1] - conveccao_y[1:-1, 1:-1])

            u_[:], v_[:] = condicoes_contorno_velocidades_duto(u_, v_)

            fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt
            p_novo = np.copy(p)
            if caso == 'bfs':
                #fonte[1:-1, 1:-1] = ((u_[2:, 1:-1] - u_[:-2, 1:-1])/(2*dx) + (v_[1:-1, 2:] - v_[1:-1, :-2])/(2*dy))/dt
                num_pontos = (nx-2)*(ny-2) - (bfs_x)*(bfs_y)
                
                int_up = fonte[1:-1, bfs_y+1:-1]
                int_down = fonte[bfs_x+1:-1, 1:bfs_y+1]
                mediaf = (np.sum(int_down)+np.sum(int_up))/num_pontos
                fonte[1:-1, bfs_y+1:-1] += -mediaf
                fonte[bfs_x+1:-1, 1:bfs_y+1] += -mediaf

                p_novo[:], normal2, residuo[:] = pressao_bfs(fonte, p, p_novo, dx, dy, it_pressao, tol, residuo)
                
                dpdx[1:-1, bfs_y+1:-1] = (p_novo[2:, bfs_y+1:-1] - p_novo[:-2, bfs_y+1:-1])/(2*dx)
                dpdx[bfs_x+1:-1, 1:bfs_y+1] = (p_novo[bfs_x+2:, 1:bfs_y+1] - p_novo[bfs_x:-2, 1:bfs_y+1])/(2*dx)
                
                dpdy[1:-1, bfs_y+1:-1] = (p_novo[1:-1, bfs_y+2:] - p_novo[1:-1, bfs_y:-2])/(2*dy)
                dpdy[bfs_x+1:-1, 1:bfs_y+1] = (p_novo[bfs_x+1:-1, 2:bfs_y+2] - p_novo[bfs_x+1:-1, :bfs_y])/(2*dy)

            else:
                num_pontos = (nx-2)*(ny-2)
                mediaf = (np.sum(fonte[1:-1, 1:-1]))/num_pontos
                fonte[1:-1, 1:-1] += -mediaf
                p_novo[:], normal2, residuo[:] = pressao(fonte, p, p_novo, dx, dy, it_pressao, tol, residuo)

                dpdx[1:-1, 1:-1] = (p_novo[2:, 1:-1] - p_novo[:-2, 1:-1])/(2*dx)
                dpdy[1:-1, 1:-1] = (p_novo[1:-1, 2:] - p_novo[1:-1, :-2])/(2*dy)
                

            u_novo[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt
            v_novo[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt

            # Condições de Contorno Velocidades Finais
            u_novo[:], v_novo[:] = condicoes_contorno_velocidades_duto(u_novo, v_novo)

            u[:], v[:], p[:] = u_novo, v_novo, p_novo

            tempo_transcorrido += dt

            V = (u**2 + v**2)**(0.5)
            V_max = np.max(V)

            x1r = estima_x1r(u, X, bfs_x, dx)
            informacoes(it, passos_tempo, tempo_transcorrido, t_final, V_max, normal2, residuo)
            print(f"X/r={x1r}")
            
            if V_max > 2 or np.isnan(V_max):
                break
            if plotar_evolucao and it % plotar_a_cada == 0:
                if caso == 'bfs':
                    plt.plot(X[:bfs_x+1, bfs_y], Y[:bfs_x+1, bfs_y], c='w') # Horizontal
                    plt.plot(X[bfs_x, :bfs_y+1], Y[bfs_x, :bfs_y+1], c='w') #Vertical
                    evolucao(X, Y, V)
                else:
                    evolucao(X, Y, V)
    except KeyboardInterrupt:
        pass
    fim = perf_counter()
    tempo = fim - inicio
    salvar(u, v, p, it)
    return X, Y, u, v, p, V, tempo, x1r