import matplotlib.pyplot as plt
import numpy as np
from time import perf_counter

bfs_x = 0
bfs_y = 0 

def plotar_contorno(X, Y, V, Re, t_final, nx, ny, comprimento, altura, tau, it_pressao, titulo, quadrado):
    if not quadrado:
        c, a = 15, 2
    else:
        c, a = 8, 12
    plt.figure(figsize=(c, a))
    plt.contourf(X, Y, V, levels=200, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"{titulo}, Re = {Re}, t = {t_final}, ({comprimento} x {altura}), ({nx} x {ny}) pontos")
    plt.colorbar(orientation='horizontal')
    plt.show()
    return
def plotar_streamlines(X, Y, u, v, V, Re, t_final, quadrado):
    if not quadrado:
        c, a = 15, 2
    else:
        c, a = 8, 12
    plt.figure(figsize=(c, a))
    plt.streamplot(X.T, Y.T, u.T, v.T, color=V.T, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"Re = {Re}, t = {t_final}")
    plt.colorbar(orientation='horizontal')
    plt.show()
    return
def plotar_vetores(X, Y, u, v, V, Re, t_final, escala, step, quadrado):
    if not quadrado:
        c, a = 15, 2
    else:
        c, a = 8, 12
    plt.figure(figsize=(c, a))
    plt.quiver(X[::step, ::step], Y[::step, ::step], u[::step, ::step], v[::step, ::step], V[::step, ::step], scale=escala, cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"Re = {Re}, t = {t_final}")
    plt.colorbar(orientation='horizontal')
    plt.show()  
    return
def informacoes(i, passos_tempo, tempo_transcorrido, t_final, u_max, v_max, normal2):
    print('It=', i, '/', passos_tempo, f"t={tempo_transcorrido:.3}/{t_final}", '||u||=', u_max, '||v||=', v_max, 'Resíduo=', normal2, end=' ')
    return
def estima_x1r(u,X, bfs_x, dx):
    ivort=bfs_x
    for i in range(bfs_x, 5*bfs_x):
        if u[i,1]<0 and u[i+1,1]>0:
            ivort=i
            break
    return (X[ivort,1]-bfs_x*dx)
def salvar(u, v, p):
    query = input('Quer salvar? [Digite 1 para sim]')
    if query:
        nome = input('Nome do arquivo: ')
        np.savez(rf"E:\Documentos\IC\Códigos\python-fluidos\grid_staggered\dados\{nome}", u=u, v=v, p=p)
    return
def carregar(u, v, p):
    nome = input('Nome: ')
    dados = np.load(rf"E:\Documentos\IC\Códigos\python-fluidos\grid_staggered\dados\{nome}.npz")
    u[:] = dados["u"]
    v[:] = dados["v"]
    p[:] = dados["p"]
    return u, v, p
def condicoes_contorno_velocidades(u, v):
    #Sul
    u[:, 0] = -u[:, 1] # No-Slip
    v[:, 0] = 0 # No-Slip

    #Norte
    u[:, -1] = -u[:, -2] # No-slip
    v[:, -1] = 0 # No-slip

    #Leste
    v[-1, :] = v[-2, :] # Neumann
    u[-1, :] = u[-2, :] # Neumann

    #Oeste
    v[0, :] = -v[1, :] # Dirichlet = 0
    u[0, :] = 1.0 # Dirichlet

    return u, v
def condicoes_contorno_velocidades_cavidade(u, v):
    #Sul
    u[:, 0] = -u[:, 1] # No-Slip
    v[:, 0] = 0 # No-Slip
    #Norte
    u[:, -1] = 2 - u[:, -2] # Dirichlet = 1
    v[:, -1] = 0 # Dirichlet = 0
    #Leste
    v[-1, :-1] = -v[-2, :-1] # No-slip
    u[-1, :-1] = 0 # No-slip
    #Oeste
    v[0, :-1] = -v[1, :-1] # No-slip
    u[0, :-1] = 0 # No-slip
    return u, v
def condicoes_contorno_velocidades_bfs(u, v):
    #Condições de entrada
    u[0, bfs_y+1:-1] = 1.0 #Oeste
    v[0, bfs_y+1:-1] = -v[1, bfs_y+1:-1]
    #Condições na parede
    u[:, -1] = -u[:, -2] #Norte
    v[:, -1] = 0
    u[bfs_x:, 0] = -u[bfs_x+1, 1] #Sul
    v[bfs_x:, 0] = 0
    #Saída
    u[-1, 1:-1] = u[-2, 1:-1]
    v[-1, 1:-1] = v[-2, 1:-1]
    #Ressalto
    u[:bfs_x+1, bfs_y] = -u[:bfs_x+1, bfs_y+1] #Parte de cima do ressalto
    v[:bfs_x+1, bfs_y] = 0.
    u[bfs_x, :bfs_y+1] = 0. #Parte lateral do ressalto
    v[bfs_x, :bfs_y+1] = -v[bfs_x+1, :bfs_y+1]
    #Interior do ressalto
    u[:bfs_x+1, :bfs_y+1] = 0.  
    v[:bfs_x+1, :bfs_y+1] = 0.
    return u, v
def condicoes_contorno_pressao_bfs(p):
    p[bfs_x+1:, 0] = p[bfs_x+1:, 1] #Sul
    p[1:-1, -1] = p[1:-1, -2] #Norte
    p[-1, 1:-1] = p[-2, 1:-1] #Leste
    p[0, bfs_y+1:-1] = p[1, bfs_y+1:-1] #Oeste
    p[:, 0] = p[:, 1] #Sul
    p[:, -1] = p[:, -2] #Norte
    p[-1, :] = 0 #Leste
    p[0, :] = p[1, :] #Oeste
    #Ressalto
    p[1:bfs_x+1, bfs_y] = p[1:bfs_x+1, bfs_y+1]
    p[bfs_x, 1:bfs_y+1] = p[bfs_x+1, 1:bfs_y+1]
    return p
def condicoes_contorno_pressao(p):
    p[:, 0] = p[:, 1] #Sul
    p[:, -1] = p[:, -2] #Norte
    p[-1, :] = p[-2, :] #Leste
    p[0, :] = p[1, :] #Oeste
    return p
def pressao(fonte, p, p_novo, dx, dy, it_pressao, tol):
    j = 0
    normal2 = 1
    while j < it_pressao and normal2 > tol:
            #Atualiza pressão iterativamente
            p_novo[1:-1, 1:-1] = (
              (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]))
              +
              (dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]))
              -
              (dx**2 * dy**2 * fonte[1:-1, 1:-1]))/(2 * (dx**2 + dy**2))
            #Condições de Contorno Pressão
            p_novo[:] = condicoes_contorno_pressao(p_novo)
            #Calcula medidas do erro
            normal2 = np.linalg.norm(p-p_novo, ord=2)
            #Atualiza pressão e avança na iteração
            p[:] = p_novo
            j+=1
    return p, normal2
def pressao_bfs(fonte, p, p_novo, dx, dy, it_pressao, tol):
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
            p_novo[:] = condicoes_contorno_pressao(p_novo)
            
            #Calcula medidas do erro
            normal2 = (np.linalg.norm(p[1:-1, bfs_y+1:-1] - p_novo[1:-1, bfs_y+1:-1], ord=np.inf) 
            + 
            np.linalg.norm(p[bfs_x+1:-1, 1:bfs_y+1] - p_novo[bfs_x+1:-1, 1:bfs_y+1], ord=np.inf))

            #Atualiza pressão e avança na iteração
            p[:] = p_novo
            j+=1
    return p, normal2
def malha(comprimento, altura, nx, ny):
    x = np.linspace(0.0, comprimento, nx)
    y = np.linspace(0.0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X, Y = np.transpose(X), np.transpose(Y)
    return X, Y
def simulacao(comprimento, altura, nx, ny, Re, tol, u_max, v_max, tau, t_final, it_pressao, plotar_a_cada, u0, v0, p0, forma='duto', caso='bfs'):
    dx = comprimento/(nx-1)
    dy = altura/(ny-1)
    dt = tau*min(Re/2/(1/dx**2 + 1/dy**2), dx/u_max, dy/v_max)
    passos_tempo = int(t_final/dt)

    print(f"dx = {dx:.3f}, dy = {dy:.3f}, dt = {dt:.3f}")

    inicio = perf_counter()
    X, Y = malha(comprimento, altura, nx, ny)
    continuar = input("Quer carregar dados? [Digite 1 para sim]")
    if continuar:
        u, v, p = np.zeros((nx, ny+1)), np.zeros((nx+1, ny)), np.zeros((nx+1, ny+1))
        u, v, p = carregar(u, v, p)
    else:
        u = np.zeros((nx, ny+1))
        u[:, bfs_y:] = u0
        v = v0*np.ones((nx+1, ny))
        p = p0*np.ones((nx+1, ny+1))

    u, v = condicoes_contorno_velocidades(u, v)
    p = condicoes_contorno_pressao(p)

    u_, v_ = np.zeros_like(u), np.zeros_like(v)
    u_novo, v_novo = np.zeros_like(u), np.zeros_like(v)
    fonte = np.zeros_like(p)
    dpdx, dpdy = np.zeros_like(u), np.zeros_like(v)

    t=0
    plotar_evolucao = bool(input('Plotar evolução temporal? [Digite 1 para sim]'))
    if plotar_evolucao:
        if forma == 'quadrado':
            plt.figure(figsize=(10, 8))
        else:
            plt.figure(figsize=(15,2))
    try:
        for it in range(passos_tempo):
            u_maximo = np.max(u[1:-1, 1:-1])
            v_maximo = np.max(v[1:-1, 1:-1])
            gamma = 0#max(u_maximo*dt/dx, v_maximo*dt/dy)
            difusao_x = 1/Re * ( # 5x4 nesse caso
                (u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/dx**2
                +
                (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2])/dy**2
            )
            difusao_y = 1/Re * ( # 6x3 nesse caso
                (v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1])/dx**2
                +
                (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, :-2])/dy**2
            )
            du2dx = (((u[1:-1, 1:-1] + u[2:, 1:-1])**2 - (u[1:-1, 1:-1] + u[:-2, 1:-1])**2)/(4*dx) # 5x4
            + gamma*(abs(u[1:-1, 1:-1] + u[2:, 1:-1])*(u[1:-1, 1:-1] - u[2:, 1:-1]) 
                - 
                abs(u[:-2, 1:-1] + u[1:-1, 1:-1])*(u[:-2, 1:-1] - u[1:-1, 1:-1])
                )/(4*dx)
            )
            duvdy = (((v[1:-2, 1:] + v[2:-1, 1:])*(u[1:-1, 1:-1] + u[1:-1, 2:]) - (v[1:-2, :-1] + v[2:-1, :-1])*(u[1:-1, :-2] + u[1:-1, 1:-1]))/(4*dy)
            + gamma*(abs(v[1:-2, 1:] + v[2:-1, 1:])*(u[1:-1, 1:-1] - u[1:-1, 2:]) 
                - 
                abs(v[1:-2, :-1] + v[2:-1, :-1])*(u[1:-1, :-2] - u[1:-1, 1:-1])
                )/(4*dy)
            )
            conveccao_x = du2dx + duvdy
            dv2dy2 = (((v[1:-1, 1:-1] + v[1:-1, 2:])**2 - (v[1:-1, 1:-1] + v[1:-1, :-2])**2)/(4*dy) #6x3
            + gamma*(abs(v[1:-1, 1:-1] + v[1:-1, 2:]*(v[1:-1, 1:-1] - v[1:-1, 2:])) 
                - 
                abs(v[1:-1, :-2] + v[1:-1, 1:-1])*(v[1:-1, :-2] - v[1:-1, 1:-1])
                )/(4*dy)
            )
            duvdx = (((u[1:, 1:-2] + u[:-1, 2:-1])*(v[1:-1, 1:-1] + v[2:, 1:-1]) - (u[:-1, 1:-2] + u[:-1, 2:-1])*(v[:-2, 1:-1] + v[1:-1, 1:-1]))/(4*dx) 
            + gamma*(abs(u[1:, 1:-2] + u[:-1, 2:-1])*(v[1:-1, 1:-1] - v[2:, 1:-1]) 
                -
                abs(u[:-1, 1:-2] + u[:-1, 2:-1])*(v[:-2, 1:-1] - v[1:-1, 1:-1])
                )/(4*dx)
            )
            conveccao_y = dv2dy2 + duvdx

            u_[1:-1, 1:-1] = u[1:-1, 1:-1] + dt*(difusao_x - conveccao_x) #5x4 considerando os [1:-1, 1:-1]
            v_[1:-1, 1:-1] = v[1:-1, 1:-1] + dt*(difusao_y - conveccao_y) #6x3 considerando os [1:-1, 1:-1]

            u_, v_ = condicoes_contorno_velocidades(u_, v_)

            fonte[1:-1, 1:-1] = ((u_[1:, 1:-1] - u_[:-1, 1:-1])/dx + (v_[1:-1, 1:] - v_[1:-1, :-1])/dy)/dt #6x4

            p_novo = np.copy(p)
            if caso == 'bfs':
                p_novo[:], normal2 = pressao_bfs(fonte, p, p_novo, dx, dy, it_pressao, tol)
            else:    
                p_novo[:], normal2 = pressao(fonte, p, p_novo, dx, dy, it_pressao, tol)

            dpdx[1:-1, 1:-1] = (p_novo[2:-1, 1:-1] - p_novo[1:-2, 1:-1])/(dx) #7x6
            dpdy[1:-1, 1:-1] = (p_novo[1:-1, 2:-1] - p_novo[1:-1, 1:-2])/(dy) #8x5

            u_novo[1:-1, 1:-1] = u_[1:-1, 1:-1] - dpdx[1:-1, 1:-1] * dt # 7x6
            v_novo[1:-1, 1:-1] = v_[1:-1, 1:-1] - dpdy[1:-1, 1:-1] * dt # 8x5
            
            u_novo, v_novo = condicoes_contorno_velocidades(u_novo, v_novo)

            u[:], v[:], p[:] = u_novo, v_novo, p_novo
            
            t += dt
            # u_maximo = np.max(u[1:-1, 1:-1])
            # v_maximo = np.max(v[1:-1, 1:-1])

            informacoes(it, passos_tempo, t, t_final, u_maximo, v_maximo, normal2)
            print(f"X/r={estima_x1r(u, X, bfs_x, dx)}")

            u_vert, v_vert = (u[:, :-1] + u[:, 1:])/2, (v[1:, :] + v[:-1, :])/2
            u_vert[:bfs_x+1, :bfs_y+1] = 0
            v_vert[:bfs_x+1, :bfs_y+1] = 0
            V = (u_vert**2+v_vert**2)**(1/2)
            p_vert = (p_novo[:-1, 1:] + p_novo[1:, 1:] + p_novo[:-1, :-1] + p_novo[1:, :-1])/4
            if plotar_evolucao and it % plotar_a_cada == 0:
                plt.plot(X[:bfs_x+1, bfs_y], Y[:bfs_x+1, bfs_y], c='w')
                plt.plot(X[bfs_x, :bfs_y+1], Y[bfs_x, :bfs_y+1], c='w')
                plt.contourf(X, Y, V, levels=50, cmap='jet')
                plt.colorbar()
                plt.draw()
                plt.pause(0.01)
                plt.clf()
    except KeyboardInterrupt:
        pass
    fim = perf_counter()
    tempo = fim - inicio
    salvar(u, v, p)
    plt.show()
    return X, Y, u_vert, v_vert, V, p_vert, tempo