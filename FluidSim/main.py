import numpy as np
import matplotlib.pyplot as plt

# Sempre inicializar: d2udx2, d2udy2, du2dx, duvdy, d2vdx2, d2vdy2, dv2dy, duvdx, dFdx, dGdy, dpdx, dpdy = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))

def gerar_malha(comprimento, altura, nx, ny):
    x = np.linspace(0, comprimento, nx)
    y = np.linspace(0, altura, ny)
    X, Y = np.meshgrid(x, y)
    X = np.transpose(X)
    Y = np.transpose(Y)
    return X, Y
def condicoes_iniciais(u0, v0, p0, nx, ny):
    u = u0*np.ones((nx, ny))
    v = v0*np.ones((nx, ny))
    p = p0*np.ones((nx, ny))
    return u, v, p
def condicao_contorno_duto(u, v):
    u[:, 0] = 0
    u[:, -1] = 0
    u[-1, 1:-1] = u[-2, 1:-1]
    u[0, 1:-1] = 1

    v[:, 0] = 0
    v[:, -1] = 0
    v[0, :] = 0
    v[-1, 1:-1] = v[-2, 1:-1]
    return u, v
def derivada_central(campo, dx, dy, nx, ny):
    dcdx, dcdy = np.zeros((nx, ny)), np.zeros((nx, ny))
    dcdx[1:-1, 1:-1] = (campo[2:, 1:-1] - campo[:-2, 1:-1]) / (2 * dx)
    dcdy[1:-1, 1:-1] = (campo[1:-1, 2:] - campo[1:-1, :-2]) / (2 * dy)
    return dcdx, dcdy
def derivada_segunda(campo, dx, dy, nx, ny):
  d2cdx2, d2cdy2 = np.zeros((nx, ny)), np.zeros((nx, ny))
  d2cdx2[1:-1, 1:-1] = (campo[2:, 1:-1] -2*campo[1:-1, 1:-1] + campo[:-2, 1:-1]) / dx**2
  d2cdy2[1:-1, 1:-1] = (campo[1:-1, 2:] -2*campo[1:-1, 1:-1] + campo[1:-1, :-2]) / dy**2
  return d2cdx2, d2cdy2
def desenhar_grafico_u(u, X, Y):
    plt.contourf(X, Y, u)
    plt.title('Velocidade horizontal u')
    plt.colorbar()
    plt.show()
    return
def passo(u, v, p, dx, dy, Re, dt, tol, N_p, nx, ny):
    #F
    d2udx2, d2udy2 = derivada_segunda(u, dx, dy, nx, ny)
    du2dx = derivada_central(u**2, dx, dy, nx, ny)[0]
    duvdy = derivada_central(u*v, dx, dy, nx, ny)[1]

    convectivo_x = du2dx + duvdy
    difusivo_x = 1/Re * (d2udx2 + d2udy2)

    u_t = u + dt * (difusivo_x - convectivo_x)

    #G

    d2vdx2, d2vdy2 = derivada_segunda(v, dx, dy, nx, ny)
    dv2dy = derivada_central(v**2, dx, dy, nx, ny)[1]
    duvdx = derivada_central(u*v, dx, dy, nx, ny)[0]

    convectivo_y = dv2dy + duvdx
    difusivo_y = 1/Re * (d2vdx2 + d2vdy2)

    v_t = v + dt * (difusivo_y - convectivo_y)

    #fonte

    dFdx = derivada_central(u_t, dx, dy, nx, ny)[0]
    dGdy = derivada_central(v_t, dx, dy, nx, ny)[1]

    fonte = (dFdx + dGdy)*dt

    #pressao

    c = 0
    erro = 1
    while erro > tol and c < N_p:
        c += 1
        p_old = np.copy(p)
        p[1:-1, 1:-1] = (
            (
                dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1])
                +
                dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])
                -
                dx**2 * dy**2 * fonte[1:-1, 1:-1]
            )
            /
            (2 * (dx**2 + dy**2))
            )
        erro = np.linalg.norm(p - p_old, ord=np.inf) / np.linalg.norm(p)
        p[-1, :] = p[-2, :]
        p[:, -1] = p[:, -2]
        p[:, 0] = p[:, 1]
        p[0, :] = p[1, :]
    dpdx, dpdy = derivada_central(p, dx, dy, nx, ny)
    u = u_t - dt*dpdx
    v = v_t - dt*dpdy
    u, v = condicao_contorno_duto(u, v)
    return u, v, p
    