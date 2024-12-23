import numpy as np
import matplotlib.pyplot as plt

#variáveis globais

altura = 10.
comprimento = 10.
nx = 5
ny = 6
dx = comprimento / nx
dy = altura / ny
Re = 1e-6
N_p = 100000
dt = 0.01
tol = 1e-12
t_final = 1

x = np.linspace(0, comprimento, nx)
y = np.linspace(0, altura, ny)

X, Y = np.meshgrid(x, y)
X = np.transpose(X)
Y = np.transpose(Y)

# arrays iniciais

u0 = np.ones((nx, ny))
v0 = 0*np.ones((nx, ny))

p0 = 1*np.zeros((nx, ny))

def derivada_central(campo, dx, dy):
    dcdx = np.zeros_like(campo)
    dcdy = np.zeros_like(campo)
    dcdx[1:-1, 1:-1] = (campo[2:, 1:-1] - campo[:-2, 1:-1]) / (2 * dx)
    dcdy[1:-1, 1:-1] = (campo[1:-1, 2:] - campo[1:-1, :-2]) / (2 * dy)
    return dcdx, dcdy

def derivada_segunda(campo, dx, dy):
  d2cdx2 = np.zeros_like(campo)
  d2cdy2 = np.zeros_like(campo)
  d2cdx2[1:-1, 1:-1] = (campo[2:, 1:-1] -2*campo[1:-1, 1:-1] + campo[:-2, 1:-1]) / dx**2
  d2cdy2[1:-1, 1:-1] = (campo[1:-1, 2:] - 2*campo[1:-1, 1:-1] + campo[1:-1, :-2]) / dy**2
  return d2cdx2, d2cdy2

def F(u):
  dudx, dudy = derivada_central(u, dx, dy)
  d2udx2, d2udy2 = derivada_segunda(u, dx, dy)
  u[0, :] = 1.0
  u[-1, :] = u[-2, :]
  u[:, 0] = 0
  u[:, -1] = 0
  return u + dt * (1/Re * (d2udx2 + d2udy2))

def G(v):
  dvdx, dvdy = derivada_central(v, dx, dy)
  d2vdx2, d2vdy2 = derivada_segunda(v, dx, dy)
  v[0, :] = 0
  v[-1, :] = v[-2, :]
  v[:, 0] = 0
  v[:, -1] = 0
  return v + dt * (1/Re * (d2vdx2 + d2vdy2))

def f(F, G):
  dFdx, _ = derivada_central(F, dx, dy)
  _, dGdy = derivada_central(G, dx, dy)
  return (dFdx + dGdy)*dt

def pressao(p, f, u, v, F, G): # essa função vai dentro do step depois. por quê?
  c = 0
  erro = 1
  erros = np.array([])
  while (erro > tol) and (c < N_p):
      c += 1
      p_old = np.copy(p)
      p[1:-1, 1:-1] = (
          (
              dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1])
              +
              dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])
              -
              dx**2 * dy**2 * f(F(u), G(v))[1:-1, 1:-1]
          )
          /
          (2 * (dx**2 + dy**2))
          )
      p[-1, :] = p[-2, :] #  LESTE
      p[:, -1] = p[:, -2] #  NORTE
      p[:, 0] = p[:, 1] # SUL
      p[0, :] = p[1, :] # OESTE
      erro = np.linalg.norm(p - p_old, ord=np.inf)
      erros = np.append(erros, erro)
  return p

def passo(un, vn, F, G, pn):
  p_next = pressao(pn, f, un, vn, F, G) # itera e calcula a próxima pressão
  dpdx, dpdy = derivada_central(p_next, dx, dy)
  u_next = F(un) - dt*dpdx
  v_next = G(vn) - dt*dpdy
  return p_next, u_next, v_next

def init(passo, u, v, p, t_final):
  n = 0
  t = 0
  while t < t_final:
    p_next, u_next, v_next = passo(u, v, F, G, p)
    if n < 10:
      plt.contourf(X, Y, p_next)
      plt.colorbar()
      plt.show()
      #plt.contourf(X, Y, u_next)
      #plt.colorbar()
      #plt.show()
      #plt.contourf(X, Y, u_next)
      #plt.colorbar()
      #plt.show()
    p, u, v = p_next, u_next, v_next
    t += dt
    n += 1
  return

init(passo, u0, v0, p0, t_final)