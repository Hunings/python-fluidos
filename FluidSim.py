import numpy as np
import matplotlib.pyplot as plt

#variáveis globais

altura = 100
comprimento = 100
nx = 20
ny = 20
dx = comprimento / nx
dy = altura / ny
Re = 100
N_p = 10000
tol = 1e-4
u_max = 1.5
v_max = 1.5
tau = 0.3
dt = tau*min(Re/2*(1/dx**2 + 1/dy**2), dx/u_max, dy/u_max)
t_final = 100*dt

x = np.linspace(0, comprimento, nx)
y = np.linspace(0, altura, ny)

X, Y = np.meshgrid(x, y)
X = np.transpose(X)
Y = np.transpose(Y)

# arrays iniciais

u0 = 1*np.ones((nx, ny))
v0 = 0*np.ones((nx, ny))
p0 = 1*np.ones((nx, ny))

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

def convecF(u, v):
  du2dx = np.zeros_like(u)
  duvdy = np.zeros_like(u)
  du2dx[1:-1, 1:-1] = 1/dx * ((u[1:-1, 1:-1] + u[2:, 1:-1])**2 / 2 - (u[:-2, 1:-1] + u[1:-1, 1:-1])**2 /2)
  duvdy[1:-1, 1:-1] = 1/dy * ((v[1:-1, 1:-1] + v[2:, 1:-1])*(u[1:-1, 1:-1] + u[1:-1, 2:])/4 - (v[1:-1, :-2] + v[2:, :-2])*(u[1:-1, :-2] + u[1:-1, 1:-1])/4)
  return du2dx + duvdy

# duas funções ou apenas uma?

def convecG(u, v):
  dv2dy = np.zeros_like(v)
  duvdx = np.zeros_like(u)
  dv2dy[1:-1, 1:-1] = 1/dy * ((v[1:-1, 1:-1] + v[2:, 1:-1])**2 / 2 - (v[:-2, 1:-1] + v[1:-1, 1:-1])**2 /2)
  duvdx[1:-1, 1:-1] = 1/dx * ((u[1:-1, 1:-1] + u[1:-1, 2:])*(v[1:-1, 1:-1] + v[2:, 1:-1])/4 - (u[:-2, 1:-1] + u[:-2, 2:])*(v[:-2, 1:-1] + v[1:-1, 1:-1])/4)
  return dv2dy + duvdx

def F(u, v):
  d2udx2, d2udy2 = derivada_segunda(u, dx, dy)

  convectivo = convecF(u, v)
  difusivo = 1/Re * (d2udx2 + d2udy2)

  return u + dt * (difusivo - convectivo)

def G(u, v):
  d2vdx2, d2vdy2 = derivada_segunda(v, dx, dy)

  convectivo = convecG(u, v)
  difusivo = 1/Re * (d2vdx2 + d2vdy2)

  return (v + dt * (difusivo - convectivo))

def f(F, G):
  dFdx, _ = derivada_central(F, dx, dy)
  _, dGdy = derivada_central(G, dx, dy)
  return (dFdx + dGdy)*dt

def pressao(u, v, p):
  c = 0
  erro = 1
  F_ = F(u, v)
  G_ = G(u, v)
  fonte = f(F_, G_)[1:-1, 1:-1] # será melhor retornar a matriz certa na função?

  while erro > tol and c < N_p:
      c += 1
      p_old = np.copy(p)
      p[1:-1, 1:-1] = (
          (
              dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1])
              +
              dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])
              -
              dx**2 * dy**2 * fonte
          )
          /
          (2 * (dx**2 + dy**2))
          )
      p[-1, :] = p[-2, :]
      p[:, -1] = p[:, -2]
      p[:, 0] = p[:, 1]
      p[0, :] = p[1, :] # deve ser antes do erro

      erro = np.linalg.norm(p - p_old, ord=np.inf) / np.linalg.norm(p)
      print("Diferença p: ", erro)
  return p

def passo(u, v, p):
  p_next = pressao(u, v, p)
  dpdx, dpdy = derivada_central(p_next, dx, dy)

  F_ = F(u, v)
  G_ = G(u, v)

  u_next = F_ - dt*dpdx
  v_next = G_ - dt*dpdy

  u_next[:, 0] = 0.0
  u_next[:, -1] = 0.0
  u_next[-1, :] = u_next[-2, :]
  u_next[0, :] = 1

  v_next[:, 0] = 0.0
  v_next[:, -1] = 0.0
  v_next[0, :] = 0
  v_next[-1, :] = v_next[-2, :]

  return u_next, v_next, p_next

def init(u, v, p, t_final):
  n = 0
  t = 0
  c = 0
  while t < t_final:
    u_next, v_next, p_next  = passo(u, v, p)
    u, v, p = u_next, v_next, p_next
    t += dt
    n += 1

  plt.contourf(X, Y, u)
  plt.title('u')
  plt.colorbar()
  plt.show()

  plt.title('v')
  plt.contourf(X, Y, v)
  plt.colorbar()
  plt.show()

  plt.title('Vetores')
  plt.quiver(X, Y, u, v)
  plt.show()
  return

p_ = pressao(u0, v0, p0)
plt.pcolormesh(X, Y, p_)
plt.colorbar()
plt.show()
init(u0, v0, p0, t_final)