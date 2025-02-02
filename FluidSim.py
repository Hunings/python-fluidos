import numpy as np
import matplotlib.pyplot as plt

#variáveis globais

altura = 10.
comprimento = 10.
nx = 10
ny = 10
dx = comprimento / nx
dy = altura / ny
Re = 1
N_p = 10000
dt = 0.001
tol = 1e-6
t_final = 1

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

def F(u, v):
  dudx, dudy = derivada_central(u, dx, dy)
  d2udx2, d2udy2 = derivada_segunda(u, dx, dy)
  du2dx, _ = derivada_central(u**2, dx, dy) # termos não lineares
  _, duvdy = derivada_central(u*v, dx, dy) 

  convectivo = du2dx + duvdy
  difusivo = 1/Re * (d2udx2 + d2udy2)

  return u + dt * (difusivo - convectivo)

def G(u, v):
  dvdx, dvdy = derivada_central(v, dx, dy)
  d2vdx2, d2vdy2 = derivada_segunda(v, dx, dy)
  _, dv2dy = derivada_central(v**2, dx, dy) # termos não lineares
  duvdx, _ = derivada_segunda(u*v, dx, dy)

  convectivo = duvdx + dv2dy
  difusivo = 1/Re * (d2vdx2 + d2vdy2)

  return v + dt * (difusivo - convectivo)

def f(F, G):
  dFdx, _ = derivada_central(F, dx, dy)
  _, dGdy = derivada_central(G, dx, dy)
  return (dFdx + dGdy)*dt

def pressao(u, v, p):
  c = 0
  erro = 1
  erros = np.array([])
  F_ = F(u, v) # 
  G_ = G(u, v)
  fonte = f(F_, G_)[1:-1, 1:-1] # será melhor retornar a matriz certa na função?
  while (erro > tol) and (c < N_p):
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
      erro = np.linalg.norm(p - p_old, ord=np.inf)
      erros = np.append(erros, erro)
      p[-1, :] = 0 #p[-2, :] 
      p[:, -1] = 0 #p[:, -2] 
      p[:, 0] = 0 #p[:, 1]  
      p[0, :] = 0 #p[1, :]
  return p

def passo(u, v, p):
  p_next = pressao(u, v, p) # itera e calcula a próxima pressão, com cc já
  dpdx, dpdy = derivada_central(p_next, dx, dy)

  F_ = F(u, v)
  G_ = G(u, v)

  u_next = F_ - dt*dpdx
  v_next = G_ - dt*dpdy

  u_next[0, :] = 1.0 
  u_next[-1, :] = u_next[-2, :]  
  u_next[:, 0] = 0.0  
  u_next[:, -1] = 0.0  
  
  v_next[0, :] = 0.0
  v_next[-1, :] = v_next[-2, :]  
  v_next[:, 0] = 0.0
  v_next[:, -1] = 0.0
  return u_next, v_next, p_next 

#pseudo 

def init(u, v, p, t_final):
  n = 0
  t = 0
  while t < t_final:
    u_next, v_next, p_next  = passo(u, v, p)
    if n < 10:
      plt.contourf(X, Y, p_next)
      plt.colorbar()
      plt.title('p')
      plt.show()
      plt.contourf(X, Y, u_next)
      plt.title('u')
      plt.colorbar()
      plt.show()
      plt.title('v')
      plt.contourf(X, Y, u_next)
      plt.colorbar()
      plt.show()
    u, v, p = p_next, u_next, v_next
    t += dt
    n += 1
  return

p_ = pressao(u0, v0, p0)
plt.pcolormesh(X, Y, p_)
plt.colorbar()
plt.show()
init(u0, v0, p0, t_final)