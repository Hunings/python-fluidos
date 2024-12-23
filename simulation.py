import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

H = 2  # altura do tubo
h = 1.0297 # altura do ressalto
Li = 5*h # largura do ressalto
L0 = 30*h # largura do tubo depois do ressalto
L = Li + L0 # largura total

# número de colunas e linhas da malha
n = 350
m = 30
dx = L/(n-1)
dy = H/(m-1)

# cria malha

x = np.linspace(0, L, n)
y = np.linspace(0, H, m)

X0, Y0 = np.meshgrid(x, y)

# transpõe as matrizes e cria as matrizes deg, u0, v0, p0
X = np.transpose(X0)
Y = np.transpose(Y0)
deg = np.ones_like(X)
u = np.zeros((n, m))
v = np.zeros((n, m))
p = np.zeros((n, m))

# condições de contorno

u[0, :] = 8
p[:, :] = 1e5
v[0, :] = 0
u[:, -1] = 0
u[:, 0] = 0

# perfil parabólico da velocidade e condições

u[:, :] = 6*Y[:, :]*(H - Y[:, :])

# degrau

deg[0:51, 0:16] = 0
u[0:51, 0:16] = 0

# gauss seidel

def gauss_seidel(p, f, dy, dx, tol, max_iter):
  c = 0
  error = 1.
  while (error > tol) and (c < max_iter):
      c += 1
      p_old = np.copy(p)
      for i in range(1, n-):
          for j in range(1, m-):
              if deg[i, j] != 0:  # ignora células dentro do degrau
                  p[i, j] = ((dy**2 * (p[i+1, j] + p[i-1, j]) +
                              dx**2 * (p[i, j+1] + p[i, j-1]) -
                              dx**2 * dy**2 * f[i, j]) /
                            (2 * (dx**2 + dy**2)))
      error = np.linalg.norm(p - p_old, ord=np.inf)
  print(c, error)

#função qualquer 
gauss_seidel(p=p, f=np.ones_like(p), dy=dy, dx=dx, tol=1e-4, max_iter=5)


# plotando 
while True:
  resp = str(input('O que deseja plotar? '))
  if resp in 'Uu0':
    plt.figure(figsize=(18, 4))
    plt.pcolormesh(X, Y, u, cmap='viridis')
    plt.colorbar()
    plt.show()
    break
  elif resp in 'pP0':
    plt.figure(figsize=(18, 4))
    plt.pcolormesh(X, Y, p, cmap='viridis')
    plt.colorbar()
    plt.show()
    break
  else:
     pass

# plota a malha em 3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, deg, color='black', linewidth=0.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Malha computacional')
plt.scatter(0, 5, c='red', s=0.04)
plt.show()
print(dx, dy)