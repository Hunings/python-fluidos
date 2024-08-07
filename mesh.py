import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

H = 2  # altura do tubo
h = 1.0297 # altura do ressalto
Li = 5*h # largura do ressalto
L0 = 30*h # largura do tubo depois do ressalto
L = Li + L0 # largura total
bfsY = 1.0344827586206897 # altura do degrau exata no ponto i = 0, j = 15
bfsX = 5.163252148997135 # largura exata do degrau no ponto i = 50, j = 15

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
u0 = np.zeros((n, m))
v0 = np.zeros((n, m))
p = np.zeros((n, m))

#perfil parabólico da velocidade e condições

u0[1:51, :] = 23*(Y[1:51, :] - bfsY)*(H - Y[1:51, :])
u0[51:, :] = 6*Y[51:, :]*(H - Y[51:, :])

# condições em t0

u0[0, :] = 8
p[:, :] = 10e5
v0[0, :] = 0

# degrau

deg[0:51, 0:16] = 0
u0[0:51, 0:16] = 0

# plotando u0
plt.figure(figsize=(18, 4))
plt.pcolormesh(X, Y, u0, cmap='viridis')
plt.colorbar()
plt.show()

# plotando p

plt.figure(figsize=(18, 4))
plt.pcolormesh(X, Y, np.multiply(p, deg), cmap='viridis')
plt.show()

def gauss_seidel(p, f, dy, dx, tol, max_iter):
  c = 0
  error = 1.
  while (error > tol) or c==max_iter:
      c += 1
      p_old = np.copy(p)
      for i in range(1, n-1):
          for j in range(1, m-1):
              if deg[i, j] != 0:  # ignora células dentro do degrau
                  p[i, j] = ((dy**2 * (p[i+1, j] + p[i-1, j]) +
                              dx**2 * (p[i, j+1] + p[i, j-1]) -
                              dx**2 * dy**2 * f[i, j]) /
                            (2 * (dx**2 + dy**2)))
      error = np.linalg.norm(p - p_old, ord=np.inf)
gauss_seidel(p, np.ones_like(p), dy, dx, 1e-6, 1000)

plt.figure(figsize=(18, 4))
plt.pcolormesh(X, Y, p, cmap='viridis')
plt.show()

# plota a malha em 3d
'''fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z, color='black', linewidth=0.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Malha computacional')
plt.scatter(0, 5, c='red', s=0.04)
plt.show()'''