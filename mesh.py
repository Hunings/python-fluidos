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

# cria malha

x = np.linspace(0, L, n)
y = np.linspace(0, H, m)

X0, Y0 = np.meshgrid(x, y)

# transpõe as matrizes e cria as matrizes Z, u0, v0, p0
X = np.transpose(X0)
Y = np.transpose(Y0)
Z = np.ones_like(X)
u0 = np.zeros((n, m))
v0 = np.zeros((n, m))
p0 = np.zeros((n, m))

# loops para condições iniciais

for i in range(n):
  for j in range(m):
    if i == 0:
      u0[i, j] = 8
      p0[i, j] = 10e5
      v0[i, j] = 0
    elif i < 51:
      u0[i, j] = 23*(Y[i, j] - bfsY)*(H - Y[i, j])
      v0[i, j] = 0
      p0[i, j] = 10e5
    elif i > 50:
      #u0[i, j] = 4*Y[i, j]*(h - Y[i, j])/h**2 Nowruzi
      u0[i, j] = 6*Y[i, j]*(H - Y[i, j])
      v0[i, j] = 0
      p0[i, j] = 10e5
    if i < 51 and j < 16:
       Z[i, j] = 0
       u0[i, j] = 0

# plotando u0
plt.figure(figsize=(18, 4))
plt.pcolormesh(X, Y, np.multiply(u0, Z), cmap='viridis')
plt.colorbar()
plt.show()

# plotando p0

plt.figure(figsize=(18, 4))
plt.pcolormesh(X, Y, np.multiply(p0, Z), cmap='viridis')
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