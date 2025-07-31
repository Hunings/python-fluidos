import numpy as np
import matplotlib.pyplot as plt
import simulacao as s

nx = 200
ny = 100
comprimento = 20
altura = 1

u, v, p = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))

u, v, p, it = s.carregar(u, v, p)

x = np.linspace(0, comprimento, nx)
y = np.linspace(0, altura, ny)

X, Y = np.meshgrid(x, y)
X, Y, = X.T, Y.T

plt.figure(figsize=(15, 2))
plt.contourf(X, Y, u, levels=100, cmap='jet')
plt.colorbar()
plt.show()