import numpy as np
import matplotlib.pyplot as plt
import simulacao as s

nx = 80
ny = 80
comprimento = 24
altura = 2

u, v, p = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))

u, v, p, it = s.carregar(u, v, p)

x = np.linspace(0, comprimento, nx)
y = np.linspace(0, altura, ny)

X, Y = np.meshgrid(x, y)
X, Y, = X.T, Y.T

plt.figure(figsize=(15, 2))
v = (u**2 + v**2)**(0.5)
plt.contourf(X, Y, v, levels=100, cmap='jet')
plt.colorbar()
plt.show()