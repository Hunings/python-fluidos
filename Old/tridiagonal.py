import numpy as np
'''
A função gauss resolve o sistema de uma matriz tridiagonal com diagonal principal -2, e diagonais secundárias 1 de tamanho n

utiliza o método de Gauss-Seidel e realiza c iterações ou até que o erro seja menor que o desejado

'''
def gauss(b, x_old, tol, n, itmax):
    x = x_old
    eps = 1.
    c = 0 
    while eps > tol or c == itmax: # para quando eps atinge valor menor ou igual a tol ou c atinge itmax
        c += 1
        x_old = np.copy(x)
        x[0] = (-x_old[1] + b[0])/(-2) # linha 1 é diferente das demais
        for j in range(1, n-1):
            x[j] = (-x[j-1] - x_old[j+1] + b[j])/(-2) # linhas do meio da matriz
        x[n-1] = (-x[n-2] + b[n-1])/(-2) # linha final da matriz
        eps = np.linalg.norm(x - x_old, ord=np.inf) # calcula a norma infinito entre x e x_old, para calcular erro
    print(f"O código realizou {c} iterações")
    print(f"O erro é de {eps}")
    return x

n = int(input('Digite o tamanho da matriz: '))
b = np.array([])
x_old = np.array([])
for i in range(n):
    valor = float(input(f"Elemento {i+1} do vetor b: "))
    b = np.append(b, valor)

while True:
    try:
        q = str(input('Quer digitar um x inicial? (zero por padrão) [S/N] ')).strip().upper()
        if q not in ['S', 'N']:
            raise ValueError("Resposta inválida. Por favor, digite 'S' ou 'N'.")
        break  # Sai do loop se a resposta for válida
    except ValueError as ve:
        print(ve)
if q == 'S':
    for i in range(n):
        valor = float(input(f"Elemento {i+1} do x inicial: "))
        x_old = np.append(x_old, valor)
else:
    x_old = np.zeros(n)
print(b)
print(x_old)
tol = float(input('Digite a tolerância: '))
itmax = int(input('Digite um número máximo de iterações: '))
sol = gauss(b, x_old, tol, n, itmax)
print(f"A solução obtida foi: {sol}")