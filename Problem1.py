import math

def f(x , y):
    return x * y

h = 0.1 * math.pi 
k = 0.1 * math.pi 
n = 9
m = 4
nx = n + 2
ny = m + 2

u = [[0.0 for _ in range(ny)] for _ in range(nx)]

for j in range(ny):
    y = j * k
    u[0][j] = math.cos(y) 
    u[n + 1][j] = -math.cos(y)
    
for i in range(nx):
    x = i * h
    u[i][0] = math.cos(x)
    u[i][m + 1] = 0.0

N = n * m
U = [0.0 for _ in range(N)]
F = [0.0 for _ in range(N)]

for j in range(1, m + 1):
    for i in range(1, n + 1):
        l = i + n * (j - 1) - 1
        x = i * h
        y = j * k
        F[l] = h**2 * f(x , y)
        
        if i == 1:
            F[l] -= u[0][j]
        if i == n:
            F[l] -= u[n + 1][j]
        if j == 1:
            F[l] -= u[i][0]
        if j == m:
            F[l] -= u[i][m + 1]

max_iterations = 1000
tolerance = 1e-6
for iteration in range(max_iterations):
    max_diff = 0.0
    for j in range(1, m + 1):
        for i in range(1, n + 1):
            l = i + n * (j - 1) - 1
            
            sum_terms = 0.0
            if i > 1:
                sum_terms += U[l - 1]
            if i < n:
                sum_terms += U[l + 1]
            if j > 1:
                sum_terms += U[l - n]
            if j < m:
                sum_terms += U[l + n]
            new_u = (F[l] - sum_terms) / -4.0  # a_{ll} = -4
            
            max_diff = max(max_diff, abs(new_u - U[l]))
            U[l] = new_u
            
    if max_diff < tolerance:
        break

for j in range(1, m + 1):
    for i in range(1, n + 1):
        l = i + n * (j - 1) - 1
        u[i][j] = U[l]

print("x\t y\t u(x,y)")
for j in range(ny):
    for i in range(nx):
        x = i * h
        y = j * k
        print(f"{x:.3f}\t{y:.3f}\t{u[i][j]:.6f}")