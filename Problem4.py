import math

dx = 0.1
dt = 0.1
nx = 11
nt = 11
lambda_sq = (dt / dx) ** 2

p = [[0.0 for _ in range(nt)] for _ in range(nx)]
    
for i in range(nx):
    x = i * dx
    p[i][0] = math.cos(2 * math.pi * x)
    
for j in range(nt):
    p[0][j] = 1.0
    p[nx-1][j] = 2.0

#Expand u(xi , t1)
for i in range(1, nx-1):
    x = i * dx
    p[i][1] = (0.5 * lambda_sq * math.cos(2 * math.pi * (x - dx)) + 
               (1 - lambda_sq) * math.cos(2 * math.pi * x) + 
               0.5 * lambda_sq * math.cos(2 * math.pi * (x + dx)) + 
               0.1 * 2 * math.pi * math.sin(2 * math.pi * x))
    
for j in range(1, nt-1):
    for i in range(1, nx-1):
        p[i][j+1] = p[i-1][j] + p[i+1][j] - p[i][j-1]
        
print("x\t t\t p(x,t)")
for j in range(nt):
    t = j * dt
    for i in range(nx):
        x = i * dx
        print(f"{x:.1f}\t{t:.1f}\t{p[i][j]:.6f}")