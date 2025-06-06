dr = 0.1
dt = 0.5
K = 0.1
nr = 6
nt = 21
alpha_sq = 1 / (4 * K)
lambda_val = alpha_sq * dt / (dr * dr)  # λ = 20

T = [[0.0 for _ in range(nt)] for _ in range(nr)]

for i in range(nr):
    r = 0.5 + i * dr
    T[i][0] = 200 * (r - 0.5)
    
for j in range(nt):
    t = j * dt
    T[nr - 1][j] = 100 + 40 * t
 
def solve_tridiagonal(a, b, c, d):
    n = len(d)
    c_prime = [0] * (n - 1)
    d_prime = [0] * n
    x = [0] * n
    
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    for i in range(1, n - 1):
        denom = b[i] - a[i] * c_prime[i - 1]
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom
    d_prime[n - 1] = (d[n - 1] - a[n - 1] * d_prime[n - 2]) / (b[n - 1] - a[n - 1] * c_prime[n - 2])
    
    x[n - 1] = d_prime[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x

# (a) Forward-Difference Method
T_forward = [row[:] for row in T]
for j in range(nt - 1): # 0 ~ 19
    for i in range(1, nr - 1):  # 1 ~ 4
        r = 0.5 + i * dr
        dr_term = (T_forward[i+1][j] - T_forward[i-1][j]) / (2 * dr)
        d2r_term = (T_forward[i+1][j] - 2 * T_forward[i][j] + T_forward[i-1][j]) / (dr * dr)
        T_forward[i][j+1] = T_forward[i][j] + dt * (1 / alpha_sq) * (d2r_term + (1 / r) * dr_term)
    
    T_forward[0][j+1] = T_forward[1][j+1] / 1.3

# (b) Backward-Difference Method
T_backward = [row[:] for row in T]
for j in range(1, nt):  # j = 1, 2, ..., 20
    r = [0.5 + i * dr for i in range(nr)]
    
    a = [-lambda_val * (1 - dr / (2 * r[i])) for i in range(1, nr - 1)]
    b = [1 + 2 * lambda_val for _ in range(1, nr - 1)]
    c = [-lambda_val * (1 + dr / (2 * r[i])) for i in range(1, nr - 1)]
    d = [T_backward[i][j - 1] for i in range(1, nr - 1)]
    
    d[0] += lambda_val * (1 - dr / (2 * r[1])) * T_backward[0][j]
    d[-1] += lambda_val * (1 + dr / (2 * r[nr - 2])) * T_backward[nr - 1][j] 
    
    solution = solve_tridiagonal(a, b, c, d)
    for i in range(1, nr - 1):
        T_backward[i][j] = solution[i - 1]
    
    T_backward[0][j] = T_backward[1][j] / 1.3
    
# (c) Crank-Nicolson Algorithm
T_cn = [row[:] for row in T]
for j in range(nt - 1):
    r = [0.5 + i * dr for i in range(nr)]

    a = [-lambda_val / 2 * (1 - dr / (2 * r[i])) for i in range(1, nr - 1)]
    b = [1 + lambda_val for _ in range(1, nr - 1)]
    c = [-lambda_val / 2 * (1 + dr / (2 * r[i])) for i in range(1, nr - 1)]
    d = [0.0 for _ in range(1, nr - 1)]

    for i in range(1, nr - 1):
        r_i = 0.5 + i * dr
        d[i - 1] = (lambda_val / 2 * (1 - dr / (2 * r_i)) * T_cn[i - 1][j] +
                    (1 - lambda_val) * T_cn[i][j] +
                    lambda_val / 2 * (1 + dr / (2 * r_i)) * T_cn[i + 1][j])

    d[0] += lambda_val / 2 * (1 - dr / (2 * r[1])) * T_cn[0][j + 1]  # r=0.5
    d[-1] += lambda_val / 2 * (1 + dr / (2 * r[nr - 2])) * (T_cn[nr - 1][j] + T_cn[nr - 1][j + 1])  # r=1

    solution = solve_tridiagonal(a, b, c, d)
    for i in range(1, nr - 1):
        T_cn[i][j + 1] = solution[i - 1]
        
    T_cn[0][j + 1] = T_cn[1][j + 1] / 1.3
    
print("Forward-Difference Method:")
print("r\t t\t T(x,t)")
for j in range(nt):
    t = j * dt
    for i in range(nr):
        r = 0.5 + i * dr
        print(f"{r:.1f}\t{t:.1f}\t{T_forward[i][j]:.6f}")
        
print("\nBackward-Difference Method:")
print("r\t t\t T(r,t)")
for j in range(nt):
    t = j * dt
    for i in range(nr):
        r = 0.5 + i * dr
        print(f"{r:.1f}\t{t:.1f}\t{T_backward[i][j]:.6f}")
        
print("\nCrank-Nicolson Algorithm:")
print("r\t t\t T(r,t)")
for j in range(nt):
    t = j * dt
    for i in range(nr):
        r = 0.5 + i * dr
        print(f"{r:.1f}\t{t:.1f}\t{T_cn[i][j]:.6f}")