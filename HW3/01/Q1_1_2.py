from cmath import sin, cos, exp, sqrt, pi
def eigen_solver(a, b, c, f): ## 求三对角矩阵特征向量
    N = len(b)
    m = a.copy()
    for i in range(1, N): # 用追赶法求三对角矩阵方程 
        m[i] = a[i] / b[i-1]
        b[i] -= m[i] * c[i-1]
    Nloop = 0
    f_pre = f.copy()
    while True: # 不断对方程的解迭代
        for i in range(1, N): # 追
            f[i] -= m[i] * f[i-1]
        f[N-1] /= b[N-1]
        for i in range(N-2, -1, -1): # 赶
            f[i] = (f[i] - c[i] * f[i+1]) / b[i]
        # 此时 f[i] 是方程的解
        f_abs_max = -1
        for i in range(N):
            if abs(f[i]) > abs(f_abs_max):
                f_abs_max = f[i] # 找到|f|最大的项
        for i in range(N):
            f[i] /= (f_abs_max) # 规格化矢量
        if sum([(abs(f_pre[i] - f[i]))**2 for i in range(N)]) < 1e-6: # 两次迭代足够接近->判停
            print(Nloop)
            break
        else:
            f_pre = f.copy()
            Nloop += 1
#         if Nloop % 100 == 0:
#             print(bp)
    return 1/(f_abs_max), f # 返回：本征值 和 本征矢
    
x_max, dx = 2000, 0.1
Nx = int(x_max / dx)
x = [dx*i for i in range(-Nx, Nx+1)]
a = [0] + [-0.5/dx**2 for i in range(len(x)-1)]
c = [-0.5/dx**2 for i in range(len(x)-1)]+ [0]
b = [1./dx**2 -1 / sqrt(((i-Nx)*dx)**2 + 2) + 0.48  for i in range(len(x))]
lam, u = eigen_solver(a, b, c, [1 for i in range(len(x))])
print(lam)

unorm = sqrt(sum([abs(u[i])**2 for i in range(len(x))]))
for i in range(len(x)):
    u[i] /= unorm
from matplotlib import pyplot as plt
plt.plot(x, u)
plt.show()
u0 = u.copy()

##########################################
E0, Npulse = sqrt(1 / 3.5094448314), 18 
ome = 1
T = 2*Npulse*pi / ome
dt = 0.05
Nt = int(T/dt)
t_list = [i/Nt*T for i in range(Nt+1)]
t = dt / 2
N = len(x)
u = u0.copy() # 储存初始的 u
u_ = u.copy() # 用作储存演化下一刻的 u
m, be = a.copy(), a.copy()
P0t = [abs(sum([u0[i].conjugate()*u[i] for i in range(N)]))**2] # 储存 P_0(t)
u_list = [u[:]]
for it in t_list:
    for i in range(N):
        b[i] = 1./dx**2 -1 / sqrt(((i-Nx)*dx)**2 + 2) + (i-Nx)*dx*E0*sin(ome*t/(2*Npulse))**2*sin(ome*t)
    ## 先求 u_=(1-1/2iHΔt) u
    u_[0] = (1 - 0.5j*dt*b[i]) * u[0] - 0.5j*dt*c[i] * u[1]
    u_[N-1] = - 0.5j*dt*a[i] * u[N-2] + (1 - 0.5j*dt*b[i]) * u[N-1] 
    for i in range(1, N-1):
        u_[i] = - 0.5j*dt*a[i] * u[i-1] + (1 - 0.5j*dt*b[i]) * u[i] - 0.5j*dt*c[i] * u[i+1]
    ## 再用追赶法求 (1+1/2iHΔt) u_' = u_
    be[0] = (1 + 0.5j*dt*b[0])
    for i in range(1, N): # 追
        be[i] = (1 + 0.5j*dt*b[i])
        m[i] = 0.5j*dt*a[i] / be[i-1]
        be[i] -= m[i] * 0.5j*dt*c[i-1]
        u_[i] -= m[i] * u_[i-1]
    u_[N-1] /= (1 + 0.5j*dt*b[i-1])
    for i in range(N-2, -1, -1): # 赶
        u_[i] = (u_[i] - 0.5j*dt*c[i] * u_[i+1]) / be[i]
    ## 此时u_[i] 就是演化下一刻的波函数
    
    for i in range(N):
        if abs((i-Nx)*dx) > 0.75*x_max:
            u_[i] *= exp(-(abs((i-Nx)*dx)-0.75*x_max)**2/0.04)
    ## 归一化
    unorm = sqrt(sum([abs(u_[i])**2 for i in range(N)]))
    for i in range(N):
        u_[i] /= unorm
    
    u, u_ = u_, u # 直接交换两个变量引用，省空间做法
    t += dt
    P0t.append(abs(sum([u0[i].conjugate()*u[i] for i in range(N)]))**2)
    u_list.append(u[:])
    if it%20 == 0:
        print('已经 %d'%it)

plt.plot(P0t)
plt.show()