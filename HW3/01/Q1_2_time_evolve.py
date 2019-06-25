# -*- coding: utf-8 -*-
from cmath import sin, cos, exp, sqrt, pi # 导入常用函数，但要从复数库里导入哦！

######################################################################################
### 以下是同样的函数
def eigen_solver(a, b, c, f): ## 求三对角矩阵特征向量, a,b,c 为-1,0,1级对角线
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
            break
        else:
            f_pre = f.copy()
            Nloop += 1
    return 1/(f_abs_max), f # 返回：本征值 和 本征矢
### 同样的函数结束
######################################################################################

## 新函数
def time_evolve(u, u_): # 按薛定谔方程演化一步的函数 u->u_
    for i in range(N): # 此时刻 H 的对角元 b[i]，为后文做准备
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
    
    ## 吸收虚假信号
    for i in range(N):
        if abs((i-Nx)*dx) > 0.75*x_max:
            u_[i] *= exp(-(abs((i-Nx)*dx)-0.75*x_max)**2/0.04)
    ## 归一化
    unorm = sqrt(sum([abs(u_[i])**2 for i in range(N)]))
    for i in range(N):
        u_[i] /= unorm

#########################################################################
## 开始第(2)小问
if __name__ == '__main__':
    x_max, dx = 2000, 0.1
    Nx = int(x_max / dx)
    x = [dx*i for i in range(-Nx, Nx+1)] # x的格点
    # 初始化-1,0,1级对角
    a = [0] + [-0.5/dx**2 for i in range(len(x)-1)]
    c = [-0.5/dx**2 for i in range(len(x)-1)]+ [0]
    b = [1./dx**2 -1 / sqrt(((i-Nx)*dx)**2 + 2) + 0.48  for i in range(len(x))]
    lam, u = eigen_solver(a, b, c, [1 for i in range(len(x))]) # 求解！
    print('和 0.48 的偏差 ΔE = %.8f'%lam.real)

    unorm = sqrt(sum([abs(u[i])**2 for i in range(len(x))]))
    for i in range(len(x)):
        u[i] /= unorm
    u0 = u.copy() # 拷贝一下初始波函数
    ## 以上是一样的内容：求初始的氢原子 1s 态波函数
    #########################################################################
    
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
    print('即将开始演化波函数，真的很慢，作为测试可尝试将 x_max 改为200或20')
    for it in range(len(t_list)): ## 开始时间演化
        time_evolve(u, u_) # 按时间演化一步
        u, u_ = u_, u # 直接交换两个变量引用，省储存空间的做法
        t += dt
        P0t.append(abs(sum([u0[i].conjugate()*u[i] for i in range(N)]))**2)
        u_list.append(u[:])
        if it%100 == 0:
            print('正在演化波函数，已经进行 %4d / %4d 步'%(it, len(t_list)))
    
    ## 调取末态波函数，计算 k 上投影
    u_fin = u_list[-2] # 其实最后多演化了一步，调取u_list[-2]为末态波函数
    p_fin = sum([u0[i].conjugate()*u_fin[i] for i in range(N)]) # 基态投影系数
    u_fin = [u_fin[i] - p_fin*u0[i] for i in range(len(u_fin))] # 减去基态函数

    ## 下面考虑在 k 本征态投影
    k_max, dk = 2.5, 0.01
    Nk = (int)(k_max/dk)
    k = [i*dk for i in range(-Nk, Nk+1)]
    ## 下式里 φ_k*=1/sqrt(dx) exp(-ikx), 除以sqrt(dx)目的是考虑到x格点的归一化
    Pk = [abs(sum([(1/sqrt(dx)*exp(-1j*k[ik]*(i-Nx)*dx))*u_fin[i] for i in range(N)]))**2 for ik in range(len(k))]

    #########################################################################
    ## 使用 matlibplot 作图
    from matplotlib import pyplot as plt
    ## 先画 1s 态布局数随时间演化
    plt.plot(t_list, P0t[:-1], label=r'$P_0(t)$')
    plt.xlabel('time')
    plt.ylabel(r'$P_0(t)$')
    plt.xlim([t_list[0], t_list[-1]])
    plt.legend()
    plt.savefig('Q1_2_P0t.pdf')
    plt.show()
    
    ## 再画末态动量谱
    plt.plot(k, Pk, label=r'$P(k)$')
    plt.xlabel(r'$k$')
    plt.ylabel(r'$P(k)$')
    plt.xlim([-k_max, k_max])
    plt.legend()
    plt.savefig('Q1_2_Pk.pdf')
    plt.show()
    plt.plot(k, Pk, label=r'$P(k)$')
    plt.yscale('log')
    plt.xlabel(r'$k$')
    plt.ylabel(r'$P(k)$')
    plt.xlim([-k_max, k_max])
    plt.legend()
    plt.savefig('Q1_2_Pk_logy.pdf')
    plt.show()