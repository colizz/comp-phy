# -*- coding: utf-8 -*-
from cmath import sin, cos, exp, sqrt, pi # 导入常用函数，但要从复数库里导入哦！

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

#########################################################################
## 开始第(1)小问
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
       
    ## 概率密度（*1/dx是考虑到格点归一化与实际归一化不同）
    u_prob = [(abs(ux)**2/dx).real for ux in u]
    #########################################################################
    ## 使用 matlibplot 作图
    from matplotlib import pyplot as plt
    print('我们只绘制 [-20, 20] 区间里的波函数啦')
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    plt.plot(x[int(len(x)/2)-200:int(len(x)/2)+200], u_prob[int(len(x)/2)-200:int(len(x)/2)+200],
            label=r'$|\psi_0(x)|^2$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|\psi_0(x)|^2$')
    plt.legend()
    plt.text(0.05, 0.9, r'$E_0=%.8f$'%(0.48+lam.real), transform=ax.transAxes)
    plt.savefig('Q1_1_psi0.pdf')