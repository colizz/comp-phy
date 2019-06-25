# -*- coding: utf-8 -*-
from cmath import sin, cos, exp, sqrt, pi, log # 导入常用函数，但要从复数库里导入哦！

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

## 第(2)问原封不动的函数
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
### 同样的函数结束
######################################################################################

## 开始第(4)小问
if __name__ == '__main__':
    x_max, dx = 600, 0.1 # 改成600啦！
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
    
    E0, Npulse = sqrt(0.01 / 3.5094448314), 4 ## 注意参数改了！！
    ome = 45.5633525316 / 400
    T = 2*Npulse*pi / ome
    dt = 0.5
    Nt = int(T/dt)
    t_list = [i/Nt*T for i in range(Nt+1)]
    t = dt / 2
    N = len(x)
    u = u0.copy() # 储存初始的 u
    u_ = u.copy() # 用作储存演化下一刻的 u
    m, be = a.copy(), a.copy()
    # 记 a 随时间演化
    a_list = [sum([u[i].conjugate()*(-((i-Nx)*dx)/(((i-Nx)*dx)**2 + 2)**(3/2) + E0*sin(ome*t/(2*Npulse))**2*sin(ome*t))*u[i] 
                       for i in range(N)])]
    print('即将开始演化波函数，这个还比较快')
    for it in range(len(t_list)): ## 开始时间演化
        time_evolve(u, u_) # 按时间演化一步
        u, u_ = u_, u # 直接交换两个变量引用，省储存空间的做法
        t += dt
        a_list.append(sum([u[i].conjugate()*(-((i-Nx)*dx)/(((i-Nx)*dx)**2 + 2)**(3/2) + E0*sin(ome*t/(2*Npulse))**2*sin(ome*t))*u[i] 
                           for i in range(N)]))
        if it%100 == 0:
            print('正在演化波函数，已经进行 %4d / %4d 步'%(it, len(t_list)))
    
    ## 计算高次谐波功率
    print('正在计算高次谐波功率，请稍后')
    Nome, Nt0, sigma = 20, 100, 15 # 共 20*Nome=400个ω上采样点, Nt0=100个t0上采样点
    ome_list = [i*ome/Nome for i in range(20*Nome+1)]
    t0_list = [i*T/Nt0 for i in range(Nt0+1)]
    # 计算 A(t0,ω):
    logA2 = [[log(abs( sum([exp(-1j*o*(it/Nt*T))*exp(-(it/Nt*T-t0)**2/(2*sigma**2))*a_list[it] for it in range(len(t_list)-1, -1, -1)]) ).real**2).real
               for t0 in t0_list] for o in ome_list]
    
    #########################################################################
    ## 使用 matlibplot 作图
    from matplotlib import pyplot as plt
    im = plt.imshow(logA2, extent=(0, T, 0, 20), aspect=T/20)
    plt.colorbar(im)
    plt.xlabel('time')
    plt.ylabel(r'$\omega\;/\;\omega_{400\,\mathrm{nm}}$')
    plt.savefig('Q1_4_logA2tome_t0.1.pdf')
    plt.show()