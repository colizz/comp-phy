from numpy import sqrt, pi, cos, sin, exp, log
from random import random

##################################
# pol = 'linear'
pol = 'elliptical'
##################################

def calc_elec(t, pol): ## 输入给定时刻 t 和偏振模式 pol，输出此时电场 E_x, E_y
    omet = ome * t
    if pol == 'linear':
        Ex = -A0*((cos(omet/8))**2*cos(omet) - sin(omet/4)*sin(omet)/8 )*ome
        Ey = 0
    elif pol == 'elliptical':
        Ex = -A0*((cos(omet/8))**2*cos(omet) - sin(omet/4)*sin(omet)/8 )*ome*2/sqrt(5)
        Ey = A0*((cos(omet/8))**2*sin(omet) + sin(omet/4)*cos(omet)/8 )*ome/sqrt(5)
    return Ex, Ey

def gen_sample(): ## 舍选法获得一个样品
    while True:
        t = (random()-0.5) * 4*T
        Ex, Ey = calc_elec(t, pol)
        Eabs = sqrt(Ex**2 + Ey**2) # 计算当前电场大小
        if WEmax * random() < 1/Eabs**(3/2) * exp(-2/(3*Eabs)): # 生成随机数是否在概率密度分布以下？
            r = -Ip/Eabs
            u, v = random(), random()
            v_perp = sqrt(-Eabs*log(u)) * cos(2*pi*v)
            return t, [r*Ex/Eabs, r*Ey/Eabs, -v_perp*Ey/Eabs, v_perp*Ex/Eabs] # return 格式：t, [x,y,vx,vy]，生成一个样品了

def rk4_solver(t, y, dt): ## 给定微分方程的 R-K4 算法求解器
    def f(y): # y[0]=x, y[1]=y, y[2]=vx, y[3]=vy
        Ex, Ey = calc_elec(t, pol)
        r3 = sqrt(y[0]**2 + y[1]**2 + 0.04)**3  ## 软化库伦势相应的 r^3
        return [y[2], y[3], -Ex-y[0]/r3, -Ey-y[1]/r3] # 根据微分方程定出 y 的导函数
    
    stop = False
    while True: # R-K4 的每次迭代
        if 2*T-t < dt: # 达到判停条件？
            dt = 2*T-t
            stop = True
        k1 = f(y)
        k2 = f([y[i] + dt/2*k1[i] for i in range(4)])
        k3 = f([y[i] + dt/2*k2[i] for i in range(4)])
        k4 = f([y[i] + dt*k3[i] for i in range(4)])
        for i in range(4):
            y[i] += dt/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        t += dt
        if stop == True:
            break
    return t, y
def coulomb_evol(t, par): ## 激光结束后在库伦场运动，格式 par=(x,y,vx,vy)
    r = sqrt(par[0]**2 + par[1]**2)  # 不必软化库伦势
    p2 = (par[2]**2 + par[3]**2) - 2/r
    if p2 <= 0:  # 能量 < 0 则舍去
        return None
    Lz = par[0] * par[3] - par[1] * par[2]  # 角动量 z 分量
    ax = par[3] * Lz - par[0]/r  # 隆格楞次矢量 x 分量
    ay = -par[2] * Lz - par[1]/r  # 隆格楞次矢量 y 分量
    denom = 1 + p2 * Lz**2
    return [(p2*(-Lz*ay) - sqrt(p2)*ax)/denom, (p2*(Lz*ax) - sqrt(p2)*ay)/denom]

if __name__ == '__main__':
    A0 = 1.325
    WEmax = 0.010623
    Ip = 0.5

    N = 1000
    ome = 0.057
    T = 2*pi/ome
    p_infx, p_infy = [], []
    Nloop = 0
    for i in range(N):
        t, par = gen_sample() # 先生成样品
        rk4_solver(t, par, dt=0.5) # R-K4演化到激光结束，步长取为0.5
        p_inf = coulomb_evol(t, par) # 计算无穷远动量
        if p_inf == None: # 能量为负?（没有返回值），舍去
            continue
        if i % 100 == 0:
            print('已经生成 %-6d 个样品啦'%i)
        p_infx.append(p_inf[0])
        p_infy.append(p_inf[1])

    #########################################################################
    ## 使用 matlibplot 作图
    from matplotlib import pyplot as plt
    ## 首先是 px-py 2D 图
    fig, ax = plt.subplots(figsize=(6.4,5))
    h = ax.hist2d(p_infx, p_infy, bins=[150, 150],range=[[-1.5,1.5],[-1.5,1.5]], cmap='jet')
    plt.colorbar(h[3], ax=ax)
    plt.xlabel(r'$p_x$')
    plt.ylabel(r'$p_y$')
    plt.savefig('Q2_%s_n%d_pxpy2d.pdf'%(pol, N))
    plt.show()
    ## px 直方图
    plt.hist(p_infx, bins=150, range=[-1.5,1.5], histtype='step', label=r'$p_x$ (%s)'%pol)
    plt.legend()
    plt.xlabel(r'$p_x$')
    plt.ylabel('Entries')
    plt.savefig('Q2_%s_n%d_px.pdf'%(pol, N))
    plt.show()
    ## py 直方图
    plt.hist(p_infy, bins=150, range=[-1.5,1.5], histtype='step', label=r'$p_y$ (%s)'%pol)
    plt.legend()
    plt.xlabel(r'$p_y$')
    plt.ylabel('Entries')
    plt.savefig('Q2_%s_n%d_py.pdf'%(pol, N))
    plt.show()