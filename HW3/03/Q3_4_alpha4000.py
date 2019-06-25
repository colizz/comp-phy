# -*- coding: utf-8 -*-
from numpy import sqrt, pi, cos, sin, exp, cosh, tanh  # 只导入常用函数啦

######################################################################################
### 以下是同样的函数
def rk4_solver(t, y, dt, fpu='0', saveEt=4, saveq=0): ## 给定微分方程的 R-K4 算法求解器
    
    ## 参数解释：基本参数：t, y(t), dt; fpu='0', 'alpha', 'beta'选择特定求解工具, saveEt: 过程中存入(前多少个)Et, saveq: 过程中存入(前多少个)q
    
    def f_0(y): # y = [q1, q2, ..., qn, p1, p2, ..., pn] 这是无高次项的微分方程
        dp = [0 for i in range(n)]
        for i in range(1, n-1):
            dp[i] = y[i-1] + y[i+1] - 2*y[i]
        dp[0] = y[1] - 2*y[0]
        dp[n-1] = y[n-2] - 2*y[n-1]
        return [y[i+n] for i in range(n)] + dp
    
    def f_alpha(y): # y = [q1, q2, ..., qn, p1, p2, ..., pn] 这是 alpha-FPU 的微分方程
        dp = [0 for i in range(n)]
        for i in range(1, n-1):
            dp[i] = (y[i-1] + y[i+1] - 2*y[i]) * (1 + alpha*(y[i-1] - y[i+1]))
        dp[0] = (y[1] - 2*y[0]) * (1 - alpha*y[1])
        dp[n-1] = (y[n-2] - 2*y[n-1]) * (1 + alpha*y[n-2])
        return [y[i+n] for i in range(n)] + dp
    
    def f_beta(y): # y = [q1, q2, ..., qn, p1, p2, ..., pn] 这是 beta-FPU 的微分方程
        dp = [0 for i in range(n)]
        for i in range(1, n-1):
            dp[i] = y[i-1] + y[i+1] - 2*y[i] + beta*(y[i-1]-y[i])**3 - beta*(y[i]-y[i+1])**3
        dp[0] = y[1] - 2*y[0] - beta*(y[0]-y[1])**3
        dp[n-1] = y[n-2] - 2*y[n-1] + beta*(y[n-2]-y[n-1])**3
        return [y[i+n] for i in range(n)] + dp
    
    if fpu == 'alpha': # 选择用哪个方程
        f = f_alpha
    elif fpu == 'beta':
        f = f_beta
    else:
        f = f_0
    stop = False
    while True: # R-K4 的每次迭代
        if Tmax-t < dt: # 判停条件
            dt = Tmax-t
            stop = True
        k1 = f(y)
        k2 = f([y[i] + dt/2*k1[i] for i in range(2*n)])
        k3 = f([y[i] + dt/2*k2[i] for i in range(2*n)])
        k4 = f([y[i] + dt*k3[i] for i in range(2*n)])
        for i in range(2*n):
            y[i] += dt/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        t += dt
        t_list.append(t)
        if saveEt > 0: # 若要过程中储存 Et
            Et.append(calc_energy(y, lead=saveEt))
        if saveq > 0:  # 若要过程中储存 q
            q_list.append(y[:saveq])
        if stop == True:
            break
    return t, y

def calc_energy(y, lead): # 通过当前 y=(q, p) 计算能量（要先变换到 Q, P）
    Q = [sum([Aij[(i+1)*(j+1)]*y[j] for j in range(n)]) for i in range(lead)]
    P = [sum([Aij[(i+1)*(j+1)]*y[j+n] for j in range(n)]) for i in range(lead)]
    return [n/(n+1) * (0.5*P[i]**2 + 0.5*(ome[i]*Q[i])**2) for i in range(lead)]

### 同样的函数结束
######################################################################################

## 开始第(4)小问（运行要花好长时间）
if __name__ == '__main__':
    alpha = 0.25
    n = 32
    Q1_0 = 4
    TN = 4000
    ome = [2*sin(k*pi/(2*n+2)) for k in range(1, n+1)] # 简振模
    Aij = [sqrt(2/n)*sin(k*pi/(n+1)) for k in range(n**2+2)] # 变换矩阵系数提前算出来
    Tmax = TN * 2*pi / ome[0]
    Q_0 = [Q1_0] + [0 for i in range(n-1)]
    t, y = 0, [0 for i in range(2*n)]
    for i in range(n):
        y[i] = sum([Aij[(i+1)*(j+1)]*Q_0[j] for j in range(n)])

    ## 初始化时间与能量 Et 列表（存入第 1 个能量）
    t_list, Et = [], []
    t_list.append(t)
    Et.append(calc_energy(y, lead=1))

    ## 计算！(选步长 dt=0.1)
    print('开始使用 R-K4 计算，预计耗时较长')
    t, y = rk4_solver(t, y, dt=0.1, fpu='alpha', saveEt=1)

    #########################################################################
    ## 使用 matlibplot 作图
    from matplotlib import pyplot as plt
    for i in range(1):
        plt.plot([t_list[i] /(2*pi/ome[0]) for i in range(len(t_list))], [Et[m][i] for m in range(len(Et))], label='%d-th'%(i+1))
    plt.xlim([0, TN])
    plt.legend()
    plt.xlabel(r'$T\,/\,(2\pi/\omega_1)$')
    plt.ylabel(r'Energy')
    plt.savefig('Q3_4_alpha4000.pdf')
    plt.show()