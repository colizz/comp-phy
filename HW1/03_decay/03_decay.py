# -*- coding: utf-8 -*-
def func(y, x): 
    # 指定 f(x,y)形式，当然本例中f与x (也即时间t)无关
    return [-y[0]/tauA, y[0]/tauA - y[1]/tauB]

def next_step(func, y, x, tstep):
    # y的下一步演化值
    return [y[i] + func(y, x)[i]*tstep for i in range(2)]

def calculator(tauA, tauB, tstep, T): # 封装求解过程
    y = [1., 1.] # 初始函数值(NA, NB)
    Lt = [i*tstep for i in range(int(T/tstep) + 1)] # 等步长的时间演化向量
    Ly = []
    Ly.append(y)
    for t in range(len(Lt)-1):
        y = next_step(func, y, t, tstep) # y在下一时间的演化结果
        Ly.append(y)
    return Lt, Ly

# 下面开始就具体参数值计算微分方程数值解
Ly_dict = {}
Lt_dict = {}
tauA = 1.
tauB = 1.
for tauB in [0.1, 1., 10.]: # 改变3种不同的tauB
    Lt, Ly_dict['tauB=%.1f'%tauB] = calculator(tauA=1., tauB=tauB, tstep=0.01, T=10.)
    print('tauB=%.1f 计算完毕'%tauB)
for tstep in [0.2, 0.1, 0.05]: # 改变3种不同的步长tstep
    Lt_dict['tstep=%.2f'%tstep], Ly_dict['tstep=%.2f'%tstep] = calculator(tauA=1., tauB=10., tstep=tstep, T=20.)
    print('tstep=%.2f 计算完毕'%tstep)


    

########## 计算完啦，下面来画图验证一下 ##########
from matplotlib import pyplot as plt # 只是为了画图啦
import numpy as np # 为了用常数e
plt.figure(figsize=(8,4))
plt.plot(Lt, [Ly_dict['tauB=0.1'][i][0] for i in range(len(Lt))], label=r'$N_A (t)$', linestyle='--')
for tauB in [0.1, 1., 10.]:
    plt.plot(Lt, [Ly_dict['tauB=%.1f'%tauB][i][1] for i in range(len(Lt))], label=r'$N_B (t)$, $\tau_B$=%.1f'%tauB)
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$N(t)$')
plt.legend(loc=1)
plt.savefig('03_tauB.pdf')
plt.show()

plt.figure(figsize=(8,4))
for tstep in [0.2, 0.1, 0.05]: # 改变3种不同的步长tstep
    plt.plot(Lt_dict['tstep=%.2f'%tstep], [Ly_dict['tstep=%.2f'%tstep][i][1] for i in range(len(Lt_dict['tstep=%.2f'%tstep]))], 
             label=r'$N_B (t)$, $\Delta t$=%.2f'%tstep, linewidth=0.6)
def func_realB(t):
    return np.exp(-t/10)+10./(1-10)*(np.exp(-t/1) - np.exp(-t/10))
plt.plot(np.linspace(0.,20.,1000), func_realB(np.linspace(0.,20.,1000)), label=r'$N_B (t)$, accurate', linewidth=0.6, linestyle='--')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$N(t)$')
plt.legend()
plt.savefig('03_tstep.pdf')
plt.show()
