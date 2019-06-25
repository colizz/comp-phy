# -*- coding: utf-8 -*-
from numpy import sin, cos, sqrt, pi  # 导入运算基本函数
import random  # Python生成随机数的工具

random.seed(0)  # 指定seed确保每次运行产生的随机数结果相同
NMC = 20000  # 蒙卡运行20000次
N = 100  # 随机行走100步
Rlist = []
for iMC in range(NMC):
    X = 0.
    Y = 0.
    Rsingle = []
    for i in range(N):
        theta = random.random() * 2*pi  # 生成[0,2π)随机数
        X += cos(theta)
        Y += sin(theta)
        Rsingle.append(sqrt(X**2 + Y**2))  # 储存100步中每步的R_n
    Rlist.append(Rsingle)  # 将每次蒙卡的Rsingle储存进来
    if iMC % 2000 == 0:
        print('正在进行随机行走模拟：[%-5d / %-5d] 进行中...'%(iMC, NMC))

R = [sum([Rlist[i][j] for i in range(NMC)])/NMC for j in range(N)]  # 对10000次蒙卡下第n步的R_n取平均，即为期望E(R_n)的实验值


########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt  # 作图工具包
from scipy.optimize import curve_fit  # 拟合工具包
import numpy as np
def func(x, a):
    return a * sqrt(x)
a_fit, cov = curve_fit(func, np.arange(1,N+1), np.array(R))
plt.plot(range(1,N+1), R, '.', label='random walk simulation')
plt.plot(range(N+1), [a_fit[0] * sqrt(i) for i in range(N+1)], linewidth=0.6, label='fit curve')
plt.xlabel(r'Step $n$')
plt.ylabel(r'Average distance $R_n$')
plt.legend()
plt.text(0, 7.8, r'$R_n = %.5f\sqrt{n}$'%a_fit[0])
plt.ylim([0, 100])
plt.ylim([0, 10])
plt.savefig('03_2.pdf')
plt.show()