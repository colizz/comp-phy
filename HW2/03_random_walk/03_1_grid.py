# -*- coding: utf-8 -*-
from numpy import sin, cos, sqrt, pi  # 导入运算基本函数
import random  # Python生成随机数的工具

N = 50  # 步数
X = 0  # 初始位置
Y = 0
Xlist = [0]
Ylist = [0]
random.seed(0)  # 指定seed确保每次运行产生的随机数结果相同
for i in range(N):
    r = random.random()  # 产生[0,1)之间随机数
    if r <= 0.25:  # 区间四等分，分别对应X+1, X-1, Y+1, Y-1
        X += 1
    elif r <= 0.50:
        X -= 1
    elif r <= 0.75:
        Y += 1
    else:
        Y -= 1
    Xlist.append(X)
    Ylist.append(Y)
    print('第 %2d 步,  ( X, Y ) = (%3d, %3d ),  R = %f'%(i+1, X, Y, sqrt(X**2+Y**2)))

    
########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt  # 作图工具包
plt.plot(Xlist, Ylist, label='random walk simulation')
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.legend()
plt.grid()
plt.savefig('03_1.pdf')
plt.show()