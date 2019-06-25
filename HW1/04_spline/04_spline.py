# -*- coding: utf-8 -*-

def thomas(a, b, c, f): 
    # Thomas算法解三对角矩阵方程
    # a: -1对角线向量补首位, b: 主对角线向量, c: +1对角线向量补末位, f: 结果向量
    m = [a[i+1] / b[i] for i in range(len(f)-1)]
    # 追
    for i in range(len(f)-1): 
        b[i+1] = b[i+1] - m[i]*c[i]
        f[i+1] = f[i+1] - m[i]*f[i]
    x = [f[-1]/b[-1]]
    # 赶
    for i in range(len(f)-1, -1, -1):
        x.insert(0, (f[i] - c[i]*x[0]) / b[i])
    return x

def spline(x, y, Lx): # 计算插值函数在采样Lx上的值Ly
    # 计算N-2阶三对角矩阵
    h = [x[i+1] - x[i] for i in range(len(x)-1)]
    sh = [h[i+1] + h[i] for i in range(len(x)-2)]
    mu = [h[i]/sh[i] for i in range(len(x)-2)]
    lam = [1 - mu[i] for i in range(len(x)-2)]
    d = [6 * (y[i]/h[i]/sh[i] + y[i+2]/h[i+1]/sh[i] - y[i+1]/h[i]/h[i+1]) for i in range(len(x)-2)]
    b = [2 for i in range(len(x)-2)]
    
    # 追赶法解三对角方程组
    M = thomas(mu, b, lam, d)
    M.insert(0, 0)
    M.append(0)
    
    # 以指定步长遍历
    Ly = []
    hsqure6 = [h[i]**2 / 6 for i in range(len(h))]
    for i in range(len(Lx)):
        ii = 0
        for j in range(1, len(x)): # 找到相应分段
            if Lx[i] <= x[j]:
                ii = j - 1
                break
        h1 = (x[ii+1] - Lx[i]) / h[ii]
        h2 = (Lx[i] - x[ii]) / h[ii]
        # 计算插值函数值
        yy = M[ii] * h1**3 * hsqure6[ii] + M[ii+1] * h2**3 * hsqure6[ii] + (y[ii] - M[ii] * hsqure6[ii]) * h1 + (y[ii+1] - M[ii+1] * hsqure6[ii]) * h2 
        Ly.append(yy)
    return Ly

# 题中所给参数
x = [0., 3., 5., 7., 9., 11., 12., 13., 14., 15.]
y = [0, 1.2, 1.7, 2.0, 2.1, 2.0, 1.8, 1.2, 1.0, 1.6]
Lx = [i*0.1 for i in range(int(15/0.1))]
Ly = spline(x, y, Lx)
# 逐点打印
print('样条插值函数所得结果：')
for i in range(len(Lx)):
    print('x = %4.1f         y = %2.10f'%(Lx[i], Ly[i]))



########## 计算完啦，下面来画图验证一下 ##########
from matplotlib import pyplot as plt # 只是为了画图啦
from scipy.interpolate import interp1d # 用python自带spline来检验是否正确
import numpy as np # 只是用来画图啦
plt.figure(figsize=(8,4))
plt.plot(x, y, 'x', label='data')
# 我的spline
plt.plot(Lx, Ly, label='my spline', linewidth=0.6)
# python自带spline
plt.plot(np.linspace(0,15,500), interp1d(x, y, kind='cubic')(np.linspace(0,15,500)), '--', label='python\'s spline', linewidth=0.6)
plt.xlabel(r'$x$'); plt.ylabel(r'$y$')
plt.legend()
plt.savefig('04_spline_result.pdf')
plt.show()