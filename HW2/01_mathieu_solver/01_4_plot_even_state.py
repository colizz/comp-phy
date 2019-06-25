# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix
from func01 import mat_initialize, find_eigenvalues

from numpy import sin, cos, sqrt, exp, sign, pi


def f(m, k, phi):  # 定义基函数 f_k(phi)
    return 1/(m*2*pi) * sum([exp(-1j*(n-(m-1)/2)*(2*pi*k/m))*exp(1j*(n-(m-1)/2)*phi) for n in range(m)]).real


print('第(4)问：大概耗时 ~10s \n')
M = 50
q = 10.
m = 2*M + 1
H = mat_initialize(M, q)  # M=5 初始化矩阵
eigendic_M50 = find_eigenvalues(H, need_eigenstates=True) # 将计算出的本征值和本征矢都存入eigen_list
print('计算完成，下面画图')


########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt
fig, ax1 = plt.subplots(figsize=(10,4.5)) # 调节画布大小等
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width* 0.65 , box.height])
Nphi = 40  # 40个 φ 的格点
phi_list = [i*0.5*pi/Nphi for i in range(0, Nphi+1)]  # φ 区间[0,π/2]
linestyle = ['-', '--', '--', '-.', '-.', ':']
iplot = 0
for i, flag in enumerate(eigendic_M50['is_even']):
    if flag == True:  # 如果是偶宇称解，画图！
        plt.plot([phi*180/pi for phi in phi_list], [-sum([eigendic_M50['Q'].elem[n][i]*f(m, n+1, phi) for n in range(m)]) for phi in phi_list],\
                 linestyle=linestyle[iplot], label='the %d-th $P$-even eigenstate (M=50)'%(iplot+1))
        iplot += 1
        if iplot == 6:  # 画满 6 个为止
            break
print('画图完成')
plt.xlim(0,90.)
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.xlabel(r'$\phi$ [degree]')
plt.ylabel(r'$P$-even eigenstate $\Psi(\phi)$')
plt.savefig('01_4.pdf')
plt.show()
