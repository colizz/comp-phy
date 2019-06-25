# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix
from func01 import mat_initialize, find_eigenvalues

from numpy import sin, cos, sqrt, exp, sign, pi


print('第(2)问：\n')
M = 50
eigen_list = []  #储存每个q下的一系列本征值
q_list = [i for i in range(0,21,1)]  # q的各参数点
for i, q in enumerate(q_list):
    H = mat_initialize(M, q)  # 先初始化原矩阵 H
    eigen_list.append(find_eigenvalues(H, need_eigenstates=False)[0:11]) # 只计算本征值，并储存头11个本征值
    print('正在对不同的 q 计算本征值：[%-3d / %-3d] 进行中...'%(i, len(q_list)))
print('计算完成，下面画图')


########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt
fig, ax1 = plt.subplots(figsize=(9,4.5)) # 调节画布大小等
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width* 0.7 , box.height])
for i in range(11):
    plt.plot(q_list, [eigen_list[iq][i] for iq in range(len(q_list))], label='the %d-th eigenvalue'%(i+1))
print('画图完成')
plt.xlim(0,20)
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.xlabel('q')
plt.ylabel(r'Eigenvalue $A$')
plt.savefig('01_2.pdf')
plt.show()