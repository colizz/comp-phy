# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix
from func01 import mat_initialize, find_eigenvalues

from numpy import sin, cos, sqrt, exp, sign, pi


def f(m, k, phi):  # 定义基函数 f_k(phi)
    return 1/(m*2*pi) * sum([exp(-1j*(n-(m-1)/2)*(2*pi*k/m))*exp(1j*(n-(m-1)/2)*phi) for n in range(m)]).real


print('第(3)问：大概要耗时 1.5min, 请稍等~ \n')
eigen_list_M5 = []
eigen_list_M40 = []
q_list = [iq*1. for iq in range(0,21,1)]
for iq, q in enumerate(q_list):
    print('M = 5, 40, 正在对不同的 q 计算“本征值”和“本征矢（花费时间较多）”：  [%-3d / %-3d] 进行中...'%(iq, len(q_list)))
    H = mat_initialize(5, q)  # M=5
    eigen_list_M5.append(find_eigenvalues(H, need_eigenstates=True)) # 将计算出的本征值和本征矢都存入eigen_list
    H = mat_initialize(40, q)  # M=40
    eigen_list_M40.append(find_eigenvalues(H, need_eigenstates=True)) # 将计算出的本征值和本征矢都存入eigen_list

# 下面从is_even列表中提取出偶宇称
eigen_even_list_M5 = []  # M=5 偶宇称态的本征值列表
eigen_even_list_M40 = []  # M=40 偶宇称态的本征值列表
for iq in range(len(q_list)):
    iqe = iq if iq!=0 else 1  # 因奇偶宇称的简并，需对 q=0 宇称态情形的特殊处理（详见“算法提纲”一节）
    eigens_M5 = []
    eigens_M40 = []
    for i, flag in enumerate(eigen_list_M5[iqe]['is_even']):
        if flag == True:  # 如果是奇宇称
            eigens_M5.append(eigen_list_M5[iq]['eigens'][i])
    for i, flag in enumerate(eigen_list_M40[iqe]['is_even']):
        if flag == True:
            eigens_M40.append(eigen_list_M40[iq]['eigens'][i])
    eigen_even_list_M5.append(eigens_M5)
    eigen_even_list_M40.append(eigens_M40)
print('计算完成，下面画图')

    
########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt
fig, ax1 = plt.subplots(figsize=(10,4.5)) # 调节画布大小等
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width* 0.65 , box.height])
color = ['cyan', 'orange', 'green', 'red', 'purple']
for i in range(5):
    plt.plot(q_list, [eigen_even_list_M5[iq][i] for iq in range(len(q_list))], color=color[i], linestyle='--', label=r'the %d-th $P$-even eigenvalue (M=5)'%(i+1))
    plt.plot(q_list, [eigen_even_list_M40[iq][i] for iq in range(len(q_list))], color=color[i], label='the %d-th $P$-even eigenvalue (M=40)'%(i+1))
print('画图完成')
plt.xlim(0,20.)
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.xlabel(r'$q$')
plt.ylabel(r'Eigenvalue $A$')
plt.savefig('01_3.pdf')
plt.show()