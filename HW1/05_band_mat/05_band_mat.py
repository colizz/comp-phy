# -*- coding: utf-8 -*-

def convert_diag(C, op='sym'): # 将矩形数 C 转化为带状矩阵，可选 上对角、下对角、对称矩阵三个选项
    A = [[0 for i in range(n)] for i in range(n)]
    for j in range(n):
        for i in range(j, min(j+m+1, n)):
            if op=='low' or op=='sym':
                A[i][j] = C[i+1][j-i+m+1]
            if op=='up' or op=='sym':
                A[j][i] = C[i+1][j-i+m+1]
    return matrix(A)

def initialize(): # 初始化矩形数 C 和向量 b
    C = [[0 for i in range(m+2)] for i in range(n+1)]
    for i in range(1, n+1):
        C[i][m+1] = 6.
    for i in range(2, n+1):
        C[i][m] = 4.
    for i in range(3, n+1):
        C[i][m-1] = 1.
    C[1][m+1] = 5.
    C[n][m+1] = 5.

    b = [120. for i in range(n+1)]
    b[1] = 60; b[n] = 60
    return C, b

n = int(input('输入矩阵的阶数 n = ')) # 矩阵阶数
m = 2  # 带宽 2m+1
C, b = initialize()
# 拷贝一份原始的C和b（注意!! 实际算法中不需要这一步哒，只是为了最后能检验下结果这里才这样做）
b0 = [bi for bi in b] 
C0 = [[Cij for Cij in Ci] for Ci in C]

for j in range(1, n+1): # j=1,..,n
    for i in range(j, min(j+m+1, n+1)): # i=j,j+1,..,j+m or n
        for k in range(max(1, i-m), j): # k=r,..,j-1
            C[i][j-i+m+1] -= C[i][k-i+m+1] * C[j][k-j+m+1] / C[k][m+1] # 也就是 a_ij = a_ij - Σ_r^(j-1) l_ik l_jk/ l_kk

for i in range(1, n+1):  # i=1,..,n
    for j in range(max(1, i-m), i): # j=r,..,i-1
        b[i] -= C[i][j-i+m+1] * b[j] / C[j][m+1] # 也就是 b'_i = b_i - Σ_r^(i-1) l_ij b'_j/ l_jj

x = [bi for bi in b]  # 初始时 x=b'
for i in range(n, 0, -1): # i=n,n-1,..,1
    for j in range(i+1, min(i+m+1, n+1)): # j=i+1,..,i+m or n
        x[i] -= C[j][i-j+m+1] * x[j]
    x[i] /= C[i][m+1]
print('x = ', x[1:])


########## 计算完啦，下面来画图验证一下 ##########
from matplotlib import pyplot as plt # 只是为了画图啦
plt.figure(figsize=(8,4))
plt.plot(range(n), x[1:], label=r'$x_i$')
plt.xlabel(r'$i$'); plt.ylabel(r'$x_i$')
plt.legend()
plt.savefig('06_n%d_xi.pdf'%n)
plt.show()

########## 最后验证一下结果的正确性 ##########
import numpy as np # 只用来计算向量的模长
print('\n计算完毕，下面验证结果是否正确')
print('‖ A x - b ‖ = ', np.linalg.norm([sum([(C0[i][j-i+m+1] if i>j else C0[j][i-j+m+1])*x[j] for j in range(max(1,i-m), min(i+m+1,n+1))])-b0[i] for i in range(1, n+1)]))