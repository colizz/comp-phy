# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix

def initialize():
    # 初始化矩阵A
    A = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        A[i][n-1-i] = 0.5
    for i in range(n-1):
        A[i][i] = 3.
        A[i][i+1] = -1.
        A[i+1][i] = -1.
    A[n-1][n-1] = 3.
    
    # 初始化向量b
    b = [[1.5] for i in range(n)]
    b[0][0] = 2.5
    b[n-1][0] = 2.5
    b[int(n/2)-1][0] = 1.0
    b[int(n/2)][0] = 1.0
    return matrix(A), matrix(b)

def conj_grad(A, b, x):
    r = b - A*x
    p = r
    k = 0
    while r.mod() > 0.000001:
        # 优化的共轭梯度法
        alpha = r.mod2() / (p.T()*A*p).row[0][0]
        x = x + p.mul(alpha)
        r1 = r
        r = r - (A*p).mul(alpha)
        beta = r.mod2() / r1.mod2()
        p = r + p.mul(beta)
        k += 1
    return x, k

n = int(input('输入矩阵的阶数 n = ')) # 矩阵阶数
A, b = initialize()
x0 = matrix([[5] for i in range(n)]) # 给定x的初始值，其实具有任意性
x, k = conj_grad(A, b, x0)
print('第 %d 步后得到结果:\nx ='%k, x.column[0])
    