# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix

import numpy as np  # 只用numpy里的sprt，sin，cos简单函数啊喵
def specify_wfunc(wfunc='sqrt(x)', node=2): 
    # 通过给定的权函数 ω(x) 与 node 的值，直接给出 ∫x^p ρ(x) dx 的值作为求解的 input
    node += 1
    if wfunc == 'sqrt(x)':
        b = [2/(2*i +1) for i in range(2*node-2)]
    elif wfunc == '1+x^2':
        b = [2/(2*i + 1) + 2/(2*i + 3) for i in range(2*node-2)]
    return b

def poly_eval(L, x): 
    # 求多项式值
    y = 0
    for i in range(len(L)):
        y = y + L[i] * x ** i
    return y 
def ortho_poly(node, b): 
    # 利用格拉姆-施密特正交化求解多项式值
    ortho_coef = [0 for i in range(node-1)]
    ortho_coef.append(1)
    for i in range(node-1):
        for j in range(i+1,node):
            a = 0 # 这是要计算的分母<x^i, φ_j(x)>
            c = 0 # 这是要计算的分子<φ_j(x), φ_j(x)>
            for l in range(j): 
                a = a + b[l+node-1] * ortho_poly(j, b)[l]
            for l in range(j):
                for m in range(l+1):
                    c = c + b[l] * ortho_poly(j, b)[m] * ortho_poly(j, b)[l-m]
            for l in range(j, 2*j-1):
                for m in range(l+1-j, j):
                    c = c + b[l] * ortho_poly(j, b)[m] * ortho_poly(j, b)[l-m]
            ortho_coef[i] = ortho_coef[i] - a / c * ortho_poly(j, b)[i] # 正交多项式的各项系数
    return ortho_coef

def newton(L): 
    # 牛顿法求根，区间[0,1]，由于只需解决题目只对根个数<=3有效
    nzeros = len(L) - 1
    Lx = []
    if nzeros >= 1 : # 找第一个根
        x1 = -1
        x2 = -0.85
        while abs(poly_eval(L, x2)) > 0.00000001 or abs(x2-x1) > 0.00000001: # 是否足够接近零点
            s = x2 - poly_eval(L, x2) / (poly_eval(L, x2) - poly_eval(L, x1)) * (x2-x1)
            x1 = x2
            x2 = s
        Lx.append(x2)
    if nzeros >= 2: # 找第二个根
        x1 = 1
        x2 = 0.85
        while abs(poly_eval(L, x2)) > 0.00000001 or abs(x2-x1) > 0.00000001: # 是否足够接近零点
            s = x2 - poly_eval(L, x2) / (poly_eval(L, x2) - poly_eval(L, x1)) * (x2-x1)
            x1 = x2
            x2 = s
        Lx.append(x2)
    if nzeros == 3: # 第三个根用韦达定理找啦，更多的就不再本例里实现了。。
        Lx.insert(1, -L[nzeros-1] / L[nzeros]-Lx[0]-Lx[-1])
    
    return Lx

def func_max(L, func): # 返回最大值点
    maxvalue = func(L[0])
    s = 0
    for i in range(len(L)):
        if func(L[i]) > maxvalue:
            maxvalue = func(L[i])
            s = i
    maxvalue = L[s]
    return maxvalue

def up_solve(A,b): # 解上三角方程组
    x = []
    for i in range(A.nrow):
        x.append(0)
    for i in range(A.nrow-1, -1, -1):
        x[i] = b[i]
        for j in range(A.nrow-1, i, -1):
            x[i] = x[i] - A.row[i][j] * x[j]
        x[i] = x[i] / A.row[i][i]
    return x

def low_solve(A,b): # 解下三角方程组
    x = []
    for i in range(A.nrow):
        x.append(0)
    for i in range(A.nrow):
        x[i] = b[i]
        for j in range(i):
            x[i] = x[i] - A.row[i][j] * x[j]
        x[i] = x[i] / A.row[i][i]
    return x

def LU_decomp(M, B): # 线性方程组的列选主元LU分解解法
    A = M.copy()
    p = range(A.nrow)
    for k in range(A.nrow-1):
        # 列主元
        L = A.column[k][k:]
        s = k + L.index(func_max(L, abs))
        if s != k:
            exc = p[s]
            p[s] = p[k]
            p[k] = exc
            exc = A.row[s][:]
            A.row[s] = A.row[k][:]
            A.row[k] = exc[:]
        # LU分解
        for i in range(k+1, A.nrow):
            A.row[i][k] = A.row[i][k] / A.row[k][k]
            for j in range(k+1, A.nrow):
                A.row[i][j] = A.row[i][j] - A.row[i][k] * A.row[k][j]
    b0 = []
    for i in range(A.nrow):
        b0.append(B[p[i]])
    A1 = matrix.copy(A)
    for i in range(A1.nrow):
        A1.row[i][i] = 1
    # 解方程
    b1 = low_solve(A1, b0) 
    x = up_solve(A, b1)
    return x

def cal_gauss_int(wfunc, node=2):
    node += 1
    if wfunc == 'sqrt(x)':
        b = [2/(2*i + 3) for i in range(2*node-2)]
    elif wfunc == '1+x^2':
        b = [8/3]
        for i in range(node-1):
            b.append(0)
            b.append(2/(2*i + 3) + 2/(2*i + 5))

    opoly = ortho_poly(node, b) # 求正交多项式
    Lx = newton(opoly) # 正交多项式求根，得到节点表Lx
    MA = [[Lx[j]**i for j in range(len(Lx))] for i in range(len(Lx))]  # 利用节点表制成矩阵
    bA = b[:len(Lx)]
    LA = LU_decomp (matrix(MA), bA) # 解出求积系数表
    print('权函数: ', wfunc, ';  N =', node-1)
    print('节点坐标 x_i: ', Lx)
    print('对应系数 A_i: ', LA)
    print('\n')
    return Lx, LA

Lx_1, LA_1 = cal_gauss_int(wfunc='sqrt(x)', node=2)
Lx_2, LA_2 = cal_gauss_int(wfunc='sqrt(x)', node=3)
Lx_3, LA_3 = cal_gauss_int(wfunc='1+x^2', node=2)

from scipy import integrate # 使用scipy只是为了检验下答案是否正确
def myfunc(x):
    return np.sin(x)**3
def myfunc_w(x):
    return np.sqrt(x) * myfunc(x)
print('对 权函数:  sqrt(x) ;  N = 3 的验证')
print('我的函数: sin(x)^3 ')
true_value  = integrate.quad(myfunc_w, 0, 1)[0]
gauss_value = sum([LA_2[i] * myfunc(Lx_2[i]) for i in range(len(Lx_2))])
print('--> scipy积分值:   ', true_value)
print('    Gauss求积结果: ', gauss_value)
print('    相对误差:  ', (gauss_value - true_value)/true_value*100, "%")