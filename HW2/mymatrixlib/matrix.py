# -*- coding: utf-8 -*-
from numpy import sqrt  

# 定义 矩阵类
class matrix(object):
    def __init__(self, element):
        self.nrow = len(element) # 行数
        self.ncolumn = len(element[0]) # 列数
        for i in range(self.nrow):
            if len(element[i]) != self.ncolumn:
                print('矩阵维数不对啊亲')
                del self
        self.elem = element[:] # 行向量, [[],[],[]]格式
        
    def __add__(self, other): # 重载+为矩阵加法
        if self.nrow != other.nrow or self.nrow != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        return matrix([[self.elem[i][j] + other.elem[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def __sub__(self, other): # 重载-为矩阵减法
        if self.nrow != other.nrow or self.nrow != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        return matrix([[self.elem[i][j] - other.elem[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def __mul__(self, other): # 重载*为矩阵乘法
        if self.ncolumn != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        C = matrix([[sum([self.elem[i][m] * other.elem[m][j] for m in range(self.ncolumn)]) for j in range(other.ncolumn)] for i in range(self.nrow)])
        # python就是妙，一行完成
        return C

    def mul(self, n): # 矩阵数乘
        return matrix([[n * self.elem[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def T(self): # 矩阵转置
        return matrix([[self.elem[i][j] for i in range(self.nrow)] for j in range(self.ncolumn)])
    
    def mod2(self): # 向量的模方
        if self.nrow != 1 and self.ncolumn != 1:
            print('向量不是 1 * n 格式的啊亲')
            return 'error'
        if self.nrow == 1:
            return sum([self.elem[0][i]**2 for i in range(self.ncolumn)])
        if self.ncolumn == 1:
            return sum([self.elem[i][0]**2 for i in range(self.nrow)])
    
    def mod(self): # 向量的模
        ans = self.mod2()
        if ans == 'error':
            return 'error'
        else:
            return sqrt(ans)

    def copy(self): # 矩阵的自身复制
        return matrix([[self.elem[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def paste_on(self, other, row_st=0, column_st=0): # 将给定矩阵粘贴在原矩阵指定区域上
        for i in range(other.nrow):
            for j in range(other.ncolumn):
                self.elem[row_st+i][column_st+j] = other.elem[i][j]
    
    def slice(self, row_st, row_end, column_st, column_end): # 取原矩阵的指定区域的子阵
        return matrix([[self.elem[i][j] for j in range(column_st,column_end)] for i in range(row_st,row_end)])
        
    def zeros(nrow, ncolumn): # 全零矩阵
        return matrix([[0. for j in range(ncolumn)] for i in range(nrow)])
    
    def identity(nrow): # 对角阵
        return matrix([[1. if i==j else 0. for j in range(nrow)] for i in range(nrow)])
    
    def print(self, n=3): # 按矩阵格式打印
        for i in range(self.nrow):
            print([round(self.elem[i][j], n) for j in range(self.ncolumn)])
        print('\n')