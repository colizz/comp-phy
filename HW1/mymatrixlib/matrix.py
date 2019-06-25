# -*- coding: utf-8 -*-
import numpy as np  # 只用numpy里的sprt，sin，cos简单函数啊喵

# 定义 矩阵类
class matrix(object):
    def __init__(self, element):
        self.nrow = len(element) # 行数
        self.ncolumn = len(element[0]) # 列数
        for i in range(self.nrow):
            if len(element[i]) != self.ncolumn:
                print('矩阵维数不对啊亲')
                del self
        self.row = element[:] # 行向量, [[],[],[]]格式
        self.column = [[element[i][j] for i in range(self.nrow)] for j in range(self.ncolumn)] # 列向量, [[],[],[]]格式
        
    def __add__(self, other): # 矩阵加法
        if self.nrow != other.nrow or self.nrow != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        return matrix([[self.row[i][j] + other.row[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def __sub__(self, other): # 矩阵减法
        if self.nrow != other.nrow or self.nrow != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        return matrix([[self.row[i][j] - other.row[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def __mul__(self, other): # 矩阵乘法
        if self.ncolumn != other.nrow:
            print('行列不匹配啊亲: A: %d * %d, B: %d * %d'%(self.nrow, self.ncolumn, other.nrow, other.ncolumn))
            return 'error'
        C = matrix([[sum([self.row[i][m] * other.row[m][j] for m in range(self.ncolumn)]) for j in range(other.ncolumn)] for i in range(self.nrow)])
        # python就是妙，一行完成
        return C

    def mul(self, n):
        return matrix([[n * self.row[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
    
    def T(self): # 矩阵转置
        return matrix(self.column)
    
    def mod2(self):
        if self.nrow != 1 and self.ncolumn != 1:
            print('向量不是 1 * n 格式的啊亲')
            return 'error'
        if self.nrow == 1:
            return sum([self.row[0][i]**2 for i in range(self.ncolumn)])
        if self.ncolumn == 1:
            return sum([self.row[i][0]**2 for i in range(self.nrow)])
    
    def mod(self):
        ans = self.mod2()
        if ans == 'error':
            return 'error'
        else:
            return np.sqrt(ans)

    def copy(self): # 矩阵复制
        return matrix([[self.row[i][j] for j in range(self.ncolumn)] for i in range(self.nrow)])
