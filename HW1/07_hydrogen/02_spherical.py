# -*- coding: utf-8 -*-

from numpy import sin, cos, pi, sqrt, log10, sign # 导入函数
from cmath import exp # e指数需要从复数中导入

class strongdouble(object): # 编写“超精度”浮点运算，sig存significand；deci存幂次部分
    def normalize(self): # 标准化
        digit = int(log10(abs(self.sig)))
        self.sig /= 10**digit
        self.deci += digit
    def __init__(self, sig, deci):
        self.sig = sig
        self.deci = deci
        self.normalize()
    def __add__(self, other): # 加法
        dist = self.deci - other.deci
        if dist < 0:
            self, other = other, self
        self.sig += other.sig*10**(-abs(dist))
        return self
    def times(self, num): # 乘以一个普通double
        self.sig *= num
        self.normalize()
        return self
    def sprint(self): # 打印成字符串
        return "%fE%d"%(self.sig, self.deci)
    
def sphericalharmonicY(l, m, theta, phi): # 输入参数，输出球谐函数“超精度”字符串结果
    x = cos(theta)
    Pm, Pm1 = strongdouble(1, 0), strongdouble(1, -10000) # 作为每一步 P_l-1^m, P_l, m 进行递归，初始为 1; 0
    l0 = m
    for mul in [2*i+1 for i in range(m)]:
        Pm.times(-mul*sin(theta))
    while l0 < l: # 递归
        tmp = strongdouble(Pm.sig, Pm.deci)
        Pm, Pm1 = Pm.times(x*(2*l0+1)/(l0+1-m))+Pm1.times(-(l0+m)/(l0+1-m)), tmp
        l0 += 1
    for mul in [i for i in range(l-m+1, l+m+1)]:
        Pm.times(1/sqrt(mul))
    Pm.times(sqrt((2*l+1)/(4*pi)))
    PmRe = strongdouble(Pm.sig*cos(m*phi), Pm.deci)
    PmImabs = strongdouble(abs(Pm.sig*sin(m*phi)), Pm.deci)
    PmImsign = sign(Pm.sig*sin(m*phi))
    return "{re} {sign} {im}j".format(re=PmRe.sprint(), im=PmImabs.sprint(), sign='+' if PmImsign>=0 else '-')

for l in [100, 500, 1000]:
    for m in [1, int(l/100), int(l/10), l-1]:
        print("Y_{l},{m} (theta=pi/1000,phi=pi/5) = {ans}".format(l=l, m=m, ans=sphericalharmonicY(l, m, pi/1000, pi/5)))
        print("Y_{l},{m} (theta=3pi/10,phi=pi/5) = {ans}".format(l=l, m=m, ans=sphericalharmonicY(l, m, 3*pi/10, pi/5)))
        print("Y_{l},{m} (theta=501pi/1000,phi=pi/5) = {ans}".format(l=l, m=m, ans=sphericalharmonicY(l, m, 501*pi/1000, pi/5)))