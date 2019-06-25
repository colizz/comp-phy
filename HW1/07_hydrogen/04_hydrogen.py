# -*- coding: utf-8 -*-

from numpy import sin, cos, pi, sqrt, exp

def laguerreL(n, a, x):
    n0 = 1
    La, La1 = -x+a+1, -x+a+2 # 初始 n=1 时的 L_n^a(x), L_n^(a+1)(x)
    while n0 < n: # 不断迭代
        n0 += 1
        La, La1 = (-x*La1+(n0+a)*La)/n0, (-(x-n0)*La1+(n0+a)*La)/n0
    return La

n, l = 2, 1
Nr = 100000
dr = 60/Nr
r = [(i+1)*60/Nr for i in range(Nr)]
Psi = [sqrt(1/486)*exp(-r[i]/3)*2*r[i]/3*laguerreL(1,3,2*r[i]/3) for i in range(Nr)] # Ψ的格点函数值
HPsi = [-0.5*r[i]**2*(Psi[i+1]-2*Psi[i]+Psi[i-1])/dr**2 - r[i]*(Psi[i+1]-Psi[i-1])/(2*dr) + 
        Psi[i]*(1 - r[i]) for i in range(1, Nr-1)] # 算符作用在Ψ后的格点函数值
res = sum([Psi[i]*HPsi[i-1] for i in range(1, Nr-2)]) * dr
print("E_311 = {res}".format(res=res))