# -*- coding: utf-8 -*-

def laguerreL(n, a, x):
    n0 = 1
    La, La1 = -x+a+1, -x+a+2 # 初始 n=1 时的 L_n^a(x), L_n^(a+1)(x)
    while n0 < n: # 不断迭代
        n0 += 1
        La, La1 = (-x*La1+(n0+a)*La)/n0, (-(x-n0)*La1+(n0+a)*La)/n0
    return La
    
for n in [3, 10, 30]:
    for a in [2, 20, 40]:
        for x in [0.001, 1, 100]:
            print("L_{n}^{a}({x}) = {ans}".format(n=n, a=a, x=x, ans=laguerreL(n,a,x)))