# -*- coding: utf-8 -*-
import numpy as np

def dec2ieee(a_dec):
    a_dec = float(a_dec)
    
    # 讨论特殊a=0情况
    if a_dec == 0:
        return '0'*64
    
    # 确定正负，幂次
    a = a_dec if a_dec>0 else -a_dec
    n = np.int(np.log2(a))
    a_bin = [0 for i in range(64)]
    a_bin[0] = 0 if a_dec>0 else 1

    # 确定13-64位，即二进制下的科学计数法
    digit = 12
    for i in range(n,n-52,-1):
        if a >= 2**i:
            a -= 2**i
            a_bin[digit] = 1
        digit = digit+1

    # 确定2-12位，即幂次
    exp = 1023 + n
    if exp > 2047 or exp < 0:
        print("出界了您呐!")
    for i in range(1,12):
        if exp >= 2**(11-i):
            exp -= 2**(11-i)
            a_bin[i] = 1

    # 转化为字符串输出
    a_bin_str = ''
    for i in range(64):
        a_bin_str += '%s'%a_bin[i]
    
    return a_bin_str

a_dec = input('请输入任意十进制小数: ')
print('双精度浮点数结果为: ', dec2ieee(a_dec))