# -*- coding: utf-8 -*-
import numpy as np

def ieee2dec(a_bin_str):
    # 错误提示
    if len(a_bin_str) != 64:
        print('输入的不是 64 位啊亲！')
        return 'error'
    a_bin = [0 for i in range(64)]
    for i in range(64):
        a_bin[i] = int(a_bin_str[i])
        if a_bin[i]!= 0 and a_bin[i]!= 1:
            print('只能输入 0 或 1 啊亲！')
            return 'error'

    # 确定幂次
    exp = -1023
    for i in range(1,12):
        if a_bin[i] == 1:
            exp += 2**(11-i)

    # 将每一个小数位转化为十进制相加
    a_dec = 0
    for i in range(12,64):
        if a_bin[i] == 1:
            a_dec += 2**(12+exp-i)

    # 考察正负，得到结果
    minusQ = -1 if a_bin[0]==1 else 1
    a_dec *= minusQ
    return a_dec

a_bin_str = input('请输入任意 64 位双精度浮点编码: ')
print('转化为十进制数的结果为: ', ieee2dec(a_bin_str))