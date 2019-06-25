# -*- coding: utf-8 -*-
from cmath import sin, cos, sqrt, exp, pi  # 复数运算基本函数

def sort_real(array):  # 对数列的实部排序（冒泡）
    def swap(a, b):
        return b, a
    for i in range(len(array)-1, 0, -1):
        for j in range(i):
            if array[j].real > array[j+1].real:
                array[j], array[j+1] = swap(array[j], array[j+1])
    return array
    
def dSfunc(t): # S'函数的解析形式
    return 0.5*(E0/ome * sin(ome*t) * sin(ome*t/4)**2 + pz)**2 + Ip

def ddSfunc(t): # S''函数的解析形式: 手动求导得到
    return E0*sin(ome*t/4) * (E0/ome * sin(ome*t) * sin(ome*t/4)**2 + pz) \
            * (cos(ome*t)*sin(ome*t/4) + 0.5*sin(ome*t)*cos(ome*t/4))

def Sfunc(t): # S函数的解析形式: 手动积分得到
    return 1/(ome**3)* (3/32*E0**2 *ome*t + 1/2*(2*Ip+pz**2)*ome**3 *t \
            + 1/6*pz*E0*ome*(3*cos(ome*t/2) - 3*cos(ome*t) + cos(3*ome*t/2)) 
            + E0**2 *(-1/4*sin(ome*t/2) + 1/64*sin(ome*t) + 1/24*sin(3*ome*t/2) \
            - 3/64*sin(2*ome*t) + 1/40*sin(5*ome*t/2) - 1/192*sin(3*ome*t))) \
            - E0*pz/(6*ome**2)

def muller_findroot(x0, x1, x2): # Muller求根法，与讲义的格式相同
    it = 2
    while it <= 100:
        h0 = x0 - x2  # 以下计算二次函数相关变量
        h1 = x1 - x2
        del0 = (dSfunc(x0)-dSfunc(x2)) / h0
        del1 = (dSfunc(x1)-dSfunc(x2)) / h1
        d = (del0 - del1) / (h0 - h1)
        b = del0 - h0 * d
        D = sqrt(b**2 - 4*dSfunc(x2)*d)
        x3 = x2 - 2*dSfunc(x2)/(b+D) if abs(b-D) < abs(b+D) else x2 - 2*dSfunc(x2)/(b-D)  # 距离零点较近的根
        if abs(x3-x2) < 1e-5: # 精度够了，返回零点
            return x3
        x0 = x1  # 替换初始值
        x1 = x2
        x2 = x3
        it += 1
    else:
        return 'null'  # 循环大于100次，没找到根

#  下面开始
print('先计算，再输出各小问结果')
ome = 45.5633525316 / 3200  # 题中的常数
E0 = sqrt(5 / 3509.4448314)
Ip = 13.6 / 27.2
tmax = 4*pi/ome
pzstep = 0.01
pzlist = [0.01+pzstep*i for i in range(int(2./pzstep))]  # pz的参数点

## 下面用鞍点法（SPM）求解跃迁振幅
Mp02list_SPM = []  # 储存求得的跃迁振幅模方
N = 201 # 找 N 次以确定 6 个不同的根
ipz = 1  # pz的计数
for pz in pzlist:  # 遍历每个pz参数点
    solu = []
    for i in range(N-1):
        try:    
            ts = muller_findroot(tmax/N*i, tmax/N*(i+1), tmax/N*(i+2))  # Muller法求根！
        except:  # 如果求根过程中报错（比如除以0）直接跳过本次求根
            continue
        if ts=='null' or ts.real < 0 or ts.real > tmax or abs(dSfunc(ts))>1e-2:  # 没找到根，根出界了，或找到了错根，都直接跳过本次求根
            continue

        # 运行到此处说明找到了真实的根
        ts = ts.conjugate() if ts.imag < 0 else ts  # 如果虚部 <0 则取共轭，共轭也是真实的根
        for s in solu:
            if abs(ts.real-s.real)<1e-4 : # 和已经存入的物理的根相同吗？相同说明已经找到过它了，不存入solu
                break
        else:
            solu.append(ts) # 没有找到过，将其存入solu
    
    if len(solu) != 6: # 是不是找到了 6 个根？如果不是说明出错了
        print(pz, '不是 6 个根!')
    if pz == 1.: # 刚好把第一问结果存一下，一会输出
        solu_pz1 = sort_real(solu)
        
    Mp0 = 0  # 下面代入鞍点求解的公式
    for tsa in solu:
        Mp0 += exp(1j*Sfunc(tsa)) /ddSfunc(tsa)  # 依次求和 exp(iS(t))/S''(t)
    Mp0 *= -(2*Ip)**(5/4) / sqrt(2)
    Mp02list_SPM.append(abs(Mp0)**2)  # 求解完毕，将振幅模方存入Mp02list_SPM
    
    if ipz % 20 == 0:
        print('正在鞍点法（SPM）求解： [%-3d / %-3d]  正在进行中...'%(ipz, int(2./pzstep)))
    ipz += 1
    
## 下面直接积分（DI）求解跃迁振幅
pzstep = 0.01
pzlist = [0.01+pzstep*i for i in range(int(2./pzstep))]
Mp02list_DI = []
ipz = 1
for pz in pzlist:
    tstep = 0.1
    I = 0
    ti = tstep/2
    while ti < tmax:
        qz = pz + E0/ome * sin(ome*ti) * sin(ome*ti/4)**2
        I += qz*E0*(cos(ome*ti)*sin(ome*ti/4)**2 + 1/4*sin(ome*ti)*sin(ome*ti/2)) \
              /(pi*(qz**2 + 2*Ip)**3) * exp(1j*Sfunc(ti)) * tstep
        ti += tstep
    I *= 2**(7/2) * (2*Ip)**(5/4)
    Mp02list_DI.append(abs(I)**2)
    
    if ipz%20==0:
        print('正在直接积分法（DI）求解：  [%-3d / %-3d]  正在进行中...'%(ipz, int(2./pzstep)))
    ipz += 1

print('计算完毕，下面输出各小问结果\n--------------------------------------------\n')
print('第(1)问\np_z = 1.0 时，鞍点方程的 6 个解为')
for z in solu_pz1:
    print(z)
print('\n第(2), (3)问： 参见两个图: 两种方法的曲线比较 & 整体因子曲线')   
########################## 使用 matplotlib 画图 ##########################
from matplotlib import pyplot as plt
plt.plot([pz**2/2 for pz in pzlist], Mp02list_SPM, label='SPM method') # 画 SPM 振幅模方曲线
plt.plot([pz**2/2 for pz in pzlist], Mp02list_DI, label='DI method')  #  画 DI 振幅模方曲线
plt.xlabel(r'$E_k$')
plt.ylabel(r'$|M_\mathbf{p}^0|^2$')
plt.legend()
plt.xlim(0., 2.)
plt.savefig('03_SPM_DI.pdf')
plt.show()

# 再画下整体因子曲线看看
plt.figure()
plt.plot([pz**2/2 for pz in pzlist], [Mp02list_DI[i]/Mp02list_SPM[i] for i in range(len(pzlist))], label='factor (DI/SPM)')
plt.xlabel(r'$E_k$')
plt.ylabel(r'$|M_\mathbf{p}^0|^2_{\mathrm{DI}}\;/\;|M_\mathbf{p}^0|^2_{\mathrm{SPM}}$')
plt.legend()
plt.xlim(0., 2.)
plt.savefig('03_factor.pdf')
plt.show()