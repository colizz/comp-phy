# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix
 
from numpy import sin, cos, sqrt, exp, sign, pi

# 通过指定参数 M，q 初始化哈密顿量矩阵
def mat_initialize(M, q): 
    # 初始化矩阵: 对角元为 M(M+1)/3 + V_ii, 非对角元为较复杂的 T_ij
    return matrix([[M*(M+1)/3 + 2*q*cos(4*pi*(j+1)/(2*M+1)) if j==k else 
                 (-1)**(j-k) * cos(pi*(j-k)/(2*M+1)) / (2 * sin(pi*(j-k)/(2*M+1))**2 )
                 for k in range(0, 2*M+1)] for j in range(0, 2*M+1)])


# 输入矩阵，求解本征值、本征矢并判断宇称
def find_eigenvalues(H, need_eigenstates=False): # H 为输入矩阵，need_eigenstates为是否需要输出本征矢，如不需要可简化运算。
    
    # 有关 Givens 矩阵乘法的简易乘法运算
    def mat_update_from_given(Mat, ndim, j, angle, pos='front', bandshape=True): 
        # angle: 旋转的角度，格式为[sin,cos]；pos: Givens矩阵是前/后乘在某矩阵上；bandshape表示矩阵是否为三对角的（若是，可进一步简化）
        if bandshape==True:  # 是三对角阵
            nstart = max(0, j-1)  # 只需截取长度1*4的向量即可
            nend = min(ndim, j+3)
        else:  # 不是，需截取整行/整列向量
            nstart = 0
            nend = ndim
        if pos=='front':  # 前乘在某矩阵上
            row1 = Mat.slice(j, j+1, nstart, nend)  # 切片下一个行向量
            row2 = Mat.slice(j+1, j+2, nstart, nend)
            Mat.paste_on(row1.mul(angle[1]) + row2.mul(angle[0]), j, nstart)  # 重新组合后“贴”回原矩阵
            Mat.paste_on(row1.mul(-angle[0]) + row2.mul(angle[1]), j+1, nstart)
            return Mat
        if pos=='back':  # 前乘在某矩阵上
            col1 = Mat.slice(nstart, nend, j, j+1)  # 切片下一个列向量
            col2 = Mat.slice(nstart, nend, j+1, j+2)
            Mat.paste_on(col1.mul(angle[1]) + col2.mul(angle[0]), nstart, j)  # 重新组合后“贴”回原矩阵
            Mat.paste_on(col1.mul(-angle[0]) + col2.mul(angle[1]), nstart, j+1)
            return Mat

    m = H.nrow  # 矩阵阶数
    # 注意：本例中矩阵 index 都是从 0 开始的，与数学习惯不同
    # 以下为 Householder-Hessenberg 约化
    vkts = []  # 储存每个反射矩阵相应的法向量vkt
    for k in range(1, m-1): # k=1...m-2
        xt = matrix([[H.elem[i][k-1] for i in range(k, m)]])  # 截取当前步的列矢量 A_{k+1:m, k}，不过这里用行矢量存储了
        vkt = matrix([[sign(xt.elem[0][0]) * xt.mod() if i==0 else 0 for i in range(m-k)]]) + xt  # 反射线法矢量
        vkt = vkt.mul(1 / vkt.mod())  # 归一化

        # 第一步操作：左乘反射矩阵（用讲义中算法）
        u = [sum([vkt.elem[0][i-k] * H.elem[i][j] for i in range(k, m)]) for j in range(k-1, m)]
        for i in range(k, m):
            for j in range(k-1, m):
                H.elem[i][j] -= 2 * vkt.elem[0][i-k] * u[j-k+1]

        # 第二步操作：右乘反射矩阵转置
        u = [sum([H.elem[i][j] * vkt.elem[0][j-k] for j in range(k, m)]) for i in range(0, m)]
        for i in range(0, m):
            for j in range(k, m):
                H.elem[i][j] -= 2 * u[i] * vkt.elem[0][j-k]
        vkts.append(vkt)

    # 此时 H 已约化为上 Hessenberg 矩阵

    # 以下进行 Givens 旋转变换
    k = m  # 找本征值的指标，从 m 逐渐减小
    angles_repo = [] # 储存正交矩阵的每个旋转角度
    it = 0
    while k > 1:
        if abs(H.elem[k-1][k-2]) < 1e-6:  # 找到一个本征值的判停标准
            k = k-1

        # 原点位移
        s = H.elem[k-1][k-1]
        for j in range(k):  # j=0...k-1
            H.elem[j][j] -= s

        angles = []
        for j in range(k-1):  # j=0...k-2
            norm = sqrt(H.elem[j+1][j]**2 + H.elem[j][j]**2)
            angles.append([H.elem[j+1][j]/norm, H.elem[j][j]/norm]) # 当前的[sin, cos]
            mat_update_from_given(H, k, j, angles[-1], pos='front') # 前乘正交阵G使得A_{j+1,j}=0
        # 此时 H 已旋转至上三角矩阵，角度储存在angles中  

        angles_repo.append(angles) 
        for j in range(k-1):
            mat_update_from_given(H, k, j, angles[j], pos='back') # 再后乘矩阵G^T

        # 原点移回来
        for j in range(k):  # j=0...k-1
            H.elem[j][j] += s
        it += 1

        if it > 500:
            return 'fail!' # 循环太多，就不算了
    eigens = [H.elem[i][i] for i in range(m)]
    
    # 此时所有本征值已经存入eigens，下面对eigens进行冒泡排序
    order = [i for i in range(m)]
    def swap(a, b):
        return b, a
    for i in range(m-1, 0, -1):
        for j in range(i):
            if eigens[j] > eigens[j+1]:
                eigens[j], eigens[j+1] = swap(eigens[j], eigens[j+1])
                order[j], order[j+1] = swap(order[j], order[j+1]) #  也要将原始位置按相同方式重排，为了后面本征矢按相同规则换序
    
    if need_eigenstates == False:  # 若无需输出本征矢
        return eigens  # 只返回本征值即可啦
    
    else:  # 若需要本征矢，下面要对各个正交矩阵相乘
        Q = matrix.identity(m)  # 从单位矩阵开始
        
        # Householder-Hessenberg 约化过程的各个正交阵乘在Q上
        for k in range(1, m-1): # k=1...m-2
            u = [sum([Q.elem[i][j] * vkts[k-1].elem[0][j-k] for j in range(k, m)]) for i in range(0, m)]
            for i in range(0, m):
                for j in range(k, m):
                    Q.elem[i][j] -= 2 * u[i] * vkts[k-1].elem[0][j-k]
        
        # Givens过程的各个正交阵乘在Q上
        for angles in angles_repo:
            for j in range(len(angles)):
                mat_update_from_given(Q, m, j, angles[j], pos='back', bandshape=False)
        # 此时Q的各个列向量即为本征矢
        
        # 对Q中的列向量（本征矢）按照本征值排序时的索引调换顺序，得到Q_sort
        Q_sort = matrix.identity(m)
        for i in range(m):
            Q_sort.paste_on(Q.slice(0, m, order[i], order[i]+1), 0, i)
        
        # 有了本征矢可以判断宇称，存入is_even
        is_even = []
        for i in range(m):
            if Q_sort.elem[0][i]/Q_sort.elem[m-2][i] > 0:
                is_even.append(True)
            else:
                is_even.append(False)
        dic = {'eigens': eigens, 'Q': Q_sort, 'is_even': is_even}
        return dic  # 返回本征值，本征矢，是否偶宇称