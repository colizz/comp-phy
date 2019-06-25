# -*- coding: utf-8 -*-
import sys
sys.path.append('..') # 为导入上级目录的mymatrixlib库
from mymatrixlib.matrix import matrix
from func01 import mat_initialize, find_eigenvalues

H = mat_initialize(M=2, q=1)  # 初始化矩阵
H0 = H.copy()  # 把初始 H 记下来
result = find_eigenvalues(H, need_eigenstates=True)  # 计算本征值、本征矢
## 测试结果
('H = '); H0.print()
print('Q = '); result['Q'].print()
print('Q^T H Q = '); (result['Q'].T() * H0 * result['Q']).print()  # 验证结果为对角阵