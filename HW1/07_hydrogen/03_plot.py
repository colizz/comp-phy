# -*- coding: utf-8 -*-

from numpy import sin, cos, pi, exp, sqrt, arctan # 导入函数

xmax, Nx = 8, 200 # [-xmax, xmax]的范围，划为 2Nx+1 个点
X = [i*xmax/Nx for i in range(-Nx, Nx+1)]
Y = [i*xmax/Nx for i in range(-Nx, Nx+1)]
Z1p = [[0 for i in range(-Nx, Nx+1)] for i in range(-Nx, Nx+1)] # 对应 Ψ_311+Ψ_31-1
Z1m = [[0 for i in range(-Nx, Nx+1)] for i in range(-Nx, Nx+1)] # 对应 Ψ_311-Ψ_31-1
Z0 = [[0 for i in range(-Nx, Nx+1)] for i in range(-Nx, Nx+1)] # 对应 Ψ_310
for j in range(len(Z0)):
    for i in range(len(Z0[j])):
        phi = arctan(Y[j]/X[i]) if X[i]!=0 else pi/2
        r2 = X[i]**2 + Y[j]**2
        fac = 1/(32*pi) * exp(-sqrt(r2)) * r2
        Z1p[i][j] = fac * cos(phi)**2
        Z1m[i][j] = fac * sin(phi)**2
        Z0[i][j] = fac
Z1p = Z1p[::-1]
Z1m = Z1m[::-1]
Z0 = Z0[::-1]

#########################################################################
## 使用 matlibplot 作图
from matplotlib import pyplot as plt
im = plt.imshow(Z1p, extent=(-xmax, xmax, -xmax, xmax))
plt.colorbar(im)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('Q7_3_Psi211p.pdf')
plt.show()
im = plt.imshow(Z1m, extent=(-xmax, xmax, -xmax, xmax))
plt.colorbar(im)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('Q7_3_Psi211m.pdf')
plt.show()
im = plt.imshow(Z0, extent=(-xmax, xmax, -xmax, xmax))
plt.colorbar(im)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('Q7_3_Psi210.pdf')
plt.show()