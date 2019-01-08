import numpy as np
from numpy.linalg import *

# 5X5矩阵，rank = 5
A = np.array([[1.0, 10, 3, 6, 7],
              [-100.0, 5, 4, 1, 9], 
              [2.0, 0, 3, 1, 1],
              [4.0, 4, 9, 2, 0],
              [1.0, 5, 2, 3, 4]])

# 5X5矩阵，rank = 4
B = np.array([[1.0, 10, 3, 6, 7],
              [1.0, 10, 3, 6, 7], 
              [2.0, 0, 3, 1, 1],
              [4.0, 4, 9, 2, 0],
              [1.0, 5, 2, 3, 4]])


b = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
#####################################################################

print(matrix_rank(A))        # 秩

print(np.trace(A))           # 迹

print(det(A))                # 行列式

print(inv(A))                # 逆矩阵

print(solve(A, b))           # 线性方程组

val, vec = eig(A)            # 特征值和特征向量
print(val)
print(vec)

Q, R = qr(A)                   # QR分解
print(Q)
print(R)
