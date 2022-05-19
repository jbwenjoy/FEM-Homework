from symtable import Symbol
import numpy as np
from numpy.linalg import *
from sympy import symbols, diff, integrate
import math

"""
SymPy库求偏导和积分示例

x, y = symbols("x, y")
f = x**2 + y**2 + x * y + 2
dx = diff(f, x)
f_partial_x = dx.subs({x: 1, y: 2})

x = symbols("x")
f = x**2
I = integrate(f, (x, 1, 3))
"""
np.set_printoptions(linewidth=400)

# 形状参数 mm
x_i = 20
y_i = 20
x_j = 60
y_j = 20
x_k = 40
y_k = 60
t = 1  # 平面厚度

# 温度条件 ℃
tem_inf = 20

# 物性参数
k = 7.5  # W / (mm * K)
h = 0.15  # W / (mm^2 * K)

a_i = x_j * y_k - x_k * y_j
b_i = y_j - y_k
c_i = x_k - x_j

a_j = x_k * y_i - x_i * y_k
b_j = y_k - y_i
c_j = x_j - x_k

a_k = x_i * y_j - x_j * y_i
b_k = y_i - y_j
c_k = x_j - x_i

x, y = symbols("x, y")
T_i, T_j, T_k = symbols("T_i, T_j, T_k")

# mat = np.array([[a_i, b_i, c_i], [a_j, b_j, c_j], [a_k, b_k, c_k]])
# print("\nabc")
# print(mat)
# mat = 1 / 800 / 2 * np.dot(mat, np.array([[1], [x], [y]]))
# print("\nmat")
# print(mat)


def area_size(x_i, y_i, x_j, y_j, x_k, y_k):
    """
    三角形单元的面积
    """
    matrix_A = np.array([[1, x_i, y_i], [1, x_j, y_j], [1, x_k, y_k]])
    return 0.5 * det(matrix_A)


def N_calculation(
    x, y, x_i, y_i, a_i, b_i, c_i, x_j, y_j, a_j, b_j, c_j, x_k, y_k, a_k, b_k, c_k
):
    """
    求Ni Nj Nk
    """
    A = area_size(x_i, y_i, x_j, y_j, x_k, y_k)
    N_i = (a_i + b_i * x + c_i * y) / (2 * A)
    N_j = (a_j + b_j * x + c_j * y) / (2 * A)
    N_k = (a_k + b_k * x + c_k * y) / (2 * A)
    return N_i, N_j, N_k


def matrix_N26_calculation(N_i, N_j, N_k):
    """
    [N26] = [Ni 0 Nj 0 Nk 0
             0 Ni 0 Nj 0 Nk]
    """
    return np.array([[N_i, 0, N_j, 0, N_k, 0], [0, N_i, 0, N_j, 0, N_k]])


def matrix_N13_calculation(N_i, N_j, N_k):
    """
    [N13] = [Ni Nj Nk]
    """
    return np.array([N_i, N_j, N_k])


def matrix_B_calculation(b_i, c_i, b_j, c_j, b_k, c_k):
    """
    [B] = [bi bj bk
           ci cj ck]
    """
    return np.array([[b_i, b_j, b_k], [c_i, c_j, c_k]])


def matrix_D_calculation(k):
    """
    [D] = [k 0
           0 k]
    """
    return np.array([[k, 0], [0, k]])


def matrix_g_calculation(T, x, y):
    """
    [g] = [diff(T, x)
           diff(T, y)]
    """
    return np.array([[diff(T, x)], [diff(T, y)]])


def stiff_matrix_conduct(matrix_B, matrix_D, A):
    """
    导热刚度矩阵
    """
    tmp = np.dot(np.transpose(matrix_B), matrix_D)
    return np.dot(tmp, matrix_B) * t / 4 / A


def stiff_matrix_convect(matrix_N, h, A):
    """
    对流刚度矩阵
    包括表面和jk边
    """
    L_jk = math.sqrt(20 * 20 + 40 * 40)
    # tmp = h * np.transpose(matrix_N) * matrix_N
    # bound_convect = h * np.transpose(matrix_N) * matrix_N
    surface = h * A / 12 * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    bound_jk = h * t * L_jk / 6 * np.array([[0, 0, 0], [0, 2, 1], [0, 1, 2]])
    # return bound_ij + bound_jk
    return bound_jk + surface


"""
求传热刚度矩阵
"""

matrix_T = np.array([T_i, T_j, T_k])
A = area_size(x_i, y_i, x_j, y_j, x_k, y_k)

N_i, N_j, N_k = N_calculation(
    x, y, x_i, y_i, a_i, b_i, c_i, x_j, y_j, a_j, b_j, c_j, x_k, y_k, a_k, b_k, c_k
)
# N_i = mat[0]
# N_j = mat[1]
# N_k = mat[2]
print("\nNi =", N_i)
print("Nj =", N_j)
print("Nk =", N_k)

# N矩阵
matrix_N26 = matrix_N26_calculation(N_i, N_j, N_k)
# print("N26\n", matrix_N26, "\n")
matrix_N13 = matrix_N13_calculation(N_i, N_j, N_k)
# print("N13\n", matrix_N13, "\n")

# func_T是array类型，不能直接求偏导
func_T = matrix_N13 * np.transpose(matrix_T)
print("\nT")
print(func_T)

matrix_g = matrix_g_calculation(func_T[0], x, y)
# 温度梯度g矩阵
for i in range(1, func_T.size):
    tmp = matrix_g_calculation(func_T[i], x, y)
    matrix_g = np.hstack((matrix_g, tmp))
print("\ng")
print(matrix_g)

# D矩阵
matrix_D = matrix_D_calculation(k)
# print("D\n", matrix_D, "\n")

# B矩阵
matrix_B = matrix_B_calculation(b_i, c_i, b_j, c_j, b_k, c_k)
# print("B\n", matrix_B, "\n")

# 导热刚度阵
stiff_matrix_cond = stiff_matrix_conduct(matrix_B, matrix_D, A)
print("\nstiff matrix cond")
print(stiff_matrix_cond)
print("det =", det(stiff_matrix_cond), "\n")

# 对流换热刚度阵
stiff_matrix_conv = stiff_matrix_convect(matrix_N13, h, A)
print("\nstiff matrix conv")
print(stiff_matrix_conv)
print("det =", det(stiff_matrix_conv), "\n")

# 总传热刚度阵
stiff_matrix = stiff_matrix_cond + stiff_matrix_conv
print("\nstiff matrix")
print(stiff_matrix)
print("det =", det(stiff_matrix), "\n")

"""
求取载荷向量
"""

L_jk = math.sqrt(20 * 20 + 40 * 40)
F_surface = h * tem_inf * A / 3 * np.array([[1], [1], [1]])
F_jk = h * tem_inf * L_jk * t / 2 * np.array([[0], [1], [1]])
F = F_surface + F_jk
print("\nF")
print(F)

"""
求解单元方程
"""

matrix_Tem = np.linalg.solve(stiff_matrix, F)
print("\n三节点温度")
print(matrix_Tem, "\n")
