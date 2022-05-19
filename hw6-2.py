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
x_i = 0
y_i = 100
x_j = 300
y_j = 100
x_k = 0
y_k = 300
t = 1  # 平面厚度

# 温度条件 ℃
tem_ij = 200  # ij流体温度
tem_jk = 100  # jk流体温度

# 物性参数
tem_ref = 80  # 参考温度
k = 8.0  # W / (mm * K)
h = 0.2  # W / (mm^2 * K)
E = 1.5e5  # 杨氏模量 Pa
miu = 0.25  # 泊松比
expan = 1e-5  # 热膨胀系数 K^-1

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


def area_size(x_i, y_i, x_j, y_j, x_k, y_k):
    """
    三角形单元的面积
    """
    matrix_A = np.array([[1, x_i, y_i], [1, x_j, y_j], [1, x_k, y_k]])
    A = 0.5 * det(matrix_A)
    return A


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


def rigid_matrix_conduct(matrix_B, matrix_D, A):
    """
    导热刚度矩阵
    """
    tmp = np.dot(np.transpose(matrix_B), matrix_D)
    return np.dot(tmp, matrix_B) * t / 4 / A


def rigid_matrix_convect(matrix_N, h, A):
    """
    对流刚度矩阵
    包括表面和两条边
    """
    L_ij = 300
    L_jk = math.sqrt(300 * 300 + 200 * 200)
    # tmp = h * np.transpose(matrix_N) * matrix_N
    # bound_convect = h * np.transpose(matrix_N) * matrix_N
    # surface = h * A / 12 * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    bound_ij = h * t * L_ij / 6 * np.array([[2, 1, 0], [1, 2, 0], [0, 0, 0]])
    bound_jk = h * t * L_jk / 6 * np.array([[0, 0, 0], [0, 2, 1], [0, 1, 2]])
    return bound_ij + bound_jk
    # return bound_ij + bound_jk + surface


"""
求传热刚度矩阵
"""

matrix_T = np.array([T_i, T_j, T_k])
A = area_size(x_i, y_i, x_j, y_j, x_k, y_k)
N_i, N_j, N_k = N_calculation(
    x, y, x_i, y_i, a_i, b_i, c_i, x_j, y_j, a_j, b_j, c_j, x_k, y_k, a_k, b_k, c_k
)
# print("Ni =", N_i)
# print("Nj =", N_j)
# print("Nk =", N_k)

matrix_N26 = matrix_N26_calculation(N_i, N_j, N_k)
# print("N26\n", matrix_N26, "\n")

matrix_N13 = matrix_N13_calculation(N_i, N_j, N_k)
# print("N13\n", matrix_N13, "\n")

func_T = matrix_N13 * np.transpose(matrix_T)
func_T = func_T[0]  # 前一步的func_T是array类型，不能直接求偏导
# print("T\n", func_T, "\n")

matrix_g = matrix_g_calculation(func_T, x, y)
# print("g\n", matrix_g, "\n")

matrix_D = matrix_D_calculation(k)
# print("D\n", matrix_D, "\n")

matrix_B = matrix_B_calculation(b_i, c_i, b_j, c_j, b_k, c_k)
# print("B\n", matrix_B, "\n")

# 导热刚度阵
rigid_matrix_cond = rigid_matrix_conduct(matrix_B, matrix_D, A)
print("rigid matrix cond\n", rigid_matrix_cond, "\n")
print("det\n", det(rigid_matrix_cond), "\n")

# 对流换热刚度阵
rigid_matrix_conv = rigid_matrix_convect(matrix_N13, h, A)
print("rigid matrix conv\n", rigid_matrix_conv, "\n")
print("det\n", det(rigid_matrix_conv), "\n")

# 总刚度阵
rigid_matrix = rigid_matrix_cond + rigid_matrix_conv
print("rigid matrix\n", rigid_matrix, "\n")
print("det\n", det(rigid_matrix), "\n")

"""
求取载荷向量
"""

L_ij = 300
L_jk = math.sqrt(300 * 300 + 200 * 200)
# F_surface = h * tem_ref * A / 3 * np.array([[1], [1], [1]])
F_ij = h * tem_ij * L_ij * t / 2 * np.array([[1], [1], [0]])
F_jk = h * tem_jk * L_jk * t / 2 * np.array([[0], [1], [1]])
F = F_ij + F_jk
print("F\n", F, "\n")

"""
求解单元方程
"""

matrix_Tem = np.linalg.solve(rigid_matrix, F)
print("三节点温度\n", matrix_Tem, "\n")

"""
求解热应力
"""
print("\n求解热应力\n")
# 单元平均温升
delta_T = matrix_Tem[0][0] + matrix_Tem[1][0] + matrix_Tem[2][0] - 3 * tem_ref
# print(delta_T)

# 热应变
epsilon_0 = expan * delta_T * np.array([[1], [1], [0]])
# print(epsilon_0)


matrix_D = (
    E / (1 - miu * miu) * np.array([[1, miu, 0], [miu, 1, 0], [0, 0, (1 - miu) / 2]])
)

matrix_B_i = 1 / (2 * A) * np.array([[b_i, 0], [0, c_i], [b_i, c_i]])
matrix_B_j = 1 / (2 * A) * np.array([[b_j, 0], [0, c_j], [b_j, c_j]])
matrix_B_k = 1 / (2 * A) * np.array([[b_k, 0], [0, c_k], [b_k, c_k]])
# print(matrix_B_i)
# print(matrix_B_j)
# print(matrix_B_k)
matrix_B = np.hstack((matrix_B_i, matrix_B_j))
matrix_B = np.hstack((matrix_B, matrix_B_k))
# print(matrix_B)

# 应力矩阵[S]
matrix_S = np.dot(matrix_D, matrix_B)
print("应力矩阵\n", matrix_S, "\n")

# 单元刚度阵[Ke]
matrix_Ke = np.dot(np.dot(np.transpose(matrix_B), matrix_D), matrix_B) * t * A

print("Ke\n", matrix_Ke, "\n")
# print(det(matrix_Ke))

# 热载荷向量
F_T = (
    E
    * t
    * expan
    * delta_T
    / 2
    / (1 - miu)
    * np.array([[b_i], [c_i], [b_j], [c_j], [b_k], [c_k]])
)
print("热载荷向量\n", F_T, "\n")
# area = symbols("area")
# nf = 2 * A * np.dot(np.dot(np.transpose(matrix_B), matrix_D), epsilon_0)
# f = 1 / 2 / area
# sf = integrate(f, area)
# F_T = sf.subs({area: A}) * nf
# print("热载荷向量\n", F_T, "\n")

# 施加边界条件，k固定
# for i in range(6):
#     for j in range(6):
#         if i == 0 or j == 0 or i == 1 or j == 1:
#             matrix_Ke[i][j] = 0

# matrix_Ke[0][0] = 1
# matrix_Ke[1][1] = 1

matrix_Ke_bound = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
for i in range(4):
    for j in range(4):
        matrix_Ke_bound[i][j] = matrix_Ke[i][j]

print("Ke with boundary\n", matrix_Ke_bound, "\n")
# print(det(matrix_Ke_bound))

# matrix_u = np.linalg.solve(matrix_Ke_bound, F_T)
# print("三节点应变\n", matrix_u, "\n")
