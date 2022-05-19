%# 形状参数 mm
x_i = 0; y_i = 100; x_j = 300; y_j = 100; x_k = 0; y_k = 300;
t = 1;  % 平面厚度

%# 温度条件 ℃
tem_ij = 200; tem_jk = 100;

%# 物性参数
tem_ref = 80;
k = 8.0;
h = 0.2;
E = 1.5e5;
miu = 0.25;
expan = 1e-5;

a_i = x_j * y_k - x_k * y_j;
b_i = y_j - y_k;
c_i = x_k - x_j;

a_j = x_k * y_i - x_i * y_k;
b_j = y_k - y_i;
c_j = x_j - x_k;

a_k = x_i * y_j - x_j * y_i;
b_k = y_i - y_j;
c_k = x_j - x_i;

syms x  y T_i T_j T_k;

matrix_T = [T_i, T_j, T_k];
matrix_A = [1, x_i, y_i; 1, x_j, y_j; 1, x_k, y_k];
A = (0.5 * det(matrix_A));

mat = [a_i, b_i, c_i; a_j, b_j, c_j; a_k, b_k, c_k];
disp("abc"), disp(mat);

mat = 1 / 800 / 2 * mat * [1; x; y];
disp("mat"), disp(mat);
N_i = mat(1,1);
N_j = mat(2,1);
N_k = mat(3,1);
disp("Ni ="), disp(N_i);
disp("Nj ="), disp(N_j);
disp("Nk ="), disp(N_k);

% matrix_N26 = [N_i, 0, N_j, 0, N_k, 0; 0, N_i, 0, N_j, 0, N_k];
% disp("N26");
% disp(matrix_N26);
matrix_N13 = [N_i, N_j, N_k];
disp("N13"), disp(matrix_N13);

%# func_T是array类型，不能直接求偏导
func_T = matrix_N13 * matrix_T';
disp("T"), disp(func_T);

matrix_g = [diff(func_T, x); diff(func_T, y)];

%# D矩阵
matrix_D = [k 0; 0 k];

%# B矩阵
matrix_B = [b_i, b_j, b_k; c_i, c_j, c_k];

%# 导热刚度阵
stiff_matrix_cond = matrix_B' * matrix_D * matrix_B * t / 4 / A;
disp("stiff matrix cond"), disp(stiff_matrix_cond);

%# 对流换热刚度阵
L_ij = 300;
L_jk = sqrt(300 * 300 + 200 * 200);
bound_ij = h * t * L_ij / 6 * [2, 1, 0; 1, 2, 0; 0, 0, 0];
bound_jk = h * t * L_jk / 6 * [0, 0, 0; 0, 2, 1; 0, 1, 2];
stiff_matrix_conv = bound_ij + bound_jk;
disp("stiff matrix conv"), disp(stiff_matrix_conv);

%# 传热刚度矩阵
stiff_matrix = stiff_matrix_cond + stiff_matrix_conv;
disp("stiff matrix"), disp(stiff_matrix);

%# 热载荷
L_jk = sqrt(300 * 300 + 200 * 200);
F_ij = h * tem_ij * L_jk * t / 2 * [1; 1; 0];
F_jk = h * tem_jk * L_jk * t / 2 * [0; 1; 1];
F = F_ij + F_jk;
disp("F"), disp(F);

matrix_Tem = stiff_matrix \ F;
disp("三节点温度"), disp(matrix_Tem);

%# 单元平均温升
delta_T = matrix_Tem(1, 1) + matrix_Tem(2, 1) + matrix_Tem(3, 1) - 3 * tem_ref;
% disp(delta_T);

%# 热应变
epsilon_0 = expan * delta_T * [1; 1; 0];
disp(epsilon_0)

matrix_D = E / (1 - miu * miu) * [1, miu, 0; miu, 1, 0; 0, 0, (1 - miu) / 2];

matrix_B_i = 1 / (2 * A) * [b_i, 0; 0, c_i; c_i, b_i];
matrix_B_j = 1 / (2 * A) * [b_j, 0; 0, c_j; c_j, b_j];
matrix_B_k = 1 / (2 * A) * [b_k, 0; 0, c_k; c_k, b_k];

matrix_B = [matrix_B_i matrix_B_j matrix_B_k];
disp("B"), disp(matrix_B);

%# 应力矩阵[S]
matrix_S = matrix_D * matrix_B;
disp("应力矩阵"), disp(matrix_S);

%# 单元刚度阵[Ke]
matrix_Ke = matrix_B' * matrix_D * matrix_B * t * A;
disp("Ke"), disp(matrix_Ke);

%# 热载荷向量
F_T = E * t * expan * delta_T / 2 / (1 - miu) * [b_i; c_i; b_j; c_j; b_k; c_k];
disp("热载荷向量"), disp(F_T);

%# 施加边界条件
matrix_Ke_bound = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
F_T_bound = [0; 0; 0; 0];
for i = 1:4
    F_T_bound(i, 1) = F_T(i, 1);
    for j = 1:4
        matrix_Ke_bound(i, j) = matrix_Ke(i, j);
    end
end
% matrix_Ke(1, 1) = 1e50;
% matrix_Ke(2, 2) = 1e50;
disp("Ke with boundary"), disp(matrix_Ke_bound);
disp("F_T with boundary"), disp(F_T_bound);

matrix_u = matrix_Ke_bound \ F_T_bound;
disp("三节点应变"), disp(matrix_u);
