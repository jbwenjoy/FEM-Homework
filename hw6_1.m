%# 形状参数 mm
x_i = 20; y_i = 20; x_j = 60;
y_j = 20; x_k = 40; y_k = 60;
t = 1;  % 平面厚度

%# 温度条件 ℃
tem_inf = 20;

%# 物性参数
k = 7.5;  %# W / (mm * K)
h = 0.15;  %# W / (mm^2 * K)

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
% disp("abc");
% disp(mat);

mat = 1 / 800 / 2 * mat * [1; x; y];
% disp("mat");
% disp(mat);
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
disp("stiff matrix cond"), disp(stiff_matrix_cond)

%# 对流换热刚度阵
L_jk = sqrt(20 * 20 + 40 * 40);
surface = h * A / 12 * [2, 1, 1; 1, 2, 1; 1, 1, 2];
bound_jk = h * t * L_jk / 6 * [0, 0, 0; 0, 2, 1; 0, 1, 2];
stiff_matrix_conv = surface + bound_jk;
disp("stiff matrix conv"), disp(stiff_matrix_conv);

%# 传热刚度矩阵
stiff_matrix = stiff_matrix_cond + stiff_matrix_conv;
disp("stiff matrix"), disp(stiff_matrix);

%# 热载荷
L_jk = sqrt(20 * 20 + 40 * 40);
F_surface = h * tem_inf * A / 3 * [1; 1; 1];
F_jk = h * tem_inf * L_jk * t / 2 * [0; 1; 1];
F = F_surface + F_jk;
disp("F"), disp(F);

matrix_Tem = stiff_matrix \ F;
disp("三节点温度"), disp(matrix_Tem);
