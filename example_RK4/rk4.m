% MIT License
% 
% Copyright (c) 2025 Jecricho Knox
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% 四阶龙格库塔法求解二维开普勒问题（含偏差与相对误差分析）
clear; clc; close all;

%% 物理参数设置
mu = 398600;           % 地球引力参数，单位: km^3/s^2
t0 = 0;                % 初始时间，单位: s
tf = 86400;            % 仿真总时间，单位: s（1天）
dt = 100;              % 时间步长，单位: s
t = t0:dt:tf;          % 时间序列
n = length(t);         % 时间步数

%% 初始条件设置（近地圆轨道）
r0 = 6378 + 500;       % 初始圆轨道半径，单位: km（地球半径6378km + 500km轨道高度）
v0 = sqrt(mu / r0);    % 圆轨道理想速度，单位: km/s
z0 = [r0, 0, 0, v0];   % 初始状态向量: [x0, dx0/dt, y0, dy0/dt]

z = zeros(n, 4);       % 每一行存储一个时刻的状态向量 [x, dx/dt, y, dy/dt]
z(1, :) = z0;          % 赋初始值

%% 四阶龙格库塔法（RK4）迭代
for i = 1:n-1
    z_current = z(i, :);
    t_current = t(i);
    
    k1 = kepler_deriv(t_current, z_current, mu);
    k2 = kepler_deriv(t_current + dt/2, z_current + dt/2 * k1, mu);
    k3 = kepler_deriv(t_current + dt/2, z_current + dt/2 * k2, mu);
    k4 = kepler_deriv(t_current + dt, z_current + dt * k3, mu);
    
    % RK4迭代公式
    z(i+1, :) = z_current + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% 计算误差
x = z(:, 1);            % x坐标序列
dx = z(:, 2);           % x方向速度序列
y = z(:, 3);            % y坐标序列
dy = z(:, 4);           % y方向速度序列

% 5.1 轨道半径绝对偏差：x^2 + y^2 - r0^2
r_square = x.^2 + y.^2; % 数值解的轨道半径平方
r_dev = r_square - r0^2;% 绝对偏差（理想圆轨道下应为0）

% 5.2 机械能与相对误差（单位质量）
r = sqrt(r_square);     % 数值解的轨道半径
v = sqrt(dx.^2 + dy.^2);% 数值解的速度大小
E = v.^2/2 - mu./r;     % 单位质量机械能
E0 = E(1);              % 初始机械能（参考值）
E_rel_err = (E - E0)./abs(E0) * 100; % 机械能相对误差（%）

% 5.3 角动量与相对误差（单位质量）
L = x.*dy - y.*dx;      % 单位质量角动量
L0 = L(1);              % 初始角动量（参考值）
L_rel_err = (L - L0)./abs(L0) * 100; % 角动量相对误差（%）

%% 结果可视化
figure('Color','w','Position',[100,100,1200,800]); % 调整窗口大小

% 子图1：轨道轨迹
subplot(2,2,1);
plot(x, y, 'b-', 'LineWidth', 1.2);
hold on;
theta = 0:0.01:2*pi;
earth_x = 6378 * cos(theta);
earth_y = 6378 * sin(theta);
fill(earth_x, earth_y, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'g');
axis equal;  % 等比例坐标轴，保证轨道形状正确
xlabel('X坐标 (km)');
ylabel('Y坐标 (km)');
title('二维开普勒问题圆轨道仿真（RK4法）');
legend('卫星轨道','地球','Location','best');
grid on;

% 子图2：轨道半径绝对偏差
subplot(2,2,2);
plot(t/3600, r_dev, 'r-', 'LineWidth', 1);
xlabel('时间 (h)');
ylabel('轨道半径绝对偏差 (km^2)');
title('圆轨道半径绝对偏差（x^2 + y^2 - r_0^2）');
grid on;

% 子图3：机械能相对误差
subplot(2,2,3);
plot(t/3600, E_rel_err, 'k-', 'LineWidth', 1);
xlabel('时间 (h)');
ylabel('机械能相对误差 (%)');
title('单位质量机械能相对误差');
grid on;

% 子图4：角动量相对误差
subplot(2,2,4);
plot(t/3600, L_rel_err, 'm-', 'LineWidth', 1);
xlabel('时间 (h)');
ylabel('角动量相对误差 (%)');
title('单位质量角动量相对误差');
grid on;

sgtitle('RK4求解开普勒圆轨道：轨迹与数值精度验证','FontSize',14);
tightfig;

%% 定义开普勒问题的状态方程
function dz = kepler_deriv(t, z, mu)
    % 输入：
    % t: 当前时间
    % z: 当前状态向量 [x, dx/dt, y, dy/dt]
    % mu: 引力参数
    % 输出：
    % dz: 状态向量的导数 [dx/dt, d^2x/dt^2, dy/dt, d^2y/dt^2]
    
    x = z(1);
    dx = z(2);
    y = z(3);
    dy = z(4);
    
    r = sqrt(x^2 + y^2);  
    r3 = r^3;             
    
    % 状态方程的导数
    dz = [
        dx;
        -mu * x / r3;
        dy;
        -mu * y / r3;
    ];
    dz = dz';  % 转为行向量
end

function tightfig()
    ax = findobj(gcf,'Type','Axes');
    for i = 1:length(ax)
        pos = ax(i).Position;
        ax(i).Position = pos - [0.02, 0.02, 0, 0];
    end
end