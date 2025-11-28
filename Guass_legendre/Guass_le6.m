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

clear
% 高斯-勒让德龙格库塔法求解开普勒问题
% 绘制轨迹和能量、角动量误差（对数坐标）

% 物理参数
mu = 1.0;      % 引力参数 (GM)
e = 0.6;       % 偏心率
a = 1.0;       % 半长轴
c0 = e*a;
% 初始条件 (椭圆轨道)
q0 = [a*(1-e); 0];                    % 初始位置
p0 = [0; sqrt((mu/a)*(1+e)/(1-e))];   % 初始动量
y0 = [q0; p0];                        % 初始状态向量 [q; p]

% 时间参数
T = 2*pi*sqrt(a^3/mu); % 轨道周期
tspan = [0, 1000*T];      % 模拟1000个周期
h = T/100;             % 步长

% 高斯-勒让德 RK4 系数
[c,A,b] = Guass_legendre_butcher_6();

% 初始化
t = tspan(1):h:tspan(2);
n = length(t);
y = zeros(4, n);
y(:, 1) = y0;

% 存储能量和角动量误差
E_error = zeros(1, n);
L_error = zeros(1, n);

% 精确的能量和角动量值
E_exact = 0.5*dot(p0, p0) - mu/norm(q0);
L_exact = q0(1)*p0(2) - q0(2)*p0(1);


Y1_old = y(:,1);
Y2_old = y(:,1);
Y3_old = y(:,1);
Y4_old = y(:,1);
Y5_old = y(:,1);
Y6_old = y(:,1);
Y1_error = 0;
Y2_error = 0;
Y3_error = 0;
Y4_error = 0;
Y5_error = 0;
Y6_error = 0;
% 主循环
tic
for i = 1:n-1
    % 当前状态
    y_current = y(:, i);
    
    q_current = y_current(1:2);
    p_current = y_current(3:4);

    % Y1 = y_n + h*(a11*f(Y1) + a12*f(Y2))
    % Y2 = y_n + h*(a21*f(Y1) + a22*f(Y2))


        Y1 = Y1_old + Y1_error;
        Y2 = Y2_old + Y2_error;
        Y3 = Y3_old + Y3_error;
        Y4 = Y4_old + Y4_error;
        Y5 = Y5_old + Y5_error;
        Y6 = Y6_old + Y6_error;
        F1 = kepler_rhs(t(i)+c(1)*h, Y1, mu);
        F2 = kepler_rhs(t(i)+c(2)*h, Y2, mu);
        F3 = kepler_rhs(t(i)+c(3)*h, Y3, mu);
        F4 = kepler_rhs(t(i)+c(4)*h, Y4, mu);
        F5 = kepler_rhs(t(i)+c(5)*h, Y5, mu);
        F6 = kepler_rhs(t(i)+c(6)*h, Y6, mu);

        % 类似Guass-Seidel的迭代方式
    for index = 1:12
        Y1 = (y_current + h*( ...
                             A(1,1)*F1 + ...
                             A(1,2)*F2 + ...
                             A(1,3)*F3 + ...
                             A(1,4)*F4 + ...
                             A(1,5)*F5 + ...
                             A(1,6)*F6 ...
                             ));
        F1 = kepler_rhs(t(i)+c(1)*h, Y1, mu);
        Y2 = (y_current + h*( ...
                             A(2,1)*F1 + ...
                             A(2,2)*F2 + ...
                             A(2,3)*F3 + ...
                             A(2,4)*F4 + ...
                             A(2,5)*F5 + ...
                             A(2,6)*F6 ...
                             ));
        F2 = kepler_rhs(t(i)+c(2)*h, Y2, mu);
        Y3 = (y_current + h*( ...
                             A(3,1)*F1 + ...
                             A(3,2)*F2 + ...
                             A(3,3)*F3 + ...
                             A(3,4)*F4 + ...
                             A(3,5)*F5 + ...
                             A(3,6)*F6 ...
                             ));
        F3 = kepler_rhs(t(i)+c(3)*h, Y3, mu);
        Y4 = (y_current + h*( ...
                             A(4,1)*F1 + ...
                             A(4,2)*F2 + ...
                             A(4,3)*F3 + ...
                             A(4,4)*F4 + ...
                             A(4,5)*F5 + ...
                             A(4,6)*F6 ...
                             ));
        F4 = kepler_rhs(t(i)+c(4)*h, Y4, mu);
        Y5 = (y_current + h*( ...
                             A(5,1)*F1 + ...
                             A(5,2)*F2 + ...
                             A(5,3)*F3 + ...
                             A(5,4)*F4 + ...
                             A(5,5)*F5 + ...
                             A(5,6)*F6 ...
                             ));
        F5 = kepler_rhs(t(i)+c(5)*h, Y5, mu);
        Y6 = (y_current + h*( ...
                             A(6,1)*F1 + ...
                             A(6,2)*F2 + ...
                             A(6,3)*F3 + ...
                             A(6,4)*F4 + ...
                             A(6,5)*F5 + ...
                             A(6,6)*F6 ...
                             ));
        F6 = kepler_rhs(t(i)+c(6)*h, Y6, mu);


    end

        Y1_error = Y1-Y1_old;
        Y2_error = Y2-Y2_old;
        Y3_error = Y3-Y3_old;
        Y4_error = Y4-Y4_old;
        Y5_error = Y5-Y5_old;
        Y6_error = Y6-Y6_old;
        Y1_old = Y1;
        Y2_old = Y2;
        Y3_old = Y3;
        Y4_old = Y4;
        Y5_old = Y5;
        Y6_old = Y6;




    % 计算下一步
    y_next = y_current + h*( ...
                            b(1)*F1 + ...
                            b(2)*F2 + ...
                            b(3)*F3 + ...
                            b(4)*F4 + ...
                            b(5)*F5 + ...
                            b(6)*F6   ...
                            );
    
    % 存储结果
    y(:, i+1) = y_next;
    
    % 计算当前能量和角动量误差
    q = y_next(1:2);
    p = y_next(3:4);
    E_current = 0.5*dot(p, p) - mu/norm(q);
    L_current = q(1)*p(2) - q(2)*p(1);
    
    E_error(i+1) = abs(E_current - E_exact);
    L_error(i+1) = abs(L_current - L_exact);
end
toc
% 提取位置坐标
qx = y(1, :);
qy = y(2, :);

% 绘图

figure(Color='w');

% 子图1: 轨迹
scatter(qx,qy,10,'blue','filled','MarkerEdgeColor','none')
xlabel('x');
ylabel('y');
title('Kepler Orbit error (e=0.5)');
axis equal;
grid on;

% 子图1: 轨迹的椭圆偏差

figure(Color='w');
semilogy(t, abs((qx+c0).^2 + qy.^2./(1-e^2) - a^2), 'b-', 'LineWidth', .9);
xlim(tspan)
hold on;
plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('t');
ylabel('error');
title('Kepler Orbit error (e=0.5)');
grid on;

% 子图2: 能量和角动量误差（对数坐标）
figure(Color='w')

semilogy(t, E_error, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Energy Error');
hold on;
semilogy(t, L_error, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Angular Momentum Error');
xlim(tspan)
xlabel('Time');
ylabel('Error (log scale)');
title('Conservation Laws Error');
legend('show', 'Location','southeast');
grid on;

% 显示平均误差
fprintf('Average energy error: %.4e\n', mean(E_error));
fprintf('Average angular momentum error: %.4e\n', mean(L_error));


function dydt = kepler_rhs(t, y, mu)
    % 开普勒问题的右手边函数
    % y = [q; p] = [x; y; vx; vy]
    % dydt = [v; a] = [vx; vy; ax; ay]
    
    q = y(1:2);
    p = y(3:4);
    
    r = norm(q);
    a = -mu * q / r^3;
    
    dydt = [p; a];
end