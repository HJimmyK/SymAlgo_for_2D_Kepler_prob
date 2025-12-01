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

% 8阶RKN辛分裂方法（A19型）计算彗星椭圆运动
% 参考论文：Runge–Kutta–Nyström symplectic splitting methods of order 8
clc, clear;  close all;

%% 1. 基础参数设置（彗星-太阳二体系统）
mu = 398600;           % 引力参数(km^3/s^2)，太阳-彗星系统取mu=GM（G=6.67e-11, M=1.989e30）
e = 0.9;               % 彗星轨道偏心率（椭圆，e<1）
a = 2.5e8;             % 轨道半长轴(km)
c0 = a*e;
T = 2*pi*sqrt(a^3/mu); % 轨道周期(s)，由开普勒第三定律计算
t_final = 100*T;       % 积分总时间
h = T/320;             % 积分步长（每个周期80步）
N = round(t_final/h);  % 总迭代步数

% 初始条件
r0 = a*(1 - e);       % 近心距（初始位置距太阳距离）
v0 = sqrt(mu*(1 + e)/(a*(1 - e))); % 初始速度（近心点速度）
y0 = [r0, 0];         % 初始位置向量 (x, y) [km]
v0 = [0, v0];         % 初始速度向量 (vx, vy) [km/s]

% 存储结果的数组
y = zeros(N+1, 2);    % 位置
v = zeros(N+1, 2);    % 速度
E = zeros(N+1, 1);    % 能量
y(1, :) = y0;
v(1, :) = v0;
E(1) = (norm(v0)^2)/2 - mu/norm(y0); % 初始能量（哈密顿量）

%% 2. 8阶RKN方法（A19型）系数定义（参考论文表4，保留15位有效数字）
A = [
    0.0505805, ...  % a1
     0.149999, ... % a2
     -0.0551795510771615573511026950361, ... % a3
    0.423755898835337951482264998051, ...% a4
    -0.213495353584659048059672194633, ... % a5
    -0.0680769774574032619111630736274,  ...% a6
    0.227917056974013435948887201671, ... % a7
    -0.235373619381058906524740047732,  ...% a8
    0.387413869179878047816794031058   ... % a9
];
A(10) = 0.5-sum(A(1:9));

% 回文结构扩展a系数（A19型：a1,a2,...,a9,a9,...,a2,a1）
A = [A, fliplr(A(1:end))];
B = [
     0.129478606560536730662493794395,  ...% b1
    0.222257260092671143423043559581, ... % b2
    -0.0577514893325147204757023246320, ... % b3
    -0.0578312262103924910221345032763, ... % b4
    0.103087297437175356747933252265, ... % b5
    -0.140819612554090768205554103887, ... % b6
    0.0234462603492826276699713718626,  ...% b7
    0.134854517356684096617882205068, ... % b8
    0.0287973821073779306345172160211  ... % b9
];
B(10) = 1-2*sum(B(1:9));
% 回文结构扩展b系数（A19型：b1,b2,...,b9,b8,...,b2,b1）
B = [B, fliplr(B(1:end-1))];

s = length(B); % 阶段数（19，与论文A19一致）

%% 3. 8阶RKN辛分裂方法数值积分
tic
for n = 1:N
    % 当前状态
    y_curr = y(n, :);
    v_curr = v(n, :);
    
    % 执行A19型流组合：phi(a1*h, a) -> phi(b1*h, b) ->...-> phi(a1*h, a)
    for i = 1:s
        % 1. 应用自由运动流 phi_{h*a(i)}^{[a]}: y = y + h*a(i)*v, v不变
        y_curr = y_curr + h * A(i) * v_curr;
        
        % 2. 应用势能作用流 phi_{h*b(i)}^{[b]}: v = v + h*b(i)*g(y), y不变
        r = norm(y_curr);          % 彗星到太阳的距离
        g = -mu / (r^3) * y_curr;  % 引力加速度（g(y) = -mu*y/r^3
        v_curr = v_curr + h * B(i) * g;
    end
    % 最后一步自由运动流（A19型首尾为a流，补充a(1)项）
    y_curr = y_curr + h * A(end) * v_curr;
    
    % 更新状态与能量
    y(n+1, :) = y_curr;
    v(n+1, :) = v_curr;
    E(n+1) = (norm(v_curr)^2)/2 - mu / norm(y_curr);
end
toc

%% 4. 结果分析与可视化
% 4.1 彗星轨道图
figure(Color='w');
t = 0:h:t_final;
t = t./(365*24*60*60);
subplot(2,1,1)
semilogy(t, abs(((y(:,1)+c0).^2 + y(:,2).^2./(1-e^2))./(a^2) - 1), 'b-', 'LineWidth', .2);
hold on;
xlim([0,max(t)]);
xlabel('时间 (年) '); ylabel('椭圆偏差的对数');
title('椭圆偏差');
grid on; 

% 4.2 能量误差图（验证辛特性，能量应近似守恒）
subplot(2,1,2)
E_error = abs((E - E(1)) / E(1)); % 相对能量误差
semilogy(t, E_error, 'r-', 'LineWidth', .3);
%semilogy(t, E_error, 'r-', 'LineWidth', 1.2);
xlim([0,max(t)]);
xlabel('时间 (年)'); ylabel('总能量相对误差的对数');
title('总能量误差');
grid on;

figure(Color='w');
scatter(y(:,1),y(:,2),10,'blue','filled','MarkerEdgeColor','none')
xlabel('x'); ylabel('y');
title('基于8阶RKN方法的彗星椭圆运动轨道');
grid on; 

fprintf('积分总时间：%.2f s（%.2f 个轨道周期）\n', t_final, t_final/T);
fprintf('初始能量：%.6e km^2/s^2\n', E(1));
fprintf('最大相对能量误差：%.2e\n', max(E_error));
fprintf('8阶RKN方法阶段数：%d\n', s);







