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


% 另一个算例

clc,clear

% 定义微分方程 y'' = f(t, y)
f = @(t, y) -y./norm(y)^2;  
% 初始条件
tspan = [0, 1000];
y0 = [1; 0];
vy0 = [0; 1.2];

% 容差
tol = 1e-16;

% 求解
tic
[t, y, vy] = rkn76t_adaptive(f, tspan, y0, vy0, tol);

toc


eng = y;
for i=1:length(y)
    eng(1,i) = 0.5*norm(vy(:,i))^2+log(norm(y(:,i)));
    eng(2,i) = y(1,i) * vy(2,i) - y(2,i) * vy(1,i);
end

% 绘制结果
figure;
scatter(y(1,:), y(2,:), 5, 'filled','MarkerEdgeColor','none');
xlabel('时间');
ylabel('位置');


eng01 = eng(1,1);
eng02 = eng(2,1);
figure;
plot(t,abs((eng(1,:)-eng01)/eng01),'LineWidth',1.1,'Color','r')
hold on
plot(t,abs((eng(2,:)-eng02)/eng02),'LineWidth',1.1,'Color','b')
figure;
plot(t,log10(abs((eng(1,:)-eng01)/eng01)),'LineWidth',1.1,'Color','r')
hold on
plot(t,log10(abs((eng(2,:)-eng02)/eng02)),'LineWidth',1.1,'Color','b')

%%
opts1 = odeset(RelTol=tol,MaxStep=1e-1);
tic
[t,y]=ode89(@odefun, tspan, [y0;vy0], opts1);
toc
eng = zeros(2,length(y));
y=y';
for i=1:length(y)
    eng(1,i) = 0.5*norm(y(3:4,i))^2+log(norm(y(1:2,i)));
    eng(2,i) = y(1,i) * y(4,i) - y(2,i) * y(3,i);
end

% 绘制结果
figure;
scatter(y(1,:), y(2,:), 5, 'filled','MarkerEdgeColor','none');
xlabel('时间');
ylabel('位置');


eng01 = eng(1,1);
eng02 = eng(2,1);
figure;
plot(t,abs((eng(1,:)-eng01)/eng01),'LineWidth',1.1,'Color','r')
hold on
plot(t,abs((eng(2,:)-eng02)/eng02),'LineWidth',1.1,'Color','b')
figure;
plot(t,log10(abs((eng(1,:)-eng01)/eng01)),'LineWidth',1.1,'Color','r')
hold on
plot(t,log10(abs((eng(2,:)-eng02)/eng02)),'LineWidth',1.1,'Color','b')


function dydt = odefun(t, y)
    r = [y(1);y(2)];
    a = -r./norm(r)^2;
    dydt = [y(3);y(4);a];
end





