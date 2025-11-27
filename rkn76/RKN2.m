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

clc,clear


GM = 1;
r0 = 1;
v0 = 1.3*sqrt(GM/r0);

T = 100000;
step = 0.02;
e_step = 200*step;
z = [0;v0];
E0 = 0.5.*v0^2-GM./r0;
E = E0.*ones(1,ceil(T./e_step));
t=0;
y = [r0;0];
i=0;
tic
while t<T
    [y, z] = runge_kutta_n(t, y, z, step);
    t=t+step;

    if t >= i*e_step
        i = i+1;
        R = norm(y);
        v = norm(z);
        E(i+1) = 0.5.*v.^2-GM./R;
    end
end
toc

figure
error = abs((E-E0)./E0);
plot(1:length(E),log10(error),'LineWidth',1.1)
figure
plot(1:length(E),E,'LineWidth',1.1)

%%
clc,clear
format long g
len = 10;
x = [0.4896193895465455...
0.2698116326738158...
0.9897401420622155...
0.1836756962661240...
0.8616567249240227...
0.0326325721064562...
0.3319580195866526...
0.7487469957829571...
0.6443664581816437...
0.1692379626516646];
for i=1:len
    fprintf('%1.16f+...\n',x(i))
end
tic
sum1 = sum(x);
toc
tic
sum2 = 0.4896193895465455+...
0.2698116326738158+...
0.9897401420622155+...
0.1836756962661240+...
0.8616567249240227+...
0.0326325721064562+...
0.3319580195866526+...
0.7487469957829571+...
0.6443664581816437+...
0.1692379626516646;
toc
tic
sum3 = kahan_sum(x);
toc
sum1-sum2
sum1-sum3






function a = odefun2(t, y)
    a = -y./norm(y)^3;
end


function [y_, z_] = runge_kutta_n(t, y, z, step)

    k1 = odefun2(t, y);
    k2 = odefun2(t + step/5, y+step/5*z + step^2/50 * k1);
    k3 = odefun2(t + 2*step/3, y+step/3*z*2 + step^2/27 * (-k1+7*k2));
    k4 = odefun2(t + step, y+step*z + step^2*(3/10*k1-2/35*k2+9/35*k3));
    y_ = y + step*z + step^2*(14/336*k1 + 100/336*k2 + 54/336*k3);
    z_ = z + step*(14/336*k1 + 125/336*k2 + 162/336*k3 + 35/336*k4);
end

function [y_, z_] = runge_kutta_kahan(t, y, z, step)

    k1 = odefun2(t, y);
    k2 = odefun2(t + step/5, kahan_sum([y, step/5*z , step^2/50 * k1]));
    k3 = odefun2(t + 2*step/3,kahan_sum([ y, step/3*z*2 , step^2/27 * (-k1+7*k2)]));
    k4 = odefun2(t + step, kahan_sum([y, step*z , step^2*( 3/10*k1 - 2/35*k2 + 9/35*k3)]) );
    y_ = y + step*z + step^2*( kahan_sum([14/336*k1 , 100/336*k2 , 54/336*k3]) );
    z_ = z + step*( kahan_sum([14/336*k1 , 125/336*k2 , 162/336*k3 , 35/336*k4]));
end

function [r, v, a] = verlet_half(t, r_old, v_old, a_old, step)
    v_half = v_old + 0.5 * a_old * step;
    r = r_old + v_half * step;
    a = odefun2(t, r);
    v = v_half + 0.5 * a * step;
end

function s = kahan_sum(vec)
    sum_total = 0.0;     % 累加和
    c = 0.0;             % 误差补偿变量
    for i = 1:length(vec)
        y = vec(:,i) - c;  % 减去之前累积的误差
        t = sum_total + y;
        c = (t - sum_total) - y;  % 计算新的误差
        sum_total = t;
    end
    s = sum_total;
end

