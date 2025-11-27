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

len = 10000;
step = 0.02;
z = [0;v0];
E0 = 0.5.*v0^2-GM./r0;
E = E0.*ones(1,len);
t=0;
y_ = [r0;0];
len_ = 100*len;
for i=1:len_
    [y_, z] = runge_kutta3(t, y_, z, step);
end
y = y_.*ones(2,len);
tic
for i=1:len-1
    [y(:,i+1), z] = runge_kutta3(t, y(:,i), z, step);
    R = norm(y(:,i+1));
    v = norm(z);
    E(i+1) = 0.5.*v.^2-GM./R;
end
toc
scatter(y(1,:),y(2,:),5,'filled','MarkerEdgeColor','none')
figure
plot(1:len,log10(abs((E-E0)./E0)),'LineWidth',1.1)
figure
plot(1:len,E,'LineWidth',1.1)










function a = odefun2(t, y)
    a = -y./norm(y)^3;
end




function [y_, z_] = runge_kutta3(t, y, z, step)

    k1 = odefun2(t, y);
    k2 = odefun2(t + step/5, y+step/5*z + step^2/50 * k1);
    k3 = odefun2(t + 2*step/3, y+step/3*z*2 + step^2/27 * (-k1+7*k2));
    k4 = odefun2(t + step, y+step*z + step^2*(3/10*k1-2/35*k2+9/35*k3));
    y_ = y + step*z + step^2*(14/336*k1 + 100/336*k2 + 54/336*k3);
    z_ = z + step*(14/336*k1 + 125/336*k2 + 162/336*k3 + 35/336*k4);
end

function [r, v, a] = verlet_half(t, r_old, v_old, a_old, step)
    v_half = v_old + 0.5 * a_old * step;
    r = r_old + v_half * step;
    a = odefun2(t, r);
    v = v_half + 0.5 * a * step;
end
