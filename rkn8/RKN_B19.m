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

clc, clear
% RKN 8th-order splitting method (B19) for Kepler problem
% Based on: Blanes, Casas & Escorihuela-Tomas (2021)
% "Runge-Kutta-Nystrom symplectic splitting methods of order 8"

% Kepler problem parameters
mu = 1;         % Gravitational parameter
e = 0.8;        % Eccentricity
T = 2*pi;       % Period
tfinal = 100*T; % Final time
c0 = e;

% Initial conditions (elliptic orbit)
q0 = [(1-e); 0.0];
p0 = [0; sqrt((1+e)/(1-e))];
y0 = [q0; p0];

% Method coefficients (using B19 method from Table 5 of the paper)
[a, b] = getB19Coefficients();

% Display method info
fprintf('B19 Method:\n');
fprintf('Number of a coefficients: %d\n', length(a));
fprintf('Number of b coefficients: %d\n', length(b));
fprintf('Sum of a coefficients: %.15f\n', sum(a));
fprintf('Sum of b coefficients: %.15f\n', sum(b));

% Test different step sizes
h_values = [T/20, T/40, T/80, T/160, T/320];
error_data = zeros(length(h_values), 3);

for h_idx = 1:length(h_values)
    h = h_values(h_idx);
    
    % Integrate using RKN splitting method
    t = 0:h:tfinal;
    y = zeros(4, length(t));
    y(:, 1) = y0;
    
    tic
    for i = 2:length(t)
        y(:, i) = rknBStep(y(:, i-1), h, a, b, mu);
    end
    compute_time = toc;
    
    % Calculate energy error
    E = zeros(1, length(t));
    for i = 1:length(t)
        q = y(1:2, i);
        p = y(3:4, i);
        E(i) = 0.5*dot(p, p) - mu/norm(q);
    end
    E_error = abs(E - E(1));
    
    % Calculate position error using analytical solution for Kepler
    % (simplified check using orbit equation)
    position_error = zeros(1, length(t));
    for i = 1:length(t)
        q = y(1:2, i);
        % Check if point lies on ellipse: (x+c0)^2 + y^2/(1-e^2) = 1
        position_error(i) = abs((q(1)+c0)^2 + q(2)^2/(1-e^2) - 1);
    end
    
    error_data(h_idx, 1) = h;
    error_data(h_idx, 2) = max(E_error);
    error_data(h_idx, 3) = max(position_error);
    
    fprintf('\n时间步长 h = %.6f (T/%d):\n', h, round(T/h));
    fprintf('  最大能量误差: %.2e\n', max(E_error));
    fprintf('  最大位置误差: %.2e\n', max(position_error));
    fprintf('  计算时间: %.3f s\n', compute_time);
    
    % Plot trajectory for medium step size
    if h_idx == 3
        figure(Color='w');
        scatter(y(1,:), y(2,:), 10, 'blue', 'filled', 'MarkerEdgeColor', 'none');
        xlabel('q_1');
        ylabel('q_2');
        title(sprintf('B19 Method - Kepler Orbit (e=%.1f, h=T/%.0f)', e, T/h));
        axis equal;
        grid on;
    end
end

% Plot convergence
figure(Color='w');
subplot(2, 1, 1);
loglog(error_data(:,1), error_data(:,2), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
loglog(error_data(:,1), error_data(:,3), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Step size h','FontName','Times New Roman');
ylabel('Maximum Error','FontName','Times New Roman');
title('Convergence of B19 Method','FontName','Times New Roman');
legend('Energy Error', 'Position Error', 'Location', 'best');
set(gca,'FontName','Times New Roman')
grid on;

% Add reference slopes for order 8
ref_h = error_data(:,1);
ref_e = 100 * ref_h.^8;
loglog(ref_h, ref_e, 'k--', 'LineWidth', 1);
legend('Energy Error', 'Position Error', 'O(h^8) reference', 'Location', 'best');
set(gca,'FontName','Times New Roman')
set(gca,'FontName','Times New Roman')

subplot(2, 1, 2);
semilogy(0:h:tfinal, E_error, 'b-', 'LineWidth', 1);
xlabel('Time','FontName','Times New Roman');
ylabel('Energy Error','FontName','Times New Roman');
title(sprintf('Energy Error (h=T/%.0f)', T/h_values(3)),'FontName','Times New Roman');
set(gca,'FontName','Times New Roman')
grid on;


% Function to get B19 coefficients from Table 5
function [a, b] = getB19Coefficients()
    % Coefficients for B19 method from Table 5 of the paper
    % Note: s=19 for B19, so we have s+1=20 b coefficients and s=19 a coefficients
    
    % a coefficients (first 9 from table, then calculate a10, then symmetrize)
    a_tmp = zeros(1, 10);
    a_tmp(1) = 0.337548675291317241942440116575;
    a_tmp(2) = -0.223647977575409990331768222380;
    a_tmp(3) = 0.168949714872223740906385138015;
    a_tmp(4) = 0.171179938816205886154783136334;
    a_tmp(5) = -0.349765168067292877221144631312;
    a_tmp(6) = 0.523808861006312397712070357524;
    a_tmp(7) = -0.194208871063049124066394765282;
    a_tmp(8) = -0.323496751337931087309823477561;
    a_tmp(9) = 0.322817287614899749216601693799;
    
    % Calculate a10 using the formula from Table 5: a10 = 1 - 2*sum_{i=1}^{9} a_i
    a_tmp(10) = 1 - 2*sum(a_tmp(1:9));
    
    % For B-type methods with s=19, the full symmetric a sequence is:
    % a1, a2, ..., a10, a9, a8, ..., a1 (total 19)
    a = [a_tmp(1:10), a_tmp(9:-1:1)];
    
    % b coefficients (first 9 from table, then calculate b10, then symmetrize)
    b_tmp = zeros(1, 10);
    b_tmp(1) = 0.036132460472136313416730168194;
    b_tmp(2) = 0.012697863961074113381675193011;
    b_tmp(3) = 0.201318391240629276109068041836;
    b_tmp(4) = 0.135683350134504233201330671671;
    b_tmp(5) = -0.0579071833999963041504740663015;
    b_tmp(6) = -0.0772509501792649549463874931821;
    b_tmp(7) = -0.00264758266409925952822161203471;
    b_tmp(8) = -0.0329844384945603065320797537355;
    b_tmp(9) = 0.0476781560950366927530646289755;
    
    % Calculate b10 using the formula from Table 5: b10 = 1/2 - sum_{i=1}^{9} b_i
    b_tmp(10) = 0.5 - sum(b_tmp(1:9));
    
    % For B-type methods with s=19, the full symmetric b sequence is:
    % b1, b2, ..., b10, b10, b9, ..., b1 (total 20)
    % Note the symmetry: b_{s+2-i} = b_i, with s=19
    b = [b_tmp(1:10), b_tmp(10), b_tmp(9:-1:1)];
end

function y_next = rknBStep(y, h, a, b, mu)
    % One step of B-type RKN splitting method
    % Based on equation (2.5): 
    % B_s = φ_{h b_{s+1}}^{[b]} ∘ φ_{h a_s}^{[a]} ∘ φ_{h b_s}^{[b]} ∘ ... ∘ φ_{h a_1}^{[a]} ∘ φ_{h b_1}^{[b]}
    
    q = y(1:2);
    p = y(3:4);
    
    s = length(a);  % Should be 19 for B19
    s_b = length(b); % Should be 20 for B19
    
    if s_b ~= s+1
        error('For B-type methods, b should have one more coefficient than a');
    end
    
    % First b-step (b1)
    g = keplerAcceleration(q, mu);
    p = p + b(1) * h * g;
    
    % Main loop for i = 1 to s (s=19)
    for i = 1:s
        % a-step (a_i)
        q = q + a(i) * h * p;
        
        % b-step (b_{i+1})
        g = keplerAcceleration(q, mu);
        p = p + b(i+1) * h * g;
    end
    
    y_next = [q; p];
end

function acc = keplerAcceleration(q, mu)
    % Acceleration for Kepler problem: g(q) = -mu*q/|q|^3
    r = norm(q);
    if r < 1e-12
        error('Position vector too close to zero');
    end
    acc = -mu * q / (r^3);
end

function dy = keplerFlow(t, y, mu)
    % Kepler flow for ODE45 comparison
    q = y(1:2);
    p = y(3:4);
    
    dqdt = p;
    dpdt = keplerAcceleration(q, mu);
    
    dy = [dqdt; dpdt];
end