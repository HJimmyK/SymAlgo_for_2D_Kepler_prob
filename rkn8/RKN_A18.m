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
% RKN 8th-order splitting method for Kepler problem
% Based on: Blanes, Casas & Escorihuela-Tomas (2021)
% "Runge-Kutta-Nystrom symplectic splitting methods of order 8"

% Kepler problem parameters
mu = 1;  % Gravitational parameter
e = 0.9; % Eccentricity
T = 2*pi; % Period
tfinal = 10*T; % Final time
c0 = e;

% Initial conditions (elliptic orbit)
q0 = [(1-e); 0.0];
p0 = [0; sqrt((1+e)/(1-e))];
y0 = [q0; p0];

% Method coefficients (using A18 method from the paper)
[a, b] = getA18Coefficients();

% Step size
h = T/200;

% Integrate using RKN splitting method
t = 0:h:tfinal;
y = zeros(4, length(t));
y(:, 1) = y0;
tic
for i = 2:length(t)
    y(:, i) = rknStep(@keplerFlow, y(:, i-1), h, a, b, mu);
end
toc

% Calculate energy error
E = zeros(1, length(t));
for i = 1:length(t)
    q = y(1:2, i);
    p = y(3:4, i);
    E(i) = 0.5*dot(p, p) - mu/norm(q);
end
E_error = abs(E - E(1)); % Theoretical energy is -0.5

figure(Color='w');
scatter(y(1,:),y(2,:),10,'blue','filled','MarkerEdgeColor','none')
xlabel('q_1');
ylabel('q_2');
title(sprintf('Kepler Orbit (e=%f)',e));
set(gca,'FontName','Times New Roman')

figure(Color='w');
subplot(2, 1, 1);
semilogy(t, abs((y(1,:)+c0).^2 + y(2,:).^2./(1-e^2) - 1), 'b-', 'LineWidth', .1);
xlim([0 tfinal])
xlabel('t');
ylabel('error');
title(sprintf('Kepler Orbit (e=%f)',e));
set(gca,'FontName','Times New Roman')

subplot(2, 1, 2);
semilogy(t, E_error);
xlim([0 tfinal])
xlabel('Time');
ylabel('Energy Error');
title('Energy Conservation');
set(gca,'FontName','Times New Roman')
% Display final error
fprintf('Final energy error: %e\n', E_error(end));
fprintf('Average energy error: %e\n', mean(E_error));

function [a, b] = getA18Coefficients()
    % Coefficients for A18 method from Table 4 of the paper
    % Note: s=18 for A18 (18 stages, 19 a coefficients)
    
    % a coefficients (19 total, including symmetric ones)
    a_tmp = zeros(1, 10);
    a_tmp(1) = 0.0866003822712445920135805954462;
    a_tmp(2) = -0.0231572735424388070228714693753;
    a_tmp(3) = 0.191410576083774088999564416369;
    a_tmp(4) = 0.378895558692931579545387584925;
    a_tmp(5) = -0.0467359566364556111599485526051;
    a_tmp(6) = -0.156198111997810415438979605642;
    a_tmp(7) = 0.156025836895094823718831871041;
    a_tmp(8) = 0.252844012473796333586850465807;
    a_tmp(9) = -0.640644212172254239866860564270;
    a_tmp(10) = 1 - 2*sum(a_tmp(1:9));
    
    % Build full symmetric a sequence: a1..a10, a9..a1 (total 19)
    a = [a_tmp(1:10), a_tmp(9:-1:1)];
    
    % b coefficients (18 total, including symmetric ones)
    b_tmp = zeros(1, 9);
    b_tmp(1) = -0.08;
    b_tmp(2) = 0.209460550048243262121199483001;
    b_tmp(3) = 0.274887805875735483503233064415;
    b_tmp(4) = -0.224214208870409561366168655624;
    b_tmp(5) = 0.347657740563761656321390026010;
    b_tmp(6) = -0.168783183866211679175007668385;
    b_tmp(7) = 0.144209344805460873709120777707;
    b_tmp(8) = 0.0116851121360265483381405054244;
    b_tmp(9) = 0.5 - sum(b_tmp(1:8));
    
    % Build full symmetric b sequence: b1..b10, b9..b1 (total 18)
    b = [b_tmp(1:9), b_tmp(9:-1:1)];
end

function y_next = rknStep(flowFunc, y, h, a, b, mu)
    % One step of RKN splitting method (A-type)
    % y = [q; p] state vector
    % Implementation based on equation (2.4):
    % φ_{h a_{s+1}}^{[a]} ∘ φ_{h b_s}^{[b]} ∘ ... ∘ φ_{h b_1}^{[b]} ∘ φ_{h a_1}^{[a]}
    
    q = y(1:2);
    p = y(3:4);
    
    % Number of stages (should be 18 for A18)
    s = length(b);  % b has 18 elements for A18
    
    % First a-step (a1)
    g = keplerAcceleration(q, mu);
    p = p + a(1) * h * g;
    
    % Main loop for stages 1 to s (s=18)
    for i = 1:s
        % b-step (bi)
        q = q + b(i) * h * p;
        
        % a-step (a_{i+1}) - except for the last iteration
        g = keplerAcceleration(q, mu);
        p = p + a(i+1) * h * g;
        
    end
    
    
    y_next = [q; p];
end

function acc = keplerAcceleration(q, mu)
    % Acceleration for Kepler problem: g(q) = -mu*q/|q|^3
    r = norm(q);
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



