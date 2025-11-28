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

function [A, b, c] = gauss_legendre_butcher_table_vpa(s, digit)
    % 使用高精度VPA计算s级高斯-勒让德拉格朗日龙格-库塔法的Butcher表
    % 输入: s - 级数, digits - 计算精度（有效数字位数）
    % 输出: A - s×s矩阵, b - 权重向量(1×s), c - 节点向量(s×1)
    
    % 设置计算精度
    if nargin < 2
        digit = 32; % 默认32位精度
    end
    digits(digit);
    
    % 使用符号计算和高精度算术
    syms x real
    
    % 步骤1: 计算节点c_i (勒让德多项式的根)
    % 计算s次移位勒让德多项式(定义在[0,1]区间)
    P = legendre_poly(s);
    c = sort(vpasolve(P, x, [0, 1])); % 在[0,1]区间求解根
    
    % 步骤2: 计算权重b_i (使用高斯求积公式)
    b = sym(zeros(1, s));
    for i = 1:s
        L = sym(1);
        for j = 1:s
            if j ~= i
                L = L * (x - c(j)) / (c(i) - c(j));
            end
        end
        b(i) = vpaintegral(L, x, 0, 1);
    end
    
    % 步骤3: 计算系数矩阵A (使用简化条件C(s))
    A = sym(zeros(s, s));
    for i = 1:s
        % 为每个i建立方程组
        eqns = sym(zeros(s, 1));
        vars = sym(zeros(1, s));
        
        % 创建符号变量
        for j = 1:s
            vars(j) = sym(['a', num2str(i), num2str(j)]);
        end
        
        for k = 1:s
            lhs = sym(0);
            for j = 1:s
                lhs = lhs + vars(j) * c(j)^(k-1);
            end
            eqns(k) = lhs == c(i)^k / k;
        end
        
        % 求解方程组
        sol = solve(eqns, vars);
        
        % 提取解
        for j = 1:s
            A(i, j) = sol.(char(vars(j)));
        end
    end
    
    % 验证阶条件 (可选)
    verify_order_conditions_vpa(A, b, c, s, digit);

    
    
    % 显示结果
    fprintf('高斯-勒让德 Butcher 表 (%d级, %d位精度):\n', s, digit);
    disp('c = ');
    disp(vpa(c, digit));
    disp('A = ');
    disp(vpa(A, digit));
    disp('b = ');
    disp(vpa(b, digit));
    
end

function P = legendre_poly(n)
    % 生成n次移位勒让德多项式(定义在[0,1]区间)
    syms x real
    
    % 标准勒让德多项式(定义在[-1,1])
    P_standard = legendreP(n, 2*x-1);
    
    % 归一化系数
    P = simplify(P_standard);
end

function verify_order_conditions_vpa(A, b, c, s, digits)
    % 使用VPA验证阶条件
    fprintf('\n验证阶条件 (应接近0, %d位精度):\n', digits);
    
    % 转换为高精度
    A_vpa = vpa(A, digits);
    b_vpa = vpa(b, digits);
    c_vpa = vpa(c, digits);
    
    % 验证正交条件 ∑b_i = 1
    err_b1 = abs(sum(b_vpa) - 1);
    fprintf('∑b_i - 1 = %s\n', vpa(err_b1, digits));
    
    % 验证正交条件 ∑b_i * c_i^{k-1} = 1/k
    for k = 1:2*s
        lhs = vpa(0, digits);
        for i = 1:s
            lhs = lhs + b_vpa(i) * c_vpa(i)^(k-1);
        end
        err = abs(lhs - 1/k);
        fprintf('∑b_i * c_i^%d - 1/%d = %s\n', k-1, k, vpa(err, digits));
    end
    
    % 验证简化条件C(s) ∑a_ij * c_j^{k-1} = c_i^k / k
    for i = 1:s
        for k = 1:s
            lhs = vpa(0, digits);
            for j = 1:s
                lhs = lhs + A_vpa(i, j) * c_vpa(j)^(k-1);
            end
            err = abs(lhs - c_vpa(i)^k / k);
            fprintf('A(%d,:) * c^%d - c_%d^%d/%d = %s\n', i, k-1, i, k, k, vpa(err, digits));
        end
    end
end