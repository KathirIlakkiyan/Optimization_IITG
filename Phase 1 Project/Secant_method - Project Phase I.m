clc
clear

a_ini = input('Enter lower bound: ');
b_ini = input('Enter upper bound: ');
opt = input('Enter min/max: ','s');

[a, b] = Bounding_phase(a_ini, b_ini, opt);
Secant_method(a, b, opt);

function obj_fun = f(x)
    obj_fun = 4*x*sin(x);
end

function Secant_method(x0, x1, opt)
    tolerance = input('Epsilon: ');
    max_iterations = input('Maximum iterations: ');
    ite = 1;
    delta = (x1 - x0)/max_iterations;
    f_dash_0 = Diff(x0, delta);
    f_dash_1 = Diff(x1, delta);
    % Output file
    out = fopen('Secant_method.out', 'w');
    fprintf(out,'x0 = %4.3f\tx1 = %4.3f\n\n',x0,x1);
    fprintf(out, 'Iter\tx0\tx1\t\tz\t\tf_dash(x0)\tf_dash(x1)\n');
    while abs(x1 - x0) > tolerance && ite < max_iterations+1
        slope = (Diff(x1,delta)-Diff(x0,delta))/(x1-x0);
        z = x1 - Diff(x1,delta)/slope;
        fprintf(out, '%d\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t\t%4.3f\n',ite,x0,x1,z,f_dash_0,f_dash_1);
        % Update for the next iteration
        x0 = x1;
        x1 = z;
        ite = ite + 1;
    
        f_dash_0 = Diff(x0, delta);
        f_dash_1 = Diff(x1, delta);
    end
    fprintf(out, '%d\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t\t%4.3f\n',ite,x0,x1,z,f_dash_0,f_dash_1);
    
    if opt == "min"
        fprintf('\nMinimum x: %f\n', x1);
        fprintf('Minimum value: %f\n', f(x1));
        fprintf('Number of iterations: %d\n', ite);
    elseif opt == "max"
        fprintf('\nMaximum x: %f\n', x1);
        fprintf('Maximum value: %f\n', f(x1));
        fprintf('Number of iterations: %d\n', ite);
    end
end



function f_dash = Diff(x, delta)
    f_dash = (f(x-delta)-f(x+delta))/(2*delta);
end

function [a, b] = Bounding_phase(a_ini, b_ini, opt)
    a = a_ini;
    b = b_ini;
    % Parameters for the Bounding Phase method
    delta = input('Delta: ');     % Initial step size for moving towards the minimum
    x0 = a + (b-a).*rand(1,1);
    x1 = x0;
    k = 0;
    iter = 0;
    
    if opt == "min"
        % Step size direction is towards the minimum
        if f(x0 - abs(delta)) > f(x0)
            delta = abs(delta);
        else
            delta = -abs(delta);
        end
        % Output file
        out = fopen('Bounding_Phase_iterations.out', 'w');
        fprintf(out, 'Iter\tdelta\tx0\tx1\t\tf(x0)\tf(x1)\ta\t\tb\n');
    
        while f(x0 + (2^k) * delta) < f(x0)
            x0 = x1;
            x1 = x0 + (2^k) * delta;
            % Update the interval bounds
            if delta>0
                a = x0;
                b = x1;
            else
                a = x1;
                b = x0;
            end
            k = k + 1;
            iter = iter + 1;
            fprintf(out, '%d\t%1.4f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n',iter,delta, x0,x1,f(x0),f(x0 + (2^k) * delta),a,b);
        end
    
        fprintf('Final range: (%f, %f)\n', a, b);
        fprintf('No. of iterations: %f\n\n', iter);
    
    elseif opt == "max"
        % Step size direction is towards the maximum
        if f(x0 - abs(delta)) > f(x0)
            delta = -abs(delta);
        else
            delta = abs(delta);
        end
    
        out = fopen('Bounding_Phase_iterations.out', 'w'); % Output file
        fprintf(out, 'Iter\tdelta\tx0\tx1\t\tf(x0)\tf(x1)\ta\t\tb\n');
    
        while f(x0 + (2^k) * delta) > f(x0)
            x0 = x1;
            x1 = x0 + (2^k) * delta;
            % Update the interval bounds
            if delta>0
                a = x0;
                b = x1;
            else
                a = x1;
                b = x0;
            end
            k = k + 1;
            iter = iter + 1;
            fprintf(out, '%d\t%1.4f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n',iter,delta, x0,x1,f(x0),f(x1),a,b);
        end
    
        fprintf('Final range: (%f, %f)\n', a, b);
        fprintf('No. of iterations: %f\n\n', iter);
    end
end