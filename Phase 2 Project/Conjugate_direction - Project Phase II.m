clc
clear

N = 3; % No. of variables
% Define the function to minimize
%f = @(x) sum((1:N).*x(1:N).^2);
%f = @(x) sum(100*(x(2:N) - x(1:N-1).^2).^2 + (x(1:N-1) - 1).^2);
%f = @(x) (x(1) - 1).^2 + sum((1:N).*(2.*x(2:end).^2 - x(1:end-1)).^2);
%f = @(x) sum((x - 1).^2) - sum(x(2:end).*x(1:end-1));
f = @(x) sum(x.^2) + (sum(0.5.*(1:N).*x)).^2 + (sum(0.5.*(1:N).*x)).^4;

a = -5 * ones(N, 1);
b = 10 * ones(N, 1);
x_current = a + (b - a) .* rand(N, 1);
tolerance = 0.5e-6;
max_iterations = 100;
iterations = 0;

% Define the identity matrix for initial directions
s = eye(N);
        
while iterations < max_iterations
    x_initial = x_current;
    
    for i = 1:length(x_current)
        % Unidirectional search
        [a, b] = Bounding_phase(f, x_current, s(:,i));
        alpha = Secant_method(a, b, f, x_current, s(:,i), tolerance, 100);
        
        % Update the current point
        x_current = x_current + alpha * s(:,i);
    end
    d = x_current - x_initial;
    fprintf('Iteration: %d  f(x) = %f\n', iterations+1, sum(f(x_current)))
    % Check for termination
    if norm(d) < tolerance
        iterations = iterations + 1;
        break;
    end
    for k = length(x_initial):2
    s(:,k) = s(:,k-1);
    end
    s(:,1) = d./norm(d);
    iterations = iterations + 1;
end

% The minimum is at x_current
min_x = x_current;
min_value = sum(f(min_x));

fprintf('Minimum x: %s\n', mat2str(min_x));
fprintf('Minimum value: %f\n', min_value);
fprintf('Number of iterations: %d\n', iterations);

% Secant Method for unidirectional search
function alpha = Secant_method(x0, x1, f, x, d, tol, max_iterations)
    ite = 0;
    delta = (x1 - x0)/max_iterations;
    F = @(alpha) f(x + alpha * d);
    while abs(x1 - x0) > tol && ite < max_iterations
        slope = (Diff(F,x1,delta)-Diff(F,x0,delta))/(x1-x0);
        z = x1 - Diff(F,x1,delta)/slope;
        
        % Update for the next iteration
        x0 = x1;
        x1 = z;
        ite = ite + 1;
    end
    alpha = x1;
    % fprintf('\nMinimum x: %f\n', x1);
    % fprintf('Minimum value: %f\n', f(x1));
    % fprintf('Number of iterations: %d\n', ite);
end

% Unidirectional search bracketting
function [a, b] = Bounding_phase(f, x, d)
    a = -20;
    b = 20;
    F = @(alpha) f(x + alpha * d);
    % Parameters for the Bounding Phase method
    delta = 0.15;
    x0 = a + (b-a).*rand(1,1);
    x1 = x0;
    k = 0;
    iter = 0;
    % Step size direction is towards the minimum
    if F(x0 - abs(delta)) > F(x0)
        delta = abs(delta);
    else
        delta = -abs(delta);
    end

    while F(x0 + (2^k) * delta) < F(x0)
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
     end

    % fprintf('Final range: (%f, %f)\n', a, b);
    % fprintf('No. of iterations: %f\n\n', iter);
end

function f_dash = Diff(fun, x, delta)
    f_dash = (fun(x-delta)-fun(x+delta))/(2*delta);
end