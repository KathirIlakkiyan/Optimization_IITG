clc
clear
% Minimize function
f = @(x) x.^2 - 10.*exp(0.1.*x);

% Initial interval [a, b]
a = input('Enter lower bound: ');
b = input('Enter upper bound: ');

% Parameters for the Bounding Phase method
delta = input('Delta: ');     % Initial step size for moving towards the maximum
x0 = input('Enter the initial guess: ');
x1 = x0;
k = 0;
iter = 0;
% Step size direction is towards the maximum
if f(x0 - abs(delta)) > f(x0)
    delta = abs(delta);
else
    delta = -abs(delta);
end

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
end

% if delta<0
%     a = x0 + (2^(k-1)) * delta;
%     b = x0;
% else
%     a = x0;
%     b = x0 + (2^(k-1)) * delta;
% end

fprintf('Final range: (%f, %f)\n', a, b);
fprintf('No. of iterations: %f\n', iter);
