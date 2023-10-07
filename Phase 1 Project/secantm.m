% Define the function you want to minimize
f = @(x) x.^2 - 10.*exp(0.1.*x);

% Initial guesses for the two points
x0 = -4;
x1 = -2;

% Set the tolerance for convergence
tolerance = 1e-6;

% Maximum number of iterations
max_iterations = 100;

% Initialize iteration counter
iterations = 0;

while abs(x1 - x0) > tolerance && iterations < max_iterations
    % Calculate the secant slope
    slope = (f(x1) - f(x0)) / (x1 - x0);
    
    % Calculate the next point using the secant method
    x_new = x1 - f(x1) / slope;
    
    % Update points for the next iteration
    x0 = x1;
    x1 = x_new;
    
    % Update the iteration counter
    iterations = iterations + 1;
end

% The minimum is at x1 (or x_new)
min_x = x1;
min_value = f(min_x);

fprintf('Minimum x: %f\n', min_x);
fprintf('Minimum value: %f\n', min_value);
fprintf('Number of iterations: %d\n', iterations);
