% Solve Quadratic Program
%
%   minimize    0.5*x'*A*x + b'*x
%      x
%   subject to  -a <= x <= a
%
% with proximal (projected) gradient method


% ====================================================================== 
% Parameters

n = 20;    % Dimension of the variable

a = 10;    % Box constraint on the variable

MAX_ITER = 10000; % Maximum number of iterations
eps_stop = 1e-4;  % Criterion for stopping

error_std = 0; % Standard deviation of error 
% ====================================================================== 

% ====================================================================== 
% Generate data

[Q, R] = qr(randn(n, n));
A = Q*diag([50*rand(floor(0.5*n), 1); zeros(n - floor(0.5*n), 1)])*Q';
b = rand(n, 1);

ones_n  = ones(n, 1);

L = max(eigs(A));
% ====================================================================== 


% ====================================================================== 
% Solve with CVX

cvx_begin quiet
    variable x_opt(n, 1);
    minimize(0.5*x_opt'*A*x_opt + b'*x_opt);
    subject to
        -a*ones_n <= x_opt <= a*ones_n;
cvx_end
% ====================================================================== 


% ====================================================================== 
% Proximal gradient descent

x_k = zeros(n, 1);

error_hist = zeros(MAX_ITER, 1);

for k = 1 : MAX_ITER

    % Save previous point
    x_k_1 = x_k;

    % Compute gradient
    grad_g = A*x_k + b;

    new_point = x_k - (1/L)*grad_g;

    % Project onto [-a, a]
    x_k = max(min(new_point, a*ones_n), -a*ones_n) + error_std*randn(n,1);

    error_hist(k) = norm(x_k - x_k_1)/norm(x_k);

    if norm(x_k - x_k_1) <= eps_stop
        break;
    end
end
% ====================================================================== 

% Check error
fprintf('||x_k - x_opt||/||x_opt|| = %f\n', norm(x_k - x_opt)/norm(x_opt))
fprintf('Number of iterations      = %d (out of %d)\n', k, MAX_ITER)

% Visualize error
figure(1);clf;
semilogy(1:k, error_hist(1:k));
xlabel('Iteration')
title('||x_k - x_{opt}||/||x_{opt}||')
