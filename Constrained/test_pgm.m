close all
clear
clc

n = 1e4;
x_min = -0.682327803828019*ones(n,1); % actual minumum
x = -2*rand(n,1)+2;
alpha0 = 1;
kmax = 50;
c1 = 1e-4;
rho = 0.8;
btmax = 2;
tollgrad = 1e-8;
gamma = 0.1;
tolx = 1e-6;
left = 1;
right = 2;
Pi_X = @(x) projection(x, left, right);

% function handle for the function and the gradient
f = @(x) sum((1/4)*x.^4 + (1/2)*x.^2 + x);
gradf = @(x) [x.^3 + x + 1];
fg = @(x) [(1/4)*x.^4 + (1/2)*x.^2 + x];


num_iter = 100;
times = zeros(num_iter,1);
steps = zeros(num_iter,1);

for kh=2:2:14
    h = 10^(-kh);
    for i=1:num_iter
        tic;
        [xk, fk, gradfk_norm, deltaxk_norm, k, xseq] = ...
        projected_gradient_method(x, f, gradf, alpha0, ...
        kmax, tollgrad, gamma, tolx, Pi_X, h, 'Jc');
        times(i) = toc;
        steps(i) = k;
        xseq = [x xseq];
    end 
    kh
    result = xk(end-1)
    avg_steps = mean(steps)
    avg_time = mean(times)
end