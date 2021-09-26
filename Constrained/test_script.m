close all
clear
clc

n = 1e5;  % value of i
x_min = -0.682327803828019*ones(n,1); % actual minumum
x = rand(n,1);
alpha = 1;
kmax = 70;
tollgrad = 1e-8;

% function handle for the function, the gradient, and
% the hessian
f = @(x) sum((1/4)*x.^4 + (1/2)*x.^2 - x);
gradf = @(x) [x.^3 + x + 1];
Hf = @(x) sparse(1:n,1:n, 3*x.^2+1);
fg = @(x) [(1/4)*x.^4 + (1/2)*x.^2 - x];

h = 10^-4;
Hessian = Hf(x);
HessianJacobian = fd_hess(gradf, x, h, 'Jc');
HessianFromf = fd_hess(fg, x, h, 'c');
M = [diag(Hessian) diag(HessianJacobian) diag(HessianFromf)];
M = full(M);

typeg = '';
typeh = '';
% for the iteration results
for kh=2:2:14
    h = 10^(-kh);
    [xk, fk, gradfk_norm, k, xseq] = newton_method(x, fg, gradf, Hf, alpha, kmax, tollgrad, h, typeg, typeh); 
    iter = kh
    k
    xk;
end


% for time results
for kh=2:2:14
    h = 10^(-kh);
    num_iter = 5;
    times = zeros(num_iter, 1);
    for j = 1:num_iter
        tic;
        [xk, fk, gradfk_norm, k, xseq] = newton_method(x, fg, gradf, Hf, alpha, kmax, tollgrad, h, typeg, typeh); 
        times(j,1) = toc;
    end
    figure(1)
    avg_time = mean(times)
    plot(kh,mean(times),'o')
    grid on
    hold on
end
