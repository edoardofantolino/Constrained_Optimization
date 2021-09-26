function [xk, fk, gradfk_norm, k, xseq] = ...
    newton_method(x0, f, gradf, Hf, alpha, kmax, tollgrad, h, typeg, typeh)

switch typeg
    case 'Jfw'
        gradf = @(x) fd_grad(f, x, h, 'Jfw');
        disp('Gradient approximatino with forward method')
    case 'Jbw'
        gradf = @(x) fd_grad(f, x, h, 'Jbw');
        disp('Gradient approximatino with backward method')
    case 'Jc'
        gradf = @(x) fd_grad(f, x, h, 'Jc');
        disp('Gradient approximatino with centered method')
    otherwise
end

switch typeh
    case 'Jfw'
        Hf = @(x) fd_hess(gradf, x, h, 'Jfw');
        disp('Hessian approximatino with forward method')
    case 'Jbw'
        Hf = @(x) fd_hess(gradf, x, h, 'Jbw');
        disp('Hessian approximatino with backward method')
    case 'Jc'
        Hf = @(x) fd_hess(gradf, x, h, 'Jc');
        disp('Hessian approximatino with centered method')
    case 'c'
        Hf = @(x) fd_hess(f, x, h, 'c');
        disp('Hessian approximatino with centered method no(J)')
    otherwise
end

% Initializations
xseq = zeros(length(x0), kmax);

xk = x0;
k = 0;
gradfk_norm = norm(gradf(xk));

while k < kmax && gradfk_norm >= tollgrad

    pk = -Hf(xk)\gradf(xk);
    xk = xk + alpha * pk;
    
    gradfk_norm = norm(gradf(xk));
    k = k + 1;
    xseq(:, k) = xk;
end

fk = f(xk);
xseq = xseq(:, 1:k);

end