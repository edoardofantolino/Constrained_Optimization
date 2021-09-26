function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq] = ...
    projected_gradient_method(x0, f, gradf, alpha0, ...
    kmax, tollgrad, gamma, tolx, Pi_X, h, typeg)

switch typeg
    case 'Jfw'
        gradf = @(x) fd_grad(f, x, h, 'Jfw');
%         disp('Gradient approximatino with forward method')
    case 'Jbw'
        gradf = @(x) fd_grad(f, x, h, 'Jbw');
%         disp('Gradient approximatino with backward method')
    case 'Jc'
        gradf = @(x) fd_grad(f, x, h, 'Jc');
%         disp('Gradient approximatino with centered method')
    otherwise
end

% Initializations
xseq = zeros(length(x0), kmax);

xk = x0; % Project the starting point if outside the constraints
fk = f(xk);
k = 0;
gradfk_norm = norm(gradf(xk));
deltaxk_norm = tolx + 1;

while k < kmax && gradfk_norm >= tollgrad && deltaxk_norm >= tolx
    xk = Pi_X(xk);
    
    % Compute the descent direction
    pk = -gradf(xk);
    
    xbark = xk + gamma * pk;
    xhatk = Pi_X(xbark);    
    
    % Reset the value of alpha
    alpha = alpha0;
    
    % Compute the candidate new xk
    pik = xhatk - xk;
    xnew = xk + alpha * pik;
    
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    % Update xk, fk, gradfk_norm, deltaxk_norm
    deltaxk_norm = norm(xnew - xk);
    xk = xnew;
    fk = fnew;
    gradfk_norm = norm(gradf(xk));
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
end

% "Cut" xseq to the correct size
xseq = xseq(:, 1:k);

end