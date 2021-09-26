function [Hessfx] = j_hess(f, x, h, type)

h = h*norm(x);

switch type
    case 'Jfw'
        Hessfx = spdiags((f(x+h*(ones))-f(x))/h, 0, length(x), length(x));
    case 'Jbw'
        Hessfx = spdiags((f(x)-f(x-h*(ones)))/h, 0, length(x), length(x));
    case 'Jc'
        Hessfx = spdiags((f(x+h*(ones))-f(x-h*(ones)))/(2*h), 0, length(x), length(x));
    case 'c'
        Hessfx = spdiags((f(x+h*(ones))-(2*f(x))+f(x-h*(ones)))/(h^2), 0, length(x), length(x));
end

end
