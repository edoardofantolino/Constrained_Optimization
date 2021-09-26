function [gradfx] = fd_grad(f, x, h, type)

h = h*norm(x);

switch type
    case 'Jfw'
        gradfx = (f(x+h*(ones))-f(x))/h;
    case 'Jbw'
        gradfx = (f(x)-f(x-h*(ones)))/h;
    case 'Jc'
        gradfx = (f(x+h*(ones))-f(x-h*(ones)))/(2*h);
end

end