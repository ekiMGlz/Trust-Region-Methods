function [grad_f] = gradient(f, x)
%GRADIENT Approximate the gradient of f at x
%   Input:
%       f: Anonymous function
%       x: Function arguments
%   Output:
%       grad_f: Gradient of f at x
%

    if ~iscolumn(x)
        x = x';
    end

    n = length(x);
    grad_f = zeros(n, 1);
    E = spdiags(eps^(1/3) * (abs(x) + 1), 0, n, n);

    for i = 1:n
        grad_f(i) = ( f(x + E(:, i)) - f(x - E(:, i)) ) * (0.5 / E(i, i));
    end
    
end

