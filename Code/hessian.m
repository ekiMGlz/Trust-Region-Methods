function [hess_f] = hessian(f,x)
%HESSIAN Approximate the hessian of f at x
%   Input:
%       f: Anonymous function
%       x: Function arguments
%   Output:
%       hess_f: Hessian of f at x
%

    if ~iscolumn(x)
        x = x';
    end

    n = length(x);
    hess_f = zeros(n);
    E = spdiags(eps^(1/4) * (abs(x) + 1), 0, n, n);

    for i = 1:n
        for j = 1:n
            hess_f(i, j) = (   f(x + E(:, i) + E(:, j)) ...
                             - f(x - E(:, i) + E(:, j)) ...
                             - f(x + E(:, i) - E(:, j)) ...
                             + f(x - E(:, i) - E(:, j)) ) * 0.25 / (E(i, i) * E(j, j));
        end
    end
end

