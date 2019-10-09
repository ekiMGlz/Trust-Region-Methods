function [pC] = pCauchy(B, g, delta)
%PCAUCHY Find the Cauchy Point for the trust region model of f
%   Input:
%       B: (Symmetric matrix) Approximated hessian of f at x_k
%       g: (Vector) Approximated gradient of f at x_k
%       delta: (Possitive real number) Trust region radius
%   Output:
%       pC: Cauchy point for the model
%   
    % Use pk = g/norm(g) as a first point
    g_norm = norm(g);
    p = - g/g_norm;
    pBp = p'*B*p;
    
    alpha = delta;
    
    % If pBp <= 0, then the minimum is at the edge of the trust region
    % Otherwise, check if the minimum is located inside the trust region
    if pBp > 0
        alpha = min(delta, g_norm/pBp);
    end
    
    pC = 0.99 * alpha * p;
end