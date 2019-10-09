function [p] = pDogLeg(B, g, delta)
%PCAUCHY Find the Cauchy Point for the trust region model of f
%   Input:
%       B: (Symmetric matrix) Approximated hessian of f at x_k
%       g: (Vector) Approximated gradient of f at x_k
%       delta: (Possitive real number) Trust region radius
%   Output:
%       p: Dogleg point for the model
%   
    
    % First, calculate the Cauchy Point and its norm
    pC = pCauchy(B, g, delta);
    norm_pC = norm(pC);
    
    % If the norm of pC is approximately equal to delta, p = pC
    if (delta - norm_pC)/delta < 0.01
        p = pC;
    else
        % Calculate Newton's Point, and check its norm to see if its within
        % the trust region
        pN = -B\g;
        norm_pN = norm(pN);
        
        if norm_pN < delta
            p = pN;
        else
            % Find alpha in (0, 1) such that norm(pC+alpha(pN - pC)) =
            % delta
            
            b = dot(pC, pN-pC);
            a = dot(pN - pC, pN - pC);
            alpha = (-b + sqrt(b^2 - (a)*(norm_pC^2 - delta^2))) / a;
            
            p = pC + alpha*(pN - pC);
        end
    end
end
