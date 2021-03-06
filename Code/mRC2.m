function [x, msg] = mRC2(f, x0, itmax)
%MRC1 Trust Region Method using the Dogleg Point and parameters 
%    Min. Quality (eta) = 0.1,
%    Tolerance (tol) = 1e-5 
%    Max. Trust Region Radius (delta_max) = 1.5
%   Input:
%       f: (Anonymous function) Function to optimize
%       x0: (Vector) Initial point
%       itmax: (Possitive integer) maximum number of iterations
%   Output:
%       x: (Vector) Last approximation of a stationary point
%       msg: (String) 

    % Initialize some parameters used in the alg.
    eta = 0.1;
    tol = 1e-5;
    delta_max = 1.5;
    delta = 1;
    n = length(x0);
    
    % Obtain the first values of x_k, g_k and B_k
    i = 0;
    
    if iscolumn(x0)
        x_k = x0;
    else
        x_k = x0';
    end
    
    g_k = gradient(f, x_k);
    
    % Ensure B_k is s.p.d.
    B_k = hessian(f, x_k);
    l_k = min(eigs(B_k));
    if l_k <= 0
       B_k = B_k + (1e-12 - 1.125 * l_k)*speye(n); 
    end
    
    while norm(g_k, inf) > tol && i < itmax
        
        % Find the Dogleg Point
        p_k = pDogLeg(B_k, g_k, delta);
        
        % Calculate the quality of the approximation
        quality = ( f(x_k + p_k) - f(x_k) ) ...
                  / ( dot(g_k, p_k) + 0.5 * p_k'*B_k*p_k);
        
        % Adjust the trust region radius based on the quality
        if quality < 0.25
            delta = 0.25 * delta;
        elseif quality > 0.75 && (delta - norm(p_k))/delta < 0.01    
            delta = min(delta_max, 2*delta);
        end
        
        % If the quality is above the minimum quality, then advance towards
        % p_k, and calculate the new gradient and hessian for x_k+1
        if quality > eta
            
            x_k = x_k + p_k;
            g_k = gradient(f, x_k);
            
            B_k = hessian(f, x_k);
            l_k = min(eigs(B_k));
            if l_k <= 0
               B_k = B_k + (1e-12 - 1.125 * l_k)*speye(n); 
            end
            
            i = i + 1;
        end
        
    end

    x = x_k;
    if norm(g_k, inf) < tol
        % Verify if the Hessian is (semi)definite positive to see if we've found a local minimum
        msg = "Convergio en " + i + '" iteraciones- "';
        l_k = min(eigs(B_k));
        if l_k >= 0
            msg = msg + "Hessiano es positivo (semi)definido, por lo tanto, se encontro minimo local.";
        end
    else
        msg = "El metodo no convergio";
    end
    
end

