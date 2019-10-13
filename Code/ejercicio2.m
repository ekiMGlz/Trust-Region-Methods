%% Excercise 2. Number of iterations and error.

%We define the Beale function
beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
        
%We define x_0 as the point where ||x_0 - x_opt|| > delta
x_0 = [-3; -3];
x_opt = [3; 0.5];

%Calculate the number of itrations, the solution and the result of
%convergence of the first method starting in x_0.
[x, msg, i] = mRC1(beale, x_0, 10000);

%Calculate the error between the result of the first method and the optima
%solution.
e = norm(x-x_opt);
sprintf("[Cauchy] Valor de x: ( %g, %g), el error de aproximación es: %d", x(1),x(2), e)

%Calculate the number of itrations, the solution and the result of
%convergence of the second method starting in x_0.
[x, msg, i] = mRC2(beale, x_0, 10000);

%Calculate the error between the result of the second method and the optima
%solution.
e = norm(x-x_opt);
sprintf("[DogLeg] Valor de x: ( %g, %g), el error de aproximación es: %d", x(1),x(2), e)