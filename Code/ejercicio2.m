beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
        
x_0 = [-3; -3];
x_opt = [3; 0.5];

[x, msg, i] = mRC1(beale, x_0, 10000)

norm(x-x_opt)
sprintf('RC2')
[x, msg, i] = mRC2(beale, x_0, 10000)
norm(x-x_opt)