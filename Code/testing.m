A = pascal(4);
b = -ones(4, 1);
f = @(x) 0.5 * x'*A*x + dot(b, x) + 1;
g_f = @(x) 0.5*(A + A')*x + b;
h_f = 0.5*(A+A');