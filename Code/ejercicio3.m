beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
        
x_0 = [-3; -3];


x_opt = [3; 0.5];
xk = x_0;

hold on
fcontour(@(x,y) beale([x;y]), 'LevelList', [beale(x_0), beale(x_0)/2, beale(x_0)/4, beale(x_0)/8, beale(x_0)/16, beale(x_0)/32])
colorbar

plot(x_0(1),x_0(2),'bd')
plot(x_opt(1),x_opt(2),'yd')

cont = 1;
while norm(gradient(beale,xk), inf) > 1e-5 
    [xk, ~,~] = mRC1(beale, xk, 1);
    plot(xk(1), xk(2),'ro')
    cont = cont + 1;
end
hold off
cont

xk = x_0;
hold on
cont = 1;
while norm(gradient(beale,xk), inf) > 1e-5 
    [xk, ~,~] = mRC2(beale, xk, 1);
    plot(xk(1), xk(2),'gx')
    cont = cont + 1;
end
hold off
cont

