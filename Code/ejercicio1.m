%% Ejercicio 1. Graficar el modelo cuadratico en x_0 y las direcciones Newton, Cauchy, dogleg

% Definimos la funcion de Beale
beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
            
% x_0 donde la matriz Hessiana es definida positiva
x_0 = [2; 0];

% Calculamos las aproximaciones del Hessiano y del gradiente
B = hessian(beale, x_0);
g = gradient(beale, x_0);

% Eigenvalor minimo es mayor a cero
min(eigs(B))

% Definicion del modelo cuadratico
fk = beale(x_0);
mc = @(p) fk + dot(g, p) + 0.5*p'*B*p;

% Calculamos las direcciones del punto Newton, Cauchy y dogleg
pN = -B\g;

% Definimos delta tal que el punto de Cauchy quede dentro de la region de confianza
% y el punto de Newton quede fuera.
delta = (norm(pN) + (g'*g/(g'*B*g)*norm(g)))*0.5;

pC = pCauchy(B,g,delta);
pDL = pDogLeg(B,g,delta);

%% Grafica

% Graficar modelo cuadratico en la region de confianza
fsurf(@(r,t) x_0(1)+r*cos(t), @(r,t)  x_0(2)+r*sin(t), @(r,t)  mc([r*cos(t);r*sin(t)]), [0,delta,-pi,pi], 'ShowContours', 'on')
colormap parula
hold on

% Graficar direcciones
q_C = quiver3(x_0(1), x_0(2), 0, pC(1), pC(2), 0, 0, 'Color', '#99cc04',  'LineWidth', 1.5);
q_N = quiver3(x_0(1), x_0(2), 0, pN(1), pN(2), 0, 0, 'Color', '#a55af4', 'LineWidth', 1.5);
q_DL = quiver3(x_0(1), x_0(2), 0, pDL(1), pDL(2), 0, 0, 'Color', '#ff796c', 'LineWidth', 1.5);

% Graficar la region de confianza
RC = viscircles(x_0', delta,'LineStyle', ':', 'Color', '#3c4142', 'LineWidth',1);

% Graficar trayectoria dogleg
x_N = x_0 + pN;
x_C = x_0 + pC;
D = x_N - x_C;

plot3([x_C(1),x_N(1)], [x_C(2), x_N(2)],[0,0],'LineStyle', '--','Color','#a5a391','LineWidth',0.6)

x_DL = x_0 + pDL;

plot3(x_0(1),x_0(2),0,'o','MarkerFaceColor','#f9bc08','MarkerEdgeColor','#f9bc08')
legend([q_C, q_N, q_DL, RC], {'Dirección Cauchy', 'Dirección Newton', 'Dirección Dogleg', 'Región de confianza'});

hold off
grid on

view(-147,21)