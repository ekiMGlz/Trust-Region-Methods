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
grid on
hold on

%Graficar direcciones
quiver3(x_0(1),x_0(2),0,pC(1),pC(2),0,0,'r');
quiver3(x_0(1),x_0(2),0,pN(1),pN(2),0,0,'g');
quiver3(x_0(1),x_0(2),0,pDL(1),pDL(2),0,0,'b');

% Graficar modelo cuadratico en la region de confianza
fsurf(@(r,t) x_0(1)+r*cos(t), @(r,t)  x_0(2)+r*sin(t), @(r,t)  mc([r*cos(t);r*sin(t)]), [0,delta,-pi,pi], 'ShowContours', 'on')
%fcontour(@(r,t)  mc([r*cos(t);r*sin(t)]))

% Graficar la region de confianza
%viscircles(x_0',delta,'Color','y')

view(-129,28)
hold off


