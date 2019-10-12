%% Exercise 1. Plot the quadratic model at x_0 and the Newton, Cauchy and dogleg direction

% Beale function 
beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
            
% intial point x_0, where the Hessian is positive definite
x_0 = [2; 0];

% Calculate the Hessian and the gradient approximations
B = hessian(beale, x_0);
g = gradient(beale, x_0);

% minumun eigenvalue is greater than zero
min(eigs(B))

% Quadratic model definition
fk = beale(x_0);
mc = @(p) fk + dot(g, p) + 0.5*p'*B*p;

% Calculate Newton's, Cauchy's and dogleg's direction
pN = -B\g;

% Define delta such that the Cauchy point is inside the trust region and
% the Newton point is on the outside.
delta = (norm(pN) + (g'*g/(g'*B*g)*norm(g)))*0.5;

pC = pCauchy(B,g,delta);
pDL = pDogLeg(B,g,delta);

%% Plot

% Plot the quadratic model in the trust region
fsurf(@(r,t) x_0(1)+r*cos(t), @(r,t)  x_0(2)+r*sin(t), @(r,t)  mc([r*cos(t);r*sin(t)]), [0,delta,-pi,pi], 'ShowContours', 'on')
colormap parula
hold on

% Plot directions
q_C = quiver3(x_0(1), x_0(2), 0, pC(1), pC(2), 0, 0, 'Color', '#99cc04',  'LineWidth', 1.5);
q_N = quiver3(x_0(1), x_0(2), 0, pN(1), pN(2), 0, 0, 'Color', '#a55af4', 'LineWidth', 1.5);
q_DL = quiver3(x_0(1), x_0(2), 0, pDL(1), pDL(2), 0, 0, 'Color', '#ff796c', 'LineWidth', 1.5);

% Plot the trust region's outline
RC = viscircles(x_0', delta,'LineStyle', ':', 'Color', '#3c4142', 'LineWidth',1);

% Plot dogleg's trajectory
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