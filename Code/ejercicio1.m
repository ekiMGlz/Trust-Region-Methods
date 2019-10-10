beale = @(x)  (1.5 - x(1) + x(1) * x(2))^2    ...
            + (2.25 - x(1) + x(1) * x(2)^2)^2 ...
            + (2.625 - x(1) + x(1) * x(2)^3)^2;
        
x_0 = [2; 0];

B = hessian(beale, x_0);

g = gradient(beale, x_0);

min(eigs(B))


fk = beale(x_0);

mc = @(p) fk + dot(g, p) + 0.5*p'*B*p;
pN = -B\g;
delta = (norm(pN) + (g'*g/(g'*B*g)*norm(g)))*0.5;
grid on
%fsurf(@(x,y) mc([x;y]), [x_0(1)-Delta, x_0(1)+Delta, x_0(2)-Delta, x_0(2)+Delta], 'ShowContours', 'on')
hold on

pC = pCauchy(B,g,delta);





pDL = pDogLeg(B,g,delta);

Delta = 1.1*delta;

quiver3(x_0(1),x_0(2),0,pC(1),pC(2),0,0,'r');
quiver3(x_0(1),x_0(2),0,pN(1),pN(2),0,0,'g');
quiver3(x_0(1),x_0(2),0,pDL(1),pDL(2),0,0,'b');

fmesh(@(r,t) x_0(1)+r*cos(t), @(r,t)  x_0(2)+r*sin(t), @(r,t)  mc([r*cos(t);r*sin(t)]), [0,delta,-pi,pi], 'ShowContours', 'on' )
viscircles(x_0',delta,'Color','y')
%fmesh(@(r,t) x_0(1)+r*cos(t), @(r,t)  x_0(2)+r*sin(t), @(r,t)  0, [delta,delta,-pi,pi], 'ShowContours', 'on' )
view(-129,28)
hold off


