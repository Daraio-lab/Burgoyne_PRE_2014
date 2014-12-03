function [F] = myfun(x)

E1 = 200E9;
r1 = (2/16)*(.0254);
E2 = 200E9;
r2 = (2/16)*(.0254);
v = .3;
sigy = 1896E6;

%two hemispheres
E = ((1-v^2)/E1 + (1-v^2)/E2)^-1;
r = ((1/r1)+(1/r2))^-1;

ratio = E/sigy;
deltap = (.005017*ratio^(-1) - 7.231e-6);
c1 = -5.972*ratio^(-0.1637) + 5.487;
c2 = (-4.104e-6*ratio^(-1) + 3.999e-9);

%scale for different radius
deltap = deltap*(r/(0.0625*.0254));
c2 = c2*(r/(0.0625*.0254))^2;

deltay = deltap/c1^2;
p0 = c1*sigy;

Fy = (4/3)*E*sqrt(r)*deltay^(3/2);
Fy_pr = 2*E*sqrt(r*deltay);
a(1) = c2*p0*pi;
a(2) = p0*2*pi*r;
Fp_pr = a(2);
Fp = a(1) + a(2)*deltap;

a = x(1);
b = x(2);

F = [(a*(b*Fp-b*Fy-a*Fy*deltap+a*Fp*deltay+(b*(Fp-Fy)-a*Fy*deltap+a*Fp*deltay)*log(b+a*deltay)+Fy*(b+a*deltap)*log(b+a*deltap)-b*Fp*log(b+a*deltay)-a*Fp*deltay*log(b+a*deltay)))/((b+a*deltap)*(b+a*deltay)*(log(b+a*deltap)-log(b+a*deltay))) - Fy_pr;...
    (a*(b*Fp-b*Fy-a*Fy*deltap+a*Fp*deltay+(b*(Fp-Fy)-a*Fy*deltap+a*Fp*deltay)*log(b+a*deltap)+Fy*(b+a*deltap)*log(b+a*deltap)-b*Fp*log(b+a*deltay)-a*Fp*deltay*log(b+a*deltay)))/((b+a*deltap)*(b+a*deltay)*(log(b+a*deltap)-log(b+a*deltay))) - Fp_pr]/1e6;

end




