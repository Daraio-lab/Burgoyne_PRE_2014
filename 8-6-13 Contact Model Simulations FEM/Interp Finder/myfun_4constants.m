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

scale = 1;

% F = [x(1)*deltay^x(2) + x(3)*deltay^x(4) - Fy/scale; x(1)*x(2)*deltay^(x(2)-1) + x(3)*x(4)*deltay^(x(4)-1) - Fy_pr/scale;...
%     x(1)*deltap^x(2) + x(3)*deltap^x(4) - Fp/scale; x(1)*x(2)*deltap^(x(2)-1) + x(3)*x(4)*deltap^(x(4)-1) - Fp_pr/scale];

% F = [x(1)*x(2)^deltay + x(3)*x(4)^deltay - Fy/scale; x(1)*deltay*x(2)^(deltay-1) + x(3)*deltay*x(4)^(deltay-1) - Fy_pr/scale;...
%     x(1)*x(2)^deltap + x(3)*x(4)^deltap - Fp/scale; x(1)*deltap*x(2)^(deltap-1) + x(3)*deltap*x(4)^(deltap-1) - Fp_pr/scale];

% F = [x(1)*exp(x(2)*deltay) + x(3)*exp(x(4)*deltay) - Fy/scale; x(1)*x(2)*exp(x(2)*deltay) + x(3)*x(4)*exp(x(4)*deltay) - Fy_pr/scale;...
%     x(1)*exp(x(2)*deltap) + x(3)*exp(x(4)*deltap) - Fp/scale; x(1)*x(2)*exp(x(2)*deltap) + x(3)*x(4)*exp(x(4)*deltap) - Fp_pr/scale];

% F = [x(1) + x(2)*deltay + x(3)*deltay^2 + x(4)*deltay^3 - Fy/scale; x(2) + 2*x(3)*deltay + 3*x(4)*deltay^2 - Fy_pr/scale;...
%     x(1) + x(2)*deltap + x(3)*deltap^2 + x(4)*deltap^3 - Fp/scale; x(2) + 2*x(3)*deltap + 3*x(4)*deltap^2 - Fp_pr/scale];


F = Fy*[(x(1)*deltay+x(2))*(x(3)+x(4)*log(x(1)*deltay+x(2))) - Fy/scale; x(1)*x(4)+x(1)*(x(3)+x(4)*log(x(1)*deltay+x(2))) - Fy_pr/scale;...
    (x(1)*deltap+x(2))*(x(3)+x(4)*log(x(1)*deltap+x(2))) - Fp/scale; x(1)*x(4)+x(1)*(x(3)+x(4)*log(x(1)*deltap+x(2))) - Fp_pr/scale];


end




