clear

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


%x0 = [10^8;1.5;10^5;.5];  % Make a starting guess at the solution
x0 = [1;1];
options = optimoptions('fsolve','MaxFunEvals',10000,'MaxIter',10000); % Option to display output
[x,fval] = fsolve(@myfun,x0,options); % Call solver


d = linspace(0,2*10^-4,1000);
F = zeros(length(d),1);
for i = 1:length(d)
    if d(i) <= deltay
        F(i) = (4/3)*E*sqrt(r)*d(i)^(3/2);
    elseif d(i) <= deltap && d(i) > deltay
        A = x(1);
        b = x(2);
        F(i) = ((b+A*d(i))*((b*(Fp-Fy)-A*Fy*deltap+A*Fp*deltay)*...
            log(b+A*d(i))+Fy*(b+A*deltap)*log(b+A*deltap)-Fp*(b+A*deltay)*...
            log(b+A*deltay)))/((b+A*deltap)*(b+A*deltay)*(log(b+A*deltap)-...
            log(b+A*deltay)));
    elseif d(i) > deltap
        F(i) = a(1) + a(2)*d(i);
    end
end
plot(d,F,'LineWidth',1)
hold all
plot([deltay,deltap],[Fy,Fp],'x','LineWidth',2)
%axis([0 max(d)*1.05 0 max(F)*1.05])