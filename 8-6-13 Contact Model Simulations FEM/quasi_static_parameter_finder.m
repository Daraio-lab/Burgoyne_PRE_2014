clearvars -except test_data rsave deltapsave c1save c2save
%close
load('ABAQUS_7_15/sy1000')
test_data(:,1) = -cdata(:,end);
contacts = cdata(:,2:44);
stresses = cdata(:,45:87);
test_data(:,2) = sum(contacts,2);

smoothF = smooth(test_data(:,2),50);
d = test_data(:,1);


E1 = 100E9;
r1 = (1/8)*(.0254);
E2 = 100E9;
r2 = (1/8)*(.0254);
v = .3;
%sigy = 1896E6;

%two hemispheres
E = ((1-v^2)/E1 + (1-v^2)/E2)^-1;
r = ((1/r1)+(1/r2))^-1;

figure; plot(d,smoothF); hold all

%Fit model to find min error for a given deltaP
[~,Nmax] = max(test_data(:,1));
ns = round(.05*Nmax):1:round(.75*Nmax);
error = [];
for j = ns
    d_AB = test_data(j:Nmax,1);
    F_AB = smoothF(j:Nmax);
    a = [ones(length(d_AB),1),d_AB]\F_AB;
    
    c1 = a(2)/(sigy*pi*2*r);                     %c(1) = p0/sigy
    c2 = a(1)/(c1*sigy*pi);
    
    deltap = d_AB(1);
    Fp = a(1) + a(2)*deltap;
    Fp_pr = a(2);
   
    deltay=(1/4)*(r/E^2)*(pi*1.6*sigy)^2;
    
    Fy = (4/3)*E*sqrt(r)*deltay^(3/2);
    Fy_pr = 2*E*sqrt(r*deltay);
    
    myfun = @(x)parameter_minimizer(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr);
    x0 = -100000;
    options = optimoptions('fsolve','MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-8,'Display','off'); % Option to display output
    [x,fval] = fsolve(myfun,x0,options); % Call solver
    
    d = test_data(:,1);
    F = zeros(length(d),1);
    for i = 1:length(d)
        if d(i) <= deltay
            F(i) = (4/3)*E*sqrt(r)*d(i)^(3/2);
        elseif d(i) < deltap && d(i) > deltay
            x = real(x);
            F(i) = ((-1)+d(i).*x).*((-1)+deltap.*x).^(-1).*((-1)+deltay.*x).^(-1).*(log( ...
                (-1)+deltap.*x)+(-1).*log((-1)+deltay.*x)).^(-1).*((Fy+(-1).* ...
                deltap.*Fy.*x+Fp.*((-1)+deltay.*x)).*log((-1)+d(i).*x)+Fy.*((-1)+ ...
                deltap.*x).*log((-1)+deltap.*x)+Fp.*(1+(-1).*deltay.*x).*log((-1)+ ...
                deltay.*x));
        elseif d(i) >= deltap
            F(i) = a(1) + a(2)*d(i);
        end
    end
    error = [error norm(abs(F-smoothF))];   
end


%Evaluate for the min error deltaP
[~,n] = min(error);
d_AB = test_data(ns(n):Nmax,1);
F_AB = smoothF(ns(n):Nmax);
a = [ones(length(d_AB),1),d_AB]\F_AB;

c1 = a(2)/(sigy*pi*2*r);                     %c(1) = p0/sigy
c2 = a(1)/(c1*sigy*pi);    

deltap = d_AB(1);
Fp = a(1) + a(2)*deltap;
Fp_pr = a(2);

deltay=(1/4)*(r/E^2)*(pi*1.6*sigy)^2;
 
Fy = (4/3)*E*sqrt(r)*deltay^(3/2);
Fy_pr = 2*E*sqrt(r*deltay);

myfun = @(x)parameter_minimizer(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr);
x0 = -100000;
options = optimoptions('fsolve','MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-8); % Option to display output
[x,fval] = fsolve(myfun,x0,options); % Call solver

d = test_data(:,1);
F = zeros(length(d),1);
for i = 1:length(d)
    if d(i) <= deltay
        F(i) = (4/3)*E*sqrt(r)*d(i)^(3/2);
    elseif d(i) <= deltap && d(i) > deltay
        x = real(x);
        F(i) = ((-1)+d(i).*x).*((-1)+deltap.*x).^(-1).*((-1)+deltay.*x).^(-1).*(log( ...
            (-1)+deltap.*x)+(-1).*log((-1)+deltay.*x)).^(-1).*((Fy+(-1).* ...
            deltap.*Fy.*x+Fp.*((-1)+deltay.*x)).*log((-1)+d(i).*x)+Fy.*((-1)+ ...
            deltap.*x).*log((-1)+deltap.*x)+Fp.*(1+(-1).*deltay.*x).*log((-1)+ ...
            deltay.*x));
    elseif d(i) > deltap
        F(i) = a(1) + a(2)*d(i);
    end
end
plot(d,F)
plot([deltay,deltap],[Fy,Fp],'x','MarkerSize',10)

rsave = [rsave; E/sigy];
deltapsave = [deltapsave; deltap];
c1save = [c1save; c1];
c2save = [c2save; c2];
beep; pause(.07); beep; pause(.07); beep;


return
%%
ratio = rsave;
%deltap = (.005017*ratio.^(-1) - 7.231e-6);
deltap = (.0067*ratio.^(-1) - 3.082e-6);
c1 = -5.972*ratio.^(-0.1637) + 5.487;
c2 = (-4.104e-6*ratio.^(-1) + 3.999e-9);

figure
subplot(3,1,1); plot(rsave,deltap,rsave,deltapsave,'x')
subplot(3,1,2); plot(rsave,c1,rsave,c1save,'x')
subplot(3,1,3); plot(rsave,c2,rsave,c2save,'x')

