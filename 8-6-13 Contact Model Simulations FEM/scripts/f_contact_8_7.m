function [f] = f_contact(i,d,v,var,dmax,Fmax)
% i - which contact location in the chain
% d - displacement of all contacts
% v - velocities of all particles
% var - struct with all property simulations
% dmax - maximum displacements experienced by each contact
% Fmax - maximum force experienced on each contact

% LOADING
d = max(0,d);   %only allow positive displacements, no adhesive forces
%strain rate dependence
if isequal(var.rate_dependent,'yes')
    sr=max(0,v(i))/(var.r(i)+var.r(i+1));
    sr_factor = max(1,(1+var.K_star(i)*log(sr/.001)));
else
    sr_factor = 1;
end
sy = var.sy_star(i)*sr_factor;

ratio = var.E_star(i)/sy;
%deltap = (.005017*ratio^(-1) - 7.231e-6)*(var.r_star(i)/(.0625*.0254));
deltap = (.0067*ratio.^(-1) - 3.082e-6)*(var.r_star(i)/(.0625*.0254));
c1 = -5.972*ratio^(-0.1637) + 5.487;
c2 = (-4.104e-6*ratio^(-1) + 3.999e-9)*(var.r_star(i)/(.0625*.0254))^2;

deltay=(1/4)*(var.r_star(i)/var.E_star(i)^2)*(pi*1.6*sy)^2;  %VON MISES DEFINITION
p0 = c1*sy;

if d(i) <= deltap && d(i) > deltay && isequal(var.plasticity,'yes')
    %intermediate region
    Fy = (4/3)*var.E_star(i)*sqrt(var.r_star(i))*deltay^(3/2);
    Fy_pr = 2*var.E_star(i)*sqrt(var.r_star(i)*deltay);
    a(1) = c2*p0*pi;
    a(2) = p0*2*pi*var.r_star(i);
    Fp_pr = a(2);
    Fp = a(1) + a(2)*deltap;
    
%     %want to solve for x that matches the slope of the interpolant to the
%     %slope of the linear region at deltap
%     parameter_minimizer = @(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr)...
%         x*(-1+deltap*x)^(-1)*(-1+deltay*x)^(-1)*(log(-1+deltap*x)+(-1)*...
%         log(-1+deltay*x))^(-1)*(-1*Fp+Fy+deltay*Fp*x+(-1)*deltap*Fy*x+...
%         Fp*(-1+deltay*x)*log(-1+deltap*x)+(Fp+(-1)*deltay*Fp*x)*log(-1+deltay*x)) - Fp_pr; 
%     %call non-linear solver to solve for x
%     myfun = @(x)parameter_minimizer(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr);
%     x0 = -100000;       %initial guess of x
%     options = optimoptions('fsolve','MaxFunEvals',500,'MaxIter',500,'TolFun',1e-6,'Display','off','Algorithm','levenberg-marquardt');
%     [x,~] = fsolve(myfun,x0,options); % Call solver
    
%     %Newton-Raphson Method to solve for interpolating function
%     func = @(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr)...
%         x*(-1+deltap*x)^(-1)*(-1+deltay*x)^(-1)*(log(-1+deltap*x)+(-1)*...
%         log(-1+deltay*x))^(-1)*(-1*Fp+Fy+deltay*Fp*x+(-1)*deltap*Fy*x+...
%         Fp*(-1+deltay*x)*log(-1+deltap*x)+(Fp+(-1)*deltay*Fp*x)*log(-1+deltay*x)) - Fp_pr; 
%     func_deriv = @(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr)...
%         x.^2.*((-1)+deltap.*x).^(-2).*((-1)+deltay.*x).^(-1).*(Fy+(-1).* ...
%         deltap.*Fy.*x+Fp.*((-1)+deltay.*x)).*(log((-1)+deltap.*x)+(-1).* ...
%         log((-1)+deltay.*x)).^(-1);
%     
%     maxsteps = 5000;
%     tol = 10^-1;
%     x = -100000;
%     for j = 1:maxsteps
%         x = x - func(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr)/...
%             func_deriv(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr);
%         if abs(func(x,deltay,deltap,Fy,Fp,Fy_pr,Fp_pr)) < tol
%             break
%         end
%     end
%     %j
%     
%     x = -10000000;
%     %x = real(x);
%     f = ((-1)+d(i).*x).*((-1)+deltap.*x).^(-1).*((-1)+deltay.*x).^(-1).*(log( ...
%         (-1)+deltap.*x)+(-1).*log((-1)+deltay.*x)).^(-1).*((Fy+(-1).* ...
%         deltap.*Fy.*x+Fp.*((-1)+deltay.*x)).*log((-1)+d(i).*x)+Fy.*((-1)+ ...
%         deltap.*x).*log((-1)+deltap.*x)+Fp.*(1+(-1).*deltay.*x).*log((-1)+ ...
%         deltay.*x));
%     f = real(f);
    
    f = d(i).*deltap.^(-1).*deltay.^(-1).*(log(deltap)+(-1).*log( ...
    deltay)).^(-1).*((deltay.*Fp+(-1).*deltap.*Fy).*log(d(i))+ ...
    deltap.*Fy.*log(deltap)+(-1).*deltay.*Fp.*log(deltay));
    %f = real(f);
elseif d(i) > deltap && isequal(var.plasticity,'yes')
    %plastic region
    a(1) = c2*p0*pi;
    a(2) = p0*2*pi*var.r_star(i);
    f = a(1) + a(2)*d(i);
else
    %elastic loading
    f = var.A(i)*d(i)^(3/2);
end


% UNLOADING
if f < Fmax(i)
    deltay=(1/4)*(var.r_star(i)/var.E_star(i)^2)*(pi*1.6*var.sy_star(i))^2;
    if isequal(var.plasticity,'yes') && dmax(i) > deltay
        %Plastic unloading
        Fy = (4/3)*var.E_star(i)*sqrt(var.r_star(i))*deltay^(3/2);
        Rp = 4*var.E_star(i)*((2*Fmax(i)+Fy)/(2*pi*1.6*var.sy(i)))^(3/2)/(3*Fmax(i));
        deltar = dmax(i)-(3*Fmax(i)/(4*var.E_star(i)*sqrt(Rp)))^(2/3);
        d_un = (d(i)-deltar)*((d(i)-deltar)>=0);

        f_un =(4/3)*var.E_star(i)*sqrt(Rp*(d_un).^3);
    else
        %elastic unloading
        f_un = var.A(i)*d(i)^(3/2);
    end
    f = min(Fmax(i),f_un);
end

end

