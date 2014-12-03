function [dydt]= ddefunc(t,y,Z,var)
dydt = zeros(size(y));
dydt(1:var.n) = y(var.n+1:2*var.n);

%get max displacements and forces from the previous step
dmax_0 = Z(2*var.n+1:3*var.n);
Fmax_0 = Z(3*var.n+1:4*var.n);

%relative displacements
d = zeros(var.n,1);
d(1:var.n-1) = -y(2:var.n) + y(1:var.n-1) + var.r(2:var.n) + var.r(1:var.n-1);
%d = d.*(d>=0);  %only allow positive displacement to create forces

%relative velocities
v = zeros(var.n,1);
v(1:var.n-1) = -y(var.n+2:2*var.n) + y(var.n+1:2*var.n-1);

%calculate contact forces
f = zeros(var.n,1);
for i = 1:var.n-1   
    f(i) = f_contact(i,d,v,var,dmax_0,Fmax_0);
end

%add contact forces to beads to the left
dydt(var.n+1:2*var.n) = dydt(var.n+1:2*var.n) - f./var.m(1:end);
%add contact forces to beads to the right
dydt(var.n+2:2*var.n) = dydt(var.n+2:2*var.n) + f(1:end-1)./var.m(2:end);

%add force from rigid wall
if isequal(var.wall,'yes')
    d(end) = (y(var.n) - var.xi(end))*((y(var.n) - var.xi(end))>0);
    f(var.n) = var.A(end)*d(end)^(3/2);
    dydt(2*var.n) = dydt(2*var.n) - var.A(end)*d(end)^(3/2)/var.m(end);
end

%add force from hopkinson transmission bar
if isequal(var.bar,'yes')
%     %Oxford
%     A = pi*(.01)^2;             %area of bar (10 mm radius)
%     c = 5082;                   %wave speed in bar (m/sec)
%     rho = 4430;                 %density of bar (kg/m^3)
    %Caltech
    A = pi*((3/8)*.0254)^2;             %area of bar (3/8 inch radius)
    c = 4865.5;                   %wave speed in bar (m/sec)
    rho = 8130;                 %density of bar (kg/m^3)
    if (y(var.n) - var.xi(end)) > 0
        dydt(2*var.n) = dydt(2*var.n) - A*rho*c*y(2*var.n)*(y(2*var.n)>0)/var.m(end);
    end
end

%add applied force on first bead
if isequal(var.applied_force,'yes')
    if t(end) < var.Fexternal(end,1)
        dydt(var.n+1) = dydt(var.n+1) + interp1(var.Fexternal(:,1),var.Fexternal(:,3),t(end))/var.m(1);
    end
end

%get slopes of maxima using previous step's values
d_dmax = dmax_0 - max(dmax_0,d);
d_Fmax = Fmax_0 - max(Fmax_0,f);
%apply slopes
dydt(2*var.n+1:3*var.n) = -d_dmax/var.dt;
dydt(3*var.n+1:4*var.n) = -d_Fmax/var.dt;

t
end

