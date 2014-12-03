clear
clc
addpath('Scripts');
addpath('ExperimentalData');


%Experiment Geometry
%var.particles = 'sssssaaaaaa';
var.particles = 'ssssssssss';
% for i = 1:5
%     var.particles = strcat(var.particles,'ssaaaaaaaa');
% end
var.v0 = 0;
%var.striker_mass = .01;
var.wall = 'no';
var.bar = 'no';

%Model Options
var.plasticity = 'yes';
var.rate_dependent = 'yes';

%Initialization Function
var = initialize(var);

%Applied Force
var.applied_force = 'yes';
load('force17_5');
Fexternal(:,3) = smooth(Fexternal(:,3),100);
var.Fexternal = Fexternal;

%initial conditions (position, velocity, dmax, Fmax)
z0 = [var.xi; var.vi; zeros(var.n,1); zeros(var.n,1)];
t_final = 130e-6;
tol = 10^-8;
var.dt = tol;
delays = @(t,y) max(0,t-var.dt);
options = ddeset('RelTol',tol,'AbsTol',tol,'InitialStep',tol);

%call standard DDE solver
sol = ddesd(@(t,z,Z) ddefunc(t,z,Z,var),delays,z0,[0 t_final],options);
t = sol.x;
x = sol.y';

%get data from DDE solution
u = (x(:,1:var.n) - repmat(var.xi',length(t),1));
v = x(:,var.n+1:2*var.n);
dmax = x(:,2*var.n+1:3*var.n);
Fmax = x(:,3*var.n+1:4*var.n);

%forces were not saved, recalculate forces by calling f_contact again
d = -diff(u,1,2);
v_rel = -diff(v,1,2);
f = zeros(size(d));
for i = 1:length(t)
    for j = 1:var.n-1
        f(i,j) = f_contact(j,d(i,:),v_rel(i,:),var,dmax(i,:),Fmax(i,:));
    end
end

% hold all
% plot(d,f)

save('sim_data','t','u','v','d','f');
beep; pause(.07); beep; pause(.07); beep;
