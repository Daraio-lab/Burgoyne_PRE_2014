clear
clc
addpath('Scripts');
addpath('ExperimentalData');

% Experiment Geometry
var.particles = 'ss';
var.wall = 'no';
var.bar = 'no';

% Model Options
var.plasticity = 'yes';
var.rate_dependent = 'yes';

%Applied Force
var.applied_force = 'no';
load('force17_5');
Fexternal(:,3) = smooth(Fexternal(:,3),10);
var.Fexternal = Fexternal;

%vE = .045; %.045; %.022; %.375;
vE = .045; %aa = .5018; %tt = .2702; %ss = 2.494;
%vE =  1.0474647; %aa = 0.183950; %tt = 0.0951637; %ss = 1.0474647;
v0s = logspace(log10(vE),log10(10000*vE),30);
%v0s= linspace(.2,.5,10);
%v0s= linspace(.4828,.5172,30);
%v0s = 2.42;
CoR = [];
for i = 1:length(v0s)
    var.v0 = v0s(i);
    
    %Initialization Function
    var = initialize(var);
    
    %initial conditions (position, velocity, dmax, Fmax)
    %initial conditions (position, velocity, dmax, Fmax)
    z0 = [var.xi; var.vi; zeros(var.n,1); zeros(var.n,1)];
    t_final = 70e-6;
    tol = 10^-7;
    var.dt = tol;
    delays = @(t,y) max(0,t-var.dt);
    options = ddeset('RelTol',tol,'AbsTol',tol,'InitialStep',tol);
    
    %call standard DDE solver
    sol = ddesd(@(t,z,Z) ddefunc(t,z,Z,var),delays,z0,[0 t_final],options);
    t = sol.x;
    x = sol.y';
    u = (x(:,1:var.n) - repmat(var.xi',length(t),1));
    v = x(:,var.n+1:2*var.n);
    
    %for two spheres
    CoR = [CoR (v(end,2) - v(end,1))/var.v0];
    %for sphere onto flat
    %CoR = [CoR -v(end,1)/var.v0];
end
% plot(t,v(1,:),t,v(2,:))
% ann = ['v =',num2str(var.v0),'  CoR = ',num2str(CoR)];
% annotation('textbox', [.5 .7 .35 .1], 'String', ann,'FontSize',14);
% fid = fopen('data.txt', 'at');
% fprintf(fid, ann);
% fprintf(fid, '\n');
% fclose(fid);

%figure;
hold all
%plot(v0s,CoR)
loglog(v0s/vE,CoR)
axis([10^0 10^4 10^-1 10^0])


% load('ABAQUS_7_24/CoR_data.mat')
% plot(v,cor1,'x',v,cor2,'x',v,cor3,'x')


beep; pause(.07); beep; pause(.07); beep;
return

%%
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

%% Vy calculator
ratio = var.E_star(1)/var.sy_star(1);
deltap = (.005017*ratio^(-1) - 7.231e-6)*(var.r_star(1)/(.0625*.0254));
c1 = -5.972*ratio^(-0.1637) + 5.487;

Vy = (pi/(2*var.E_star(1)))^2*(8*pi*var.r_star(1)^3*(c1*var.sy_star(1))^5/(15*var.m(1)/2))^(1/2);

