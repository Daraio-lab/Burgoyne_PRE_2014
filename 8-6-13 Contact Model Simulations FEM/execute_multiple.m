clear
clc
addpath('Scripts');
addpath('ExperimentalData');

% Experiment Geometry
var.particles = 'ss';
var.v0 = 0;
%var.striker_mass = .05;
var.wall = 'no';
var.bar = 'yes';

% Model Options
var.plasticity = 'yes';
var.rate_dependent = 'yes';

%Initialization Function
var = initialize(var);

figure;
Fs = {'force7_5','force10','force10b'};%,'force12_5'}; %,'force17_5'};
%Fs = {'force7_5','force10','force10b','force12_5'};
Fs = {'force10_oxford','force15_oxford','force17_oxford'};
%Fs = {'s302_10','s302_12','s302_15','s302_17'};
%Fs = {'force10'};
%Fs = {'al_12','al_15','al_17'};
for k = 1:length(Fs)
    clearvars -except Fs var k
    %Applied Force
    var.applied_force = 'yes';
    load(Fs{k});
    Fexternal(:,3) = smooth(Fexternal(:,3),100);
    var.Fexternal = Fexternal;
    
    %initial conditions (position, velocity, dmax, Fmax)
    z0 = [var.xi; var.vi; zeros(var.n,1); zeros(var.n,1)];
    t_final = 200e-6;
    var.dt = t_final/2499;
    tspan = 0:var.dt:t_final;
    delays = @(t,y) max(0,t-var.dt);
    options = ddeset;%('RelTol',1e-7,'AbsTol',1e-7);
    
    sol = ddesd(@(t,z,Z) ddefunc(t,z,Z,var),delays,z0,tspan,options);
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
    
    %  Plot Force-Disp on bead N
    N=1;
    plot(d(:,N),f(:,N),'b')
    hold on
    
    % Plot external loading
    plot(Fexternal(:,2),Fexternal(:,3),'r')
    drawnow
end

beep; pause(.07); beep; pause(.07); beep;
