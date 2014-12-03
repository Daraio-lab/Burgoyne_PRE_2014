clear
clc
addpath('Scripts');
addpath('ExperimentalData');
addpath('ABAQUS_8_2');
load('steel440JC_0_001')
%Experiment Geometry
var.particles = 'ss';
var.v0 = 0; %epsilon_dot*.00635;             %relative velocity, for fixed strain rate
var.wall = 'no';
var.bar = 'no';

%Model Options
var.plasticity = 'yes';
var.rate_dependent = 'yes';

%Initialization Function
var = initialize(var);


dmax = 25*10^-5;
d = [linspace(0,dmax,1000)];%,linspace(dmax,0,100)];

f = 0;
for i = 1:length(d)
    f = [f f_contact(1,d(i),var.v0(1),var,max(d(1:i)),max(f))];
end
f = f(2:end);

hold on
%figure
plot(d,f,test_data(:,1),test_data(:,2))
% 
% hold all
% plot(d,f)

