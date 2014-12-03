% Hayden Burgoyne Winter 2012-2013
%__________________________________________________________________________
% Hopkinson Bar data analysis script
%__________________________________________________________________________
% Assumes data with 20k data points
% Data requirements:  
%       time vector:                t
%       incident bar voltage:       V_inp
%       transmission bar volt.:     V_out
% Aux. Function Requirements
%       Dispersion correction:      disp_corr.m
%__________________________________________________________________________

%addpath('hopkinson_tests_6-6')
%load('hopkinson_tests_6-6/sample1')
%load('laser_tests_7-23/laser_test5')
t = t - t(1);
figure
plot(t,V_inp,t,V_out)   %plot to find timing
%%
%V_out = -V_out;

%Adjust timing for your data
t_start = .000093;      %time shortly before incident pulse arrives
t_final = .000485;      %cutoff time to avoid reflections


%subtract any initial offset, make sure that n final is before pulse
dt=mean(diff(t));   %sec, 
n_start = max(round(t_start/dt),1);
n_final = min(round(t_final/dt),length(t));
V_inp = V_inp - mean(V_inp(1:n_start));
V_out = V_out - mean(V_out(1:n_start));

%bar material properties
E=192*10^9;         %Pa
nu=0.31;
c0=4865.5;          %m/s, experimentally determined

Delta1 = (25.125 + 1*0.25)*0.0254;  %distance from first gauge to sample
Delta2 = (14.5 + 0*0.25)*0.0254; %distance to second gauge to sample
%Delta2 = (9.4 + 0.25)*0.0254;

%convert voltage to strains
Vex=15;             %external applied voltage
GF=2.12;           %calibration constant
sig_i= -E*2*V_inp./(Vex*GF*((1+nu)+V_inp*(nu-1)/Vex));
sig_t= -E*2*V_out./(Vex*GF*((1+nu)+V_out*(nu-1)/Vex));

%delays on ref/trans pulses from experimental setup
delay_r=(Delta1)/c0; %delay for reflected pulse
delay_t=(Delta2)/c0; %delay for transmitted pulse, inches to m

%translate reflected wave
sig_r=sig_i(round((2*delay_r)/dt):n_final);
t_r=(1:length(sig_r))*dt;

%translate transmitted wave
sig_t=sig_t(round((delay_r+delay_t)/dt):n_final);
t_t=(1:length(sig_t))*dt;

t_final=min(t_r(end),t_final);  %has the longest to travel, will always be the final
n_final=min(n_final,length(t_r));
T = t(n_start:n_final);

%dispersion correction
%Delta1 = 0; Delta2 = 0;        %can uncomment these to remove disp. corr.
sig_i_DC = disp_corr(T,sig_i(n_start:n_final),Delta1);
sig_r_DC = disp_corr(T,sig_r(n_start:n_final),-Delta1);
sig_t_DC = disp_corr(T,sig_t(n_start:n_final),-Delta2);

sig_i_DC = smooth(sig_i_DC,1);
sig_r_DC = smooth(sig_r_DC,1);
sig_t_DC = smooth(sig_t_DC,1);

%displacement calculations
back_disp=c0*cumtrapz(sig_t_DC)*dt/E;
front_disp=c0*cumtrapz(sig_i_DC-sig_r_DC)*dt/E;
disp=-(front_disp+back_disp);
v = diff(disp)/dt;
SR = v/(.0254/4);

%stress calculations
A=pi*(9.525*10^-3)^2;   %area of transmission bar
force=sig_t_DC*A;       %force determined from just transmission bar
front_force=(sig_i_DC+sig_r_DC)*A;
force_avg = (force-front_force)/2;   %avg of front & back force

%plot(disp,force_avg)
return

%% Various Plots

%plot(t(1:n_final),sig_i(1:n_final),t(1:n_final),sig_t(1:n_final),t(1:n_final),sig_r(1:n_final))
%plot(t(n_start:n_final),sig_i_DC,t(n_start:n_final),sig_t_DC,t(n_start:n_final),sig_r_DC)
plot(disp,force)
%plot(disp,force_avg)
hold all
%plot(t(1:length(sig_r)),sig_r,t(1:n_final),sig_i(1:n_final),t(1:length(sig_t)),sig_t)

%%
plot(t_r,-front_force)
hold on
plot(t_t,force)
plot(t_r,force_avg,'r')
%%
fd1 = [disp, force];
%save('fd1.mat','fd1')
%%
t1 = 0:dt:dt*(length(fd1)-1);
Fexternal= [t1',fd1];
%%
plot(0:dt:dt*(length(disp)-2),diff(disp)/dt/(.0254/4))

%% Output front velocity profile
front_v = -smooth(diff(front_disp)/dt,20);
front_v = front_v - front_v(1);
figure
plot(T(1:end-1)-T(1),front_v);
hold all
%%
Vexternal = [(0:dt:dt*(length(front_v)-1))',front_v];
save('Vext_exp15.mat','Vexternal')

%% Output back velocity profile
back_v = smooth(diff(back_disp)/dt,20);
back_v = back_v - back_v(1);
plot(T(1:end-1)-T(1),back_v);
%% Plot transmitted force
hold all
a = dt/5e-6;
plot(T-T(1),filter(a, [1 a-1],force),'LineWidth',1.5)
