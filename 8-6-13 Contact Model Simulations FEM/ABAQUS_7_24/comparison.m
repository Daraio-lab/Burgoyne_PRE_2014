%Force vs Time

figure;
% load('elastic50')
% color = 'm';
% plot(t,c1,color,t,c3,color,t,c5,color,t,c7,color,t,c9,color)
hold on

load('plastic50')
color = 'c';
plot(t,c1,color,t,c3,color,t,c5,color,t,c7,color,t,c9,color)


load('plasticJC50')
color = 'k';
plot(t,c1,color,t,c3,color,t,c5,color,t,c7,color,t,c9,color)


load('plasticJC50_b')
color = 'g';
plot(t,c1,color,t,c3,color,t,c5,color,t,c7,color,t,c9,color)

% Force Disp

figure;
hold all

% load('elastic50')
% plot(d1,c1)

% 
load('plastic50')
plot(d1,c1)


load('plasticJC50')
plot(d1,c1)


load('plasticJC50_b')
plot(d1,c1)

%%
load('plasticJC5')
plot(d1,c1)
hold on
load('plasticJC10')
plot(d1,c1)
load('plasticJC20')
plot(d1,c1)
load('plasticJC50')
plot(d1,c1)
