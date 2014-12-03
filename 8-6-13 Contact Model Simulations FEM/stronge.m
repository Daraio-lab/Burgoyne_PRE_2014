clearvars -except test_data
%close
load('abaqus_steel440')
smoothF = smooth(test_data(:,2),10);
d = test_data(:,1);

% %Steel 440
E1 = 200E9;
r1 = (1/8)*(.0254);
E2 = 200E9;
r2 = (1/8)*(.0254);
v = .3;
sigy = 1896E6;

%two hemispheres
E = ((1-v^2)/E1 + (1-v^2)/E2)^-1;
r = ((1/r1)+(1/r2))^-1;

figure; plot(d,smoothF); hold all

deltay=(1/4)*(r/E^2)*(pi*1.6*sigy)^2; 
Fy = (4/3)*E*sqrt(r)*deltay^(3/2);
Fy_pr = 2*E*sqrt(r*deltay);

th = 2.8;

Fp = 424*Fy;
deltap = 84*deltay;


d = test_data(:,1);
F = zeros(length(d),1);
for i = 1:length(d)
    if d(i) <= deltay
        F(i) = (4/3)*E*sqrt(r)*d(i)^(3/2);
    elseif d(i) <= deltap && d(i) > deltay
        F(i) = Fy*(2*d(i)/deltay - 1)*(1+(3*th)^-1*log(2*d(i)/deltay - 1));
    elseif d(i) > deltap
        F(i) = Fy*2.8*(2*d(i)/deltay -1)/th;
    end
end
plot(d,F)
plot([deltay,deltap],[Fy,Fp],'x','MarkerSize',10)
