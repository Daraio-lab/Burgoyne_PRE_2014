function [f_corr] = disp_corr(t,f,Delta_x)
    A=0.57344; B=0.42790; C=18.6680; D=19.6640; E=-6.6213; F=2.3118;
    a = (3/8)*.0254;    %bar radius
    c0 = 4865.5;        %m/sec
    
    N = floor(length(f)/2);
    dt = mean(diff(t));
    w0 = 2*pi/t(end);
    
    A0 = (2/t(end))*sum(f*dt);
    for k = 1:N
        Ak(k) = (2/t(end))*f'*cos(k*w0*t)*dt;
        Bk(k) = (2/t(end))*f'*sin(k*w0*t)*dt;
    end
    
    f_corr = (A0/2)*ones(length(t),1);
    for k = 1:N
        func = @(x) ((x/c0) - A - B/(1 + C*(a/(2*pi*x/(k*w0)))^4 + D*(a/(2*pi*x/(k*w0)))^3 + E*(a/(2*pi*x/(k*w0)))^2 + F*(a/(2*pi*x/(k*w0)))^1.5));
        ck = fzero(func,c0);
        phi = k*w0*(Delta_x/ck - Delta_x/c0);
        f_corr = f_corr + Ak(k)*cos(k*w0*t - phi) + Bk(k)*sin(k*w0*t - phi);
    end
end