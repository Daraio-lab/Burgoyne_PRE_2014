function [var] = initialize(var)

    var.n = length(var.particles);
    nu = .3;
    %radii and mass of each type of particle
    %Steel 440 (s)
    r_s = (1/8)*.0254;
    m_s = 7650*(4/3)*pi*r_s^3;
    E_s = 200e9;
    sy_s = 1896e6;
    %K_s = .05;  %FEM
    K_s = .035; 
    %Steel 302 (t)
    r_t = (1/8)*.0254;
    m_t = 7900*(4/3)*pi*r_t^3;
    E_t = 200e9;
    sy_t = 600e6;
    K_t = .025;
    %Aluminum (a)
    r_a = (1/8)*.0254;
    m_a = 2700*(4/3)*pi*r_a^3;
    E_a = 76e9;
    sy_a = 276e6;
    K_a = 0.005; 
    %Aluminum Oxide (o)
    r_o = (1/8)*.0254;
    m_o = 3950*(4/3)*pi*r_o^3;
    E_o = 370e9;
    sy_o = 380e6;
    K_o = 0.005; 
    %Aluminum Oxide (f)
    r_f = 10000;
    m_f = 2700*(4/3)*pi*r_f^3;
    E_f = 76e9;
    sy_f = 276e6;
    K_f = 0; 

    
    %Particle Properties (particle mass, radii, Young's moduli)
    for i = 1:var.n
        if isequal(var.particles(i),'s')
            var.r(i,1) = r_s;
            var.m(i,1) = m_s;
            var.E(i,1) = E_s;
            var.sy(i,1) = sy_s;
            var.K(i,1) = K_s;
        elseif isequal(var.particles(i),'t')
            var.r(i,1) = r_t;
            var.m(i,1) = m_t;
            var.E(i,1) = E_t;
            var.sy(i,1) = sy_t;
            var.K(i,1) = K_t;
        elseif isequal(var.particles(i),'a')
            var.r(i,1) = r_a;
            var.m(i,1) = m_a;
            var.E(i,1) = E_a;
            var.sy(i,1) = sy_a;
            var.K(i,1) = K_a;
        elseif isequal(var.particles(i),'o')
            var.r(i,1) = r_o;
            var.m(i,1) = m_o;
            var.E(i,1) = E_o;
            var.sy(i,1) = sy_o;
            var.K(i,1) = K_o;
        elseif isequal(var.particles(i),'f')
            var.r(i,1) = r_f;
            var.m(i,1) = m_f;
            var.E(i,1) = E_f;
            var.sy(i,1) = sy_f;
            var.K(i,1) = K_f;
        else
            disp('Error defining masses or radii')
        end
    end
    %if striker has different mass, change the first particle's mass
    if isfield(var,'striker_mass') == 1
        var.m(1) = var.striker_mass;
    end
    
    %Contact Properties: (Hertzian coefficients, empirical contact parameters)
    for i = 2:var.n
        %Hertzian coefficients
        var.E_star(i-1,1) = ((1-nu^2)/var.E(i-1) + (1-nu^2)/var.E(i))^-1;
        var.r_star(i-1,1) = ((1/var.r(i-1))+(1/var.r(i)))^-1;
        var.A(i-1,1) = (4/3)*var.E_star(i-1)*sqrt(var.r_star(i-1));
        
        var.sy_star(i-1,1) = min(var.sy(i-1),var.sy(i));
        var.K_star(i-1,1) = min(var.K(i-1),var.K(i));
    end
    %contact properties for bead-wall, rigid flat wall, elastic-only bead
    var.E_star(var.n) = ((1-nu^2)/var.E(end))^-1;
    var.r_star(var.n) = ((1/var.r(end)))^-1;
    var.A(var.n) = (4/3)*var.E_star(end)*sqrt(var.r_star(end));
    
    %create initial position and velocity vectors
    var.xi = zeros(var.n,1);
    for i = 1:var.n-1
        var.xi(i+1) = var.xi(i) + var.r(i) + var.r(i+1);
    end
    var.vi = zeros(var.n,1);
    var.vi(1) = var.v0;
    
    global g
    g.d = [];
    g.f = [];
    g.t = []; 
    
    
end

