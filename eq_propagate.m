%----------- Chao Huang --------------------

function dfdt = eq_propagate(t,initial_condition)

  
    u=3.986005*10^5;
    Re=6378.14;         % Earth Radius
    
    dfdt = zeros(size(initial_condition));
    
    %x=x
    %x1=dx/dt
    %dx1/dt = d^2(x)/dt
    
    x1 = initial_condition(4); % x1 = dx/dt
    y1 = initial_condition(5); % y1 = dy/dt
    z1 = initial_condition(6); % z1 = dz/dt
    x  = initial_condition(1); % x
    y  = initial_condition(2); % y
    z  = initial_condition(3); % z

    m  = initial_condition(7); % mass
    
    r=sqrt(x^2+y^2+z^2);

        
        dfdt(1) = x1;
        dfdt(2) = y1;
        dfdt(3) = z1;
        dfdt(4) = -u*x/r^3;
        dfdt(5) = -u*y/r^3;
        dfdt(6) = -u*z/r^3;
        dfdt(7) = 0;
end

