%----------- Chao Huang --------------------

function dfdt = eq_firing_lae(t,initial_condition,declination,ap_per)

    J2=1082.6*10^-6;
    u=3.986005*10^5;        % [km^3/s^2] 
    Re=6378.14;             % [km] Earth Radius
    f=490;                  % [N] engine thrust 
    Ispe=312*0.99;          % [s] specific impulse engine
    g0=9.8066;              % [m/s] gravity
    
    
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
    
    %right_asc_firing = comp_long(x,y,z,declination,ap_per);
    right_asc_firing = comp_long_new(x,y,z,declination,ap_per);
    
    
   %plot(r,right_asc_firing*180/pi,'xb')
  
%    hold on;   
%   xlabel('r [km]');
%   ylabel('x*cos(right_asc_firing)*cos(declination)+y*sin(right_asc_firing)*cos(declination)+z*sin(declination)');
%   ylim([-0.1 0.1]);
%   plot(r,x*cos(right_asc_firing)*cos(declination)+y*sin(right_asc_firing)*cos(declination)+z*sin(declination),'xb') % must be 0 ABcos(ab)

  
    
    
   
    


        dfdt(1) = x1;
        dfdt(2) = y1;
        dfdt(3) = z1;
        dfdt(4) = -u*x/r^3 + (f/m)/1000*cos(right_asc_firing)*cos(declination); %km 
        dfdt(5) = -u*y/r^3 + (f/m)/1000*sin(right_asc_firing)*cos(declination); %km
        dfdt(6) = -u*z/r^3 + (f/m)/1000*sin(declination); %km 
        dfdt(7) = -f/(Ispe*g0);
end
