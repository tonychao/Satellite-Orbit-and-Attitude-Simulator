%---------Chao Huang ---- July- 2013 ---
function [ right_asc ] = comp_long_new( x, y, z, declination,ap_per)
    %   x,y,z: represent the position of satellite in the second equatorial coordination system
    %   declination: represent the declination of engine force, which is got from your optimization region of force direction delta
    %   fm: it only has two value:0 or 1,and we can use 0 to denote it is an apogee engine fire, and use 1 to denote this is a perigee maneuvor firing.
	%   m represent the longitude of engine force that we want to get
    %   subLongi;           //longitude of sub-satellite 
	%   subLati;            //latitude of sub-satellite
    %   you need to get this two intermediate variables by yourself,using x,y,z,
    
    %syms a b c x y z;    
    %solve(x*cos(a)*cos(b) + y*sin(a)*cos(b)+z*sin(b),'a')
    
 

     b=declination;
     if ap_per == 0 % 0 = apogee
         
        right_asc = 2*atan((y*cos(b) + (x^2*cos(b)^2 + y^2*cos(b)^2 - z^2*sin(b)^2)^(1/2))/(x*cos(b) - z*sin(b)));
     else % perigee
        right_asc = 2*atan((y*cos(b) - (x^2*cos(b)^2 + y^2*cos(b)^2 - z^2*sin(b)^2)^(1/2))/(x*cos(b) - z*sin(b)));
     end
     


end