% Effectively an elliptical orbit with 0 eccentricity.
% Computation is simpler than elliptical orbit
classdef CircularOrbit < Orbit
    properties (Access = public)
        P % Period of orbit (seconds)
        a % radius (meters)
        theta_0 % (radians)
        incX % inclination x-axis (radians)
        incY % inclination y-axis (radians)
    end
    
    methods (Access = public)
        function obj = CircularOrbit(mass, r, theta_0, incX, incY)
            % constructor 
            % mass: mass of orbiting object (kg)
            % r: radius of orbit (meters)
            % theta_0: initial true anomaly (radians)
            % incX: inclination about x axis(radians) 
            % incY: inclination about y axis(radians)             
            obj.P = Orbit.calcP(mass, r);
            obj.a = r;
            obj.incX = incX;
            obj.incY = incY;
            obj.theta_0 = theta_0;
        end
        
        function theta = theta(obj, t)
            % True Anomaly
            % t: time since start (seconds)
            % returns theta (radians)
            n = 2*pi / obj.P;
            theta = n*t + obj.theta_0;
        end
        
        function r = r(obj, t)
            % current distance from center of orbit
            % although it i s constant for circular orbit, it is included as
            % a function for interoperability
            % t: time since start (seconds)
            % returns r (meters)
           r = obj.a;
        end
        
        function [x,y,z] = getXYZ(obj, t)
            % get position vector in cartesian coordinates
            % t: time since start (seconds)
            % returns position vector (meters)
            r = obj.r();
            theta = obj.theta(t);
            
            r = r*[cos(theta) sin(theta) 0];
            R = rotx(rad2deg(obj.incX))*roty(rad2deg(obj.incY));
            r = r*R;
            
            x = r(1);
            y = r(2);
            z = r(3); 
        end
    end
end