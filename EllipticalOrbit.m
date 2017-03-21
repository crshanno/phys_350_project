% https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
% higher computation time than circular orbit, since E is solved numerically
classdef EllipticalOrbit < Orbit
    properties (Access = public)
        P % Period of orbit (seconds)
        epsilon % Eccentricity
        a % Semimajor axis (meters)
        theta_0 % (radians)
        incX % inclination x-axis (radians)
        incY % inclination y-axis (radians)
    end
    
    properties (Access = private)
        % cache the last computed E value so it doesn't have to be
        % re-computed for the same t sequentially
        last_t;
        last_E;
    end
    
    methods (Access = public)
        function obj = EllipticalOrbit(mass, a, epsilon, theta_0, incX, incY)
            % constructor
            % mass: mass of orbiting object (kg)
            % epsilon: eccenticity (0 < e < 1)
            % a: semimajor axis (meters)
            % theta_0: initial true anomaly (radians)
            % inclination: (radians)
            
            if epsilon <= 0 || epsilon >= 1
                error('invalid eccentricity: %f', epsilon);
            end
            obj.a = a;
            obj.epsilon = epsilon;
            obj.theta_0 = theta_0;
            obj.P = Orbit.calcP(mass, obj.a);
            obj.incX = incX;
            obj.incY = incY;
        end
        
        function theta = theta(obj, t)
            % True Anomaly
            % t: time since start (seconds)
            % returns theta (radians)
            E = obj.E(t);
            theta = sqrt( (1 + obj.epsilon) / (1 - obj.epsilon) ) * tan (E / 2);
            theta = 2 * atan(theta);
            theta = theta + obj.theta_0;
        end
        
        function r = r(obj, t)
            % current distance from center of orbit
            % t: time since start (seconds)
            % returns r (meters)
            E = obj.E(t);
            r = obj.a*(1 - obj.epsilon*cos(E));
        end
        
        function [x,y,z] = getXYZ(obj, t)
            % get position vector in cartesian coordinates
            % t: time since start (seconds)
            % returns position vector (meters)
            r = obj.r(t);
            theta = obj.theta(t);
            
            r = r*[cos(theta) sin(theta) 0];
            R = rotx(rad2deg(obj.incX))*roty(rad2deg(obj.incY));
            r = r*R;
            
            x = r(1);
            y = r(2);
            z = r(3); 
        end
        
    end
    methods (Access = private)        
        function M = M(obj, t)
            % Mean Anomaly
            % t: time since start (seconds)
            % returns mean anomaly M
            n = 2*pi/obj.P;
            M = n*t;
        end
        
        function E = E(obj, t)
            % Eccentric Anomaly
            % t: time since start (seconds)
            % returns eccentric anomaly E
            
            % note: E is solved numerically as there is no simple solution
            % to M = E - \epsilon*sin(E)
            % and hence this computation is slow.
            % perhaps it can be optimized and/or approximated
            
            m = obj.M(t);
            function F = ef(x)
               F = x - obj.epsilon*sin(x) - m;
            end
            if t == obj.last_t
                E = obj.last_E;
            else
                eqn = @ef;
                E = fzero(eqn, 0);
                obj.last_t = t;
                obj.last_E = E;
            end
        end
    end        
end