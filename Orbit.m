% https://en.wikipedia.org/wiki/Geocentric_orbit
% https://earthobservatory.nasa.gov/Features/OrbitsCatalog/

% for kessler syndrome, consider Low Earth Orbit (LEO)
% https://en.wikipedia.org/wiki/Low_Earth_orbit

% 200–2,000 km above mean sea level

% circular orbit: 7.8–6.9 km/s (17,450–14,430 mph) respectively
% elliptic orbit: 6.5–8.2 km/s respectively

% orbital period: 1 h 29 min – 2 h 8 min

classdef (Abstract) Orbit
    properties
        % when the object is a sphere, store its radius for future
        % collision checking
        sphereRadius;  % meters
    end
    
    methods (Abstract)
        theta(obj, t)
        r(obj, t)
        getXYZ(obj, t)
    end 
    
    methods (Static)
        function P = calcP(mass, a)
            % calculate period of orbit given mass and semi major axis
            % (radius for circular orbit)
            % mass: mass of orbiting object (kg)
            % a: semimajor axis of orbit (ie. orbit radius) (meters)
            % return orbit period (seconds)
            massSum = mass + Params.Me;
            P = sqrt( (4 * pi.^2 * (a).^3) / ( Params.G * massSum ));
        end
        function mass = calcMassSphere(density, r)
            % calculate mass of sphere with given constant density and radius
            % density: (kg/m^3)
            % r: radius of sphere (meters)
            % returns mass (kg)
            mass = density * (4/3)*pi*r.^3;
        end
    end
end
