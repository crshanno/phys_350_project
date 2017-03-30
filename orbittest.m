function orbittest(ellipticaltest, N)
    tic;
    
    if nargin < 1
        ellipticaltest = 0;
    end    
    
    % number of objects to simulate
    if (nargin < 2)
        N = 5;
    end
    
    % number of orbits to compute
    nOrbits = 10;
    
    % time scale resolution (seconds)
    % ie. how many seconds in between iterations
    time_resolution = 10;
    

    function r = myRandn(min, max, dev)
        % normalized random between bounds
        if (nargin < 3)
           dev = (min-max)/10; % make up something 
        end
        avg = (max+min)/2;
        r = avg + dev*randn;
        if r < min
            r = min;
        elseif r > max
            r = max;
        end
    end

    function r = myRand(min, max)
        % uniform random between bounds
        r = rand * (max - min) + min;
    end

    % create 100 orbits
    for i = 1:N
        % random radius and density
        objectRadius = myRandn(1e-2, 10);
        density = myRandn(100, 10000);
        
        % compute mass
        mass = Orbit.calcMassSphere(density, objectRadius);
        
        if ellipticaltest ~= 1
            orbits(i) = CircularOrbit(mass, ...
                ... % random orbit radius within leo
                myRandn(Params.leoOrbitMin, Params.leoOrbitMax), ...
                myRand(0, 2*pi), ...  % random initial theta
                myRand(0, pi), myRand(0, pi));     % random inclination
        else
            orbits(i) = EllipticalOrbit(mass, ...
                ... % random orbit radius within leo
                myRandn(Params.leoOrbitMin, Params.leoOrbitMax), ...
                myRandn(0.01, 0.99), ... % random eccentricity
                myRand(0, 2*pi), ...  % random initial theta
                myRand(0, pi), myRand(0, pi));     % random inclination        
        end
        orbits(i).sphereRadius = objectRadius;
    end
        
    % iterate over orbits
    for j = 1:N
        % iterate over time
        t = 0:time_resolution:(orbits(j).P*nOrbits);
        for i = 1:length(t)
             [x, y, z] = orbits(j).getXYZ(t(i));
             r(j,i,1) = x;
             r(j,i,2) = y;
             r(j,i,3) = z;
        end
    end
   
    
    % plot an orbit    
    function msub(i)    
        l = length(r(i,:,1));
        % uniform scale
        s = ones(1, l)*10;
        % vary color intensity with time
        color = 1:l;
        color = color';
        color = color/l;
        color = [color zeros(l, 2)];
        scatter3(r(i,:,1), r(i,:,2), r(i,:,3), s, color);
        xlabel('x');
        ylabel('y');
        zlabel('z');
        d = 2E7;
        d = [-d d];
        axis([d d d]);
        grid on;
    end

    hold on;
    for i = 1:N
        msub(i);
    end
    % plot the origin
    scatter3(0, 0, 0, 100, 'blue', '*');
    hold off;
    
    fprintf('Run %d iterations on %d object(s) (total: %d) in %.1f seconds\n', ...
        length(t), N, length(t)*N, toc);
    
end