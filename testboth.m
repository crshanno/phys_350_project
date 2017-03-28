function testboth(N)
    if (nargin < 1)
        N = 5;
    end
    figure;
    subplot(1,2,1);
    orbittest(0,N);
    title('circular orbits');    
    subplot(1,2,2);
    orbittest(1,N);
    title('eliptical orbits');
end