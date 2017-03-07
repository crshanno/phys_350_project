m_earth = 5.972e24; %kg
r_earth = 6371000; %m
G = 6.67408e-11; %m^3/(kgs^2)
o = [1,200000,0, 0, 2*pi/(800)]; % mass(kg), orbital height(m), angle(deg), v_r(m/s), v_theta(deg/s)

dt = 0.1;

ag = @(r) G*m_earth/(r+r_earth)^2; %m/s^2
polar(r_earth);
p = polar(o(3),o(2)+r_earth, 'ok');
pause(dt);
for n=1:100000
    o(3) = o(3)+ o(5)*dt;
    o(2) = o(2)+ o(4)*dt;
    ag(o(2)) - o(2)*o(5)^2
    o(4) = o(4) + (ag(o(2)) - o(2)*o(5)^2)*dt
    o(5) = o(5) + (2*o(4)*o(5))*dt;
    if mod(n,10) == 0
        p = polar(o(3),o(2)+r_earth, 'ok');
        refresh;
        pause(0.01);
    end
end


       
       