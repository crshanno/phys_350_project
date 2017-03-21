classdef Params
    properties (Constant)
        G = 6.67408E-11; % Gravitational Constant (m^3 kg^-1 s^-2)
        Me = 5.9722E24; % Mass of Earth (kg)
        Re = physconst('earthradius'); % Radius of Earth (meters)

        % for Low Earth Orbit
        leoOrbitMin = Params.Re + 160E3; % Minimum orbit radius (meters)
        leoOrbitMax = Params.Re + 2E6; % Maximum orbit radius (meters)
    end
end