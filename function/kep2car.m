function [r, v] = kep2car(kep, mu)
% KEP2CAR  Convert Keplerian elements to Cartesian state vectors.
%
% PROTOTYPE:
%   [r, v] = kep2car(kep, mu)
%
% INPUTS:
%   kep : [1x6] or [6x1] vector of Keplerian elements:
%         kep = [a, e, i, Omega, omega, f]
%         a     - semi-major axis [km]
%         e     - eccentricity
%         i     - inclination [rad]
%         Omega - right ascension of ascending node [rad]
%         omega - argument of pericenter [rad]
%         f     - true anomaly [rad]
%
%   mu  : gravitational parameter [km^3/s^2]
%
% OUTPUTS:
%   r : position vector in ECI frame [km]
%   v : velocity vector in ECI frame [km/s]

% VERSIONS:
%
%  2026-01-03 Last version 
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti

    % ---------------------------------------------------------------------
    % Extract elements from vector
    % ---------------------------------------------------------------------
    a     = kep(1);
    e     = kep(2);
    i     = kep(3);
    Omega = kep(4);
    omega = kep(5);
    f     = kep(6);
    AU_km = 1.495978707e8;  % [km]

    % if abs(a) < 1.0e3
    %     % Assume a is in AU and convert to km
    %     a = a * AU_km;
    % end
    % ---------------------------------------------------------------------
    % Parameter p
    % ---------------------------------------------------------------------
    p = a * (1 - e^2);

    % ---------------------------------------------------------------------
    % Position and velocity in perifocal PQW frame
    % ---------------------------------------------------------------------
    r_pf = (p / (1 + e*cos(f))) * [cos(f); sin(f); 0];
    v_pf = sqrt(mu/p) * [-sin(f); e + cos(f); 0];

    % ---------------------------------------------------------------------
    % Rotation matrices (3-1-3 Euler sequence)
    % ---------------------------------------------------------------------
    R3_Omega = [ cos(Omega) -sin(Omega) 0;
                 sin(Omega)  cos(Omega) 0;
                 0           0          1];

    R1_i = [1 0 0;
            0 cos(i) -sin(i);
            0 sin(i)  cos(i)];

    R3_omega = [ cos(omega) -sin(omega) 0;
                 sin(omega)  cos(omega) 0;
                 0           0          1];

    % ---------------------------------------------------------------------
    % Overall transformation from PQW to ECI
    % ---------------------------------------------------------------------
    Q_pX = R3_Omega * R1_i * R3_omega;

    % ---------------------------------------------------------------------
    % Transform to ECI
    % ---------------------------------------------------------------------
    r = Q_pX * r_pf;
    v = Q_pX * v_pf;

end