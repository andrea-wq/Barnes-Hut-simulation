function kep = car2kep(s, mu)
%CAR2KEP Converts Cartesian state to Keplerian elements
%
% PROTOTYPE:
% kep = car2kep(s, mu)
%
% DESCRIPTION:
% Converts the Cartesian position and velocity vectors of a body into the
% corresponding Keplerian orbital elements, with basic handling of
% singular cases (circular and/or equatorial orbits).
%
% INPUT:
% s   State vector [6x1] = [r; v]
%     r   Position vector in inertial frame [km]
%     v   Velocity vector in inertial frame [km/s]
% mu  Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT:
% kep  Keplerian elements vector:
%      [a, e, i, Omega, omega, theta]
%      where:
%          a      semi-major axis [km]
%          e      eccentricity [-]
%          i      inclination [rad]
%          Omega  RAAN [rad]
%          omega  argument of periapsis [rad]
%          theta  true anomaly [rad]
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Bavetta
% Sophia Trivellin
% Stefano Zotii
%
% VERSIONS:
% 2026-01-03: Last version

r = s(1:3);                 % Position vector [km]
v = s(4:6);                 % Velocity vector [km/s]

% ===== Numerical tolerances for singularities =====
tol_i = 1e-10;              % Inclination tolerance
tol_e = 1e-10;              % Eccentricity tolerance

% Ensure column vectors
r = r(:);
v = v(:);

r_norm = norm(r);           % Position magnitude
v_norm = norm(v);           % Velocity magnitude

v_r = dot(r,v)/r_norm;      % Radial velocity component

% Specific angular momentum
h = cross(r,v);             % Angular momentum vector
h_norm = norm(h);           % Angular momentum magnitude

% Semi-major axis (from specific energy)
energy = v_norm^2/2 - mu/r_norm;
if abs(energy) > 1e-14
    a = -mu/(2*energy);
else
    a = Inf;                % Parabolic case (energy  0)
end

% Inclination 
i = acos(h(3)/h_norm);

% Eccentricity vector 
e_vec = (1/mu)*(cross(v,h) - mu*r/r_norm);
e = norm(e_vec);            % Eccentricity magnitude

if e > tol_e
    e_hat = e_vec/e;        % Unit vector 
else
    e = 0;                  % Circular orbit
    e_hat = [1;0;0];        % Placeholder direction
end

% Node vector
K = [0;0;1];                % Z axis
N = cross(K,h);             % Node vector
N_norm = norm(N);

% Right Ascension of Ascending Node (RAAN) 
if i > tol_i && N_norm > 0
    Omega = atan2(N(2), N(1));
else
    Omega = 0;              % Equatorial orbit: RAAN undefined
end
Omega = real(Omega);
Omega = mod(Omega, 2*pi);

% Argument of periapsis 
if e > tol_e && i > tol_i
    omega = acos(dot(N,e_hat)/N_norm);
    if e_hat(3) < 0
        omega = 2*pi - omega;
    end
elseif e > tol_e && i <= tol_i
    % Equatorial, non-circular orbit
    omega = atan2(e_vec(2), e_vec(1));
else
    % Circular orbit: periapsis undefined
    omega = 0;
end
omega = real(omega);
omega = mod(omega, 2*pi);

% True anomaly 
if e > tol_e
    theta = acos(dot(e_hat,r)/r_norm);
    if v_r < 0
        theta = 2*pi - theta;
    end
else
    % Circular orbit
    theta = atan2(r(2), r(1));
end
theta = real(theta);
theta = mod(theta, 2*pi);

% Output 
kep = [a, e, i, Omega, omega, theta];

end


