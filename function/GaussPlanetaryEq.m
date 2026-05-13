function out = GaussPlanetaryEq(t,kep,PerturbationFlag,index)
% Define the 2Bp system around index-planet subjected to perturbation (J2 and 3BP with sun) if
% selected by Perturbation flag. Defined in keplerian elements
%
%
% PROTOTYPE:
%   ode_function =@(t,y) GaussPlanetaryEq(t,y,1,index)  ---> perturbed problem analisys
%   ode_function =@(t,y) GaussPlanetaryEq(t,y,0,index)  ---> unperturbed problem analisys
%
%
% INPUT:
%   t[1]                    Time (needed for perturbation(t...) fcn)
%   kep[1x6]                Keplerian element of s/c ( a , e , i , W , w , theta)
%                               [Km , , rad,rad,rad,rad]
%   PerturbationFlag[1]     Flag to turn on or off perturbation in the problem (
%                               ( 0 = off ; 1 = on ) 
%   Index[1]                Index of the planet (mercury=1 , venus = 2 ....)
%
%
% OUTPUT:
%   out[6x1]                Derivative of the keplerian elements
%                               [Km/s ,  , rad/s , rad/s , rad/s ,rad/s]
%
% VERSIONS:
%   2025-12-        First Version (already with perturbation() fcn call)
%   2025-12-20      Bug-fix
%   2025-12-27      minor Bug-fix
%
%  CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti

    a       = kep(1);
    e       = kep(2);
    i       = kep(3);
    W       = kep(4);  % useless
    omega   = kep(5);
    f       = kep(6);

    p = a*(1-e^2);                              % Semi latum rectum calc
    h = sqrt(astroConstants(10 + index)*p);     % Angular momentum calc
    acc = perturbation(t,kep,index);            % Perturbation function call (time instant, keplerian, index of planet)
    
    [r_sat,v_sat] = kep2car(kep,astroConstants(10 + index));
    r = norm(r_sat);

% frame definition for gauss analisys

    r_cap = r_sat/r;
    w_cap = cross(r_sat,v_sat)/norm(cross(r_sat,v_sat));
    s_cap = cross(w_cap,r_cap);

% change of frame for disturbance acceleration

    a_w = dot(acc,w_cap) * PerturbationFlag;
    a_s = dot(acc,s_cap) * PerturbationFlag;
    a_r = dot(acc,r_cap) * PerturbationFlag;

% Differential Gauss problem definition

    out = zeros(6,1);

    out(1) = (2*a^2)/h * (e*sin(f)*a_r + p/r*a_s);              % da/dt

    out(2) = 1/h * (p*sin(f)*a_r + ((p+r)*cos(f) + r*e)*a_s);   % de/dt

    out(3) = (r*cos(f+omega))/h * a_w;                          % di/dt

    out(4) = (r*sin(f+omega))/(h*sin(i)) * a_w;                 % dW/dt

    out(5) = 1/(h*e) * (-p*cos(f)*a_r + (p+r)*sin(f)*a_s) - (r*sin(f+omega)*cos(i))/(h*sin(i)) * a_w;   % dw/dt

    out(6) = h/(r^2) + 1/(h*e) * (p*cos(f)*a_r - (p+r)*sin(f)*a_s);   % df/dt


end