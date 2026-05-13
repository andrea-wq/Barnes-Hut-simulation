function out= odefun(t,y,PerturbationFlag,index)
% Define the 2Bp system around index-planet subjected to perturbation (J2 and 3BP with sun) if
% selected by Perturbation flag. Defined in state vector elements [r,v]
%
%
% PROTOTYPE:
%   ode_function =@(t,y) odefun(t,y,1,index)  ---> perturbed problem analisys
%   ode_function =@(t,y) odefun(t,y,0,index)  ---> unperturbed problem analisys
%
%
% INPUT:
%   t[1]                    Time (needed for perturbation(t...) fcn)
%   y[1x6]                  state vector of s/c (rx ry rz vx vy vz)
%                               [Km , Km/s]
%   PerturbationFlag[1]     Flag to turn on or off perturbation in the problem (
%                               ( 0 = off ; 1 = on ) 
%   Index[1]                Index of the planet (mercury=1 , venus = 2 ....)
%
%
% OUTPUT:
%   out[6x1]                Derivative of the keplerian elements
%                               [Km/s , Km/s^2]
%
% VERSIONS:
%   2025-11-10        First Version
%   2025-11-18        Added Flag value for perturb-unperturbed switch
%   2025-12-20      substituted J2 perturbation with perturbation() fcn
%   2026-01-03:     Last version
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti

    r = y(1:3);
    v = y(4:6);
    
    mu = astroConstants(10+index);
    kep = car2kep(y,astroConstants(10+index));
    acc = perturbation(t,kep,index);

    out = [           v
           -mu/norm(r)^3*r + PerturbationFlag*(acc) ];
end