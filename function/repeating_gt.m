function a = repeating_gt(k,m,we,index)
%repeating_gt Semi-major axis for a repeating ground-track orbit
%
% PROTOTYPE:
% a = repeating_gt(k, m, we, index)
%
% DESCRIPTION:
% This function computes the semi-major axis of a repeating ground-track
% (RGT) orbit based on the resonance condition between the satellite mean
% motion and the planet's rotation rate. A k:m resonance is enforced such
% that the satellite completes k revolutions while the central body rotates
% m times beneath the orbital plane. The resulting mean motion is then used
% to derive the corresponding semi-major axis from Kepler's third law.
%
% INPUT:
% k[1]      Number of satellite orbital revolutions per ground-track cycle [-]
% m[1]      Number of central-body rotations per ground-track cycle [-]
% we[1]     Central-body rotation rate expressed in arcseconds per second
%           (converted internally to rad/s) [arcsec/s]
% index[1]  Index used to retrieve the gravitational parameter of the
%           central body via astroConstants [–]
%
% OUTPUT:
% a[1] Semi-major axis of the repeating ground-track orbit [km]
%
% Version
% 2025-12-22 First Version
% 2026-01-03: Last version
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti
%

    w_e = deg2rad(we/3600);
    mu = astroConstants(10+index);
    
    n = w_e * k/m;
    a = (mu/n^2)^(1/3);
end