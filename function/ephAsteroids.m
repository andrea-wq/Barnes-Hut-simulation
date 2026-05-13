function [kep, mass, M] = ephAsteroids( time, id )%#codegen

% ephAsteroids.m - Ephemerides of Asteroids considernig the data retrieved
% from JPL data (See References) 
%
% PROTOTYPE:
%	[kep, mass, M] = ephAsteroids(time, id)
%
% DESCRIPTION:
%   This function returns the orbital parameters, the mass, and the mean
%   anomaly of Asteroids in the Solar System. The id corresponds to the num 
%   given by the JPL data (See References). The keplerian parameters at date
%   are retrieved through the Kepler solution of the parameters at epoch of
%   the data. 
%   Mass is obtained considering an interpolation of the H-d diagram.
%
% INPUT:
%   time [1]	Time of the ephemeris in MJD2000 [day]
%   id[1]       Asteroid identifier
%
% OUTPUT:
%   kep[1,6]    Keplerian parameters. It is a 6 entry vector:
%               [a e i Om om wom] where:
%                   a is the semimajor axis [km];
%                   e is the eccentricity;
%                   i is the inclination [rad];
%                   Om is the anomaly of the ascending node [rad];
%                   om is the anomaly of the pericentre [rad];
%                   wom is the true anomaly (from the pericentre) [rad].
%   mass        Mass of the Asteroid [kg]. It can be read from the database, or,
%               if not available, estimated by an approximate equation.
%   M           Mean anomaly at time [rad].
% 
% CALLED FUNCTIONS:
%   astroConstants
%
% REFERENCES:
%
%   JPL Data
%   - https://ssd.jpl.nasa.gov/sb/elem_tables.html#legend
%
% AUTHOR:
%   Giacomo Borelli,  November 2024, MATLAB, ephAsteroids.m
%
% CHANGELOG:
%   05/11/2024, Giacomo Borelli: First version
%   17/12/2024, Giacomo Borelli: bug fix in Newton iteration
%
% -------------------------------------------------------------------------

% Define some useful constants
AU = astroConstants(2);
mu_sun = astroConstants(4);

persistent asteroid_data;

if isempty(asteroid_data)

    if not(exist("AsteroidsElements_num.mat", 'file') == 2 ) 
        error("Please make sure that the file AsteroidsElements_num.mat with the JPL asteroid data is in your MATLAB path");
    end
    

    data_temp = load('AsteroidsElements_num.mat');
    
    % asteroid data is organised with epoch,sma,ecc,inc,RAAN,argPer,M,H (for mass
    % interpolation) in a vector. 
    % only for the value retrieved

    % data_temp = data_temp.data(id,:);
    data_temp = cell2mat(data_temp.data(:,[3,4,5,6,7,8,9,10]));

    asteroid_data = data_temp;
    clearvars data_temp;

end


% 
% % retrieve data from the asteroid data matrix of interest
sma = asteroid_data(id,2)*AU;              % [km] 
ecc = asteroid_data(id,3);                 % [-]
inc = asteroid_data(id,4)*pi/180;          % [rad]
w = asteroid_data(id,5)*pi/180;            % [rad]
node = asteroid_data(id,6)*pi/180;         % [rad]
M0 = asteroid_data(id,7)*pi/180;           % [rad]

pi2 = 2*pi;
orbital_period = pi2*sqrt((sma*sma*sma)/mu_sun); % seconds
mean_motion = pi2/orbital_period;                % rad/s


% Epoch of the elements represented as the Modified Julian Date (MJD)
epoch0_mjd = asteroid_data(id,1);                 % mjd [TBD] (TDB)	- days
epoch0_mjd2000 = epoch0_mjd - 51544.5;  % mjd2000 [TBD] (TBD) - days
epoch_mjd2000 = time;                   % mjd2000 [TBD] (TBD) - days

% Mean anomaly at time  ( epoch_mjd2000) provided as input
M = M0 + mean_motion*(epoch_mjd2000 - epoch0_mjd2000)*86400;  % rad

% reduce M to first revolution
n_rev = fix(M/pi2); % passes by peri-apsis
M = M - n_rev*pi2;

% solve the kepler problem to find true anomaly at epoch
%
% Newton iteration performed limiting number of iteration to 5. Low
% accuracy ephemerides.

% initial guess
Ek = M + ecc*sin(M)/(1-sin(M+ecc) + sin(M));

for i = 1:5
    F1 = ecc*cos(Ek) - 1;
    Ek = Ek + (Ek - ecc*sin(Ek) - M)/(F1);
end
% compute true anomaly from Eccentric anomaly
f = 2*atan2(sqrt(1+ecc)*tan(Ek*0.5),sqrt(1-ecc));
f = mod(f,2*pi);

% writes outputs keplerian elements and mean anomaly
kep = [sma,ecc,inc,node,w,f];
M = mod(M,2*pi);


% Absolute magnitude (asteroids only and set to "99.00" when unknown).
H = asteroid_data(id,8);           % [mag]

% To retrieve the mass, a polynomial fit of the H-d diagram
% (magnitude-diameter) is used. (not cery accurate)

d = -2.522e-2*H^5 + 3.2961*H^4 - 1.7249e2*H^3 + 4.5231e3*H^2 - 5.9509e4*H + 3.1479e5;

% estimate mass from the density of 2 kg/dm3
density = 2*1000; 
mass = 4*pi/3*(0.5*d)^3*density;











