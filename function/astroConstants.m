function out = astroConstants(in)

% astroConstants.m - Returns astrodynamic-related physical constants.
%
% PROTOTYPE:
%   out = astro_constants(in)
%
% DESCRIPTION:
%   Returns a row vector of constants, in which there is the corresponding
%   constant for each element of the input vector.
%
%   List of identifiers:
%       Generic astronomical constants:
%           1   Universal gravity constant (G) (from DITAN and Horizon) [km^3/(kg*s^2)]
%           2   Astronomical Unit (AU) (from DE405) [km]
%               Note:  The value for 1 au is from the IAU 2012 Resolution B1.
%       Sun related:
%           3   Sun mean radius (from DITAN) [km]
%           4   Sun planetary constant (mu = mass * G) (from DE405)
%               [km^3/s^2]
%           81  Energy flux density of the Sun (from Wertz,SMAD)
%               [W/m2 at 1 AU]
%       Other:
%           5   Speed of light in the vacuum (definition in the SI and Horizon) [km/s]
%           6   Standard free fall (the acceleration due to gravity on the
%               Earth's surface at sea level) (from Wertz,SMAD) [m/s^2]
%           7   Mean distance Earth-Moon (from Wertz,SMAD) [km]
%           8   Obliquity (angle) of the ecliptic at Epoch 2000 (from
%               Horizon) [rad]
%           9   Gravitatonal field constant of the Earth (from Wertz,SMAD,
%               taken from JGM-2). This should be used in conjunction to
%               Earth radius = 6378.1363 km
%           82  Earth days in a Julian year y = 365.25 d  (from Horizon)
%       Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
%           11  Me      (from DE405)
%           12  V       (from DE405)
%           13  E       (from DE405)
%           14  Ma      (from DE405)
%           15  J       (from DE405)
%           16  S       (from DE405)
%           17  U       (from DE405)
%           18  N       (from DE405)
%           19  P       (from DE405)
%           20  Moon    (from DE405)
%       Mean radius of the planets [km]:
%           21  Me      (from Horizon)
%           22  V       (from Horizon)
%           23  E       (from Horizon)
%           24  Ma      (from Horizon)
%           25  J       (from Horizon)
%           26  S       (from Horizon)
%           27  U       (from Horizon)
%           28  N       (from Horizon)
%           29  P       (from Horizon)
%           30  Moon    (from Horizon)
%       J2 Gravitational Harmonic coefficent [-]:
%           31  Me      (from Nasa Fact Sheet)
%           32  V       (from Nasa Fact Sheet)
%           33  E       (from Nasa Fact Sheet)
%           34  Ma      (from Nasa Fact Sheet)
%           35  J       (from Nasa Fact Sheet)
%           36  S       (from Nasa Fact Sheet)
%           37  U       (from Nasa Fact Sheet)
%           38  N       (from Nasa Fact Sheet)
%           39  P       (from Nasa Fact Sheet)
%           40  Moon    (from Nasa Fact Sheet)
%       Planetary oblateness [-]:
%           41  Me      (from Nasa Fact Sheet)
%           42  V       (from Nasa Fact Sheet)
%           43  E       (from Nasa Fact Sheet)
%           44  Ma      (from Nasa Fact Sheet)
%           45  J       (from Nasa Fact Sheet)
%           46  S       (from Nasa Fact Sheet)
%           47  U       (from Nasa Fact Sheet)
%           48  N       (from Nasa Fact Sheet)
%           49  P       (from Nasa Fact Sheet)
%           50  Moon    (from Nasa Fact Sheet)
%       Sidereal rotation period [hours]:
%           51  Me      (from Nasa Fact Sheet)
%           52  V       (from Nasa Fact Sheet)
%           53  E       (from Nasa Fact Sheet)
%           54  Ma      (from Nasa Fact Sheet)
%           55  J       (from Nasa Fact Sheet)
%           56  S       (from Nasa Fact Sheet)
%           57  U       (from Nasa Fact Sheet)
%           58  N       (from Nasa Fact Sheet)
%           59  P       (from Nasa Fact Sheet)
%           60  Moon    (from Nasa Fact Sheet)
%       Axial tilt (inclination of the equator with respect to the
%       ecliptic) [deg]
%           61  Me      (from Nasa Fact Sheet)
%           62  V       (from Nasa Fact Sheet)
%           63  E       (from Nasa Fact Sheet)
%           64  Ma      (from Nasa Fact Sheet)
%           65  J       (from Nasa Fact Sheet)
%           66  S       (from Nasa Fact Sheet)
%           67  U       (from Nasa Fact Sheet)
%           68  N       (from Nasa Fact Sheet)
%           69  P       (from Nasa Fact Sheet)
%           70  Moon    (from Nasa Fact Sheet)
%       Solar irradiance [W/m^2]
%           71  Me      (from Nasa Fact Sheet)
%           72  V       (from Nasa Fact Sheet)
%           73  E       (from Nasa Fact Sheet)
%           74  Ma      (from Nasa Fact Sheet)
%           75  J       (from Nasa Fact Sheet)
%           76  S       (from Nasa Fact Sheet)
%           77  U       (from Nasa Fact Sheet)
%           78  N       (from Nasa Fact Sheet)
%           79  P       (from Nasa Fact Sheet)
%           80  Moon    (from Nasa Fact Sheet)
%
%   Notes for upgrading this function:
%       It is possible to add new constants.
%       - DO NOT change the structure of the function, as well as its
%           prototype.
%       - DO NOT change the identifiers of the constants that have already
%           been defined in this function. If you want to add a new
%           constant, use an unused identifier.
%       - DO NOT add constants that can be easily computed starting form
%           other ones (avoid redundancy).
%       Contact the author for modifications.
%
% INPUT:
%   in      Vector of identifiers of required constants.
%
% OUTPUT:
%   out     Vector of constants.
%
% EXAMPLE:
%   astroConstants([2, 4, 26])
%      Returns a row vector in which there is the value of the AU, the Sun
%      planetary constant and the mean radius of Saturn.
%
%   astroConstants(10 + [1:9])
%      Returns a row vector with the planetary constant of each planet.
%
% REFERENCES:
%   - DITAN (Direct Interplanetary Trajectory Analysis), Massimiliano
%       Vasile, 2006.
%	- Wertz J. R., Larson W. J., "Space Mission Analysis and Design", Third
%       Edition, Space Technology Library 2003.
%   [A]   DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
%   [B]   Explanatory Supplement to the Astronomical Almanac. 1992. K. P.
%         Seidelmann, Ed., p.706 (Table 15.8) and p.316 (Table 5.8.1),
%         University Science Books, Mill Valley, California. 
%   [C]   Tholen, D.J. and Buie, M.W. 1990. "Further Analysis of
%         Pluto-Charon Mutual Event Observations" BAAS 22(3):1129.
%   [D]   Seidelmann, P.K. et al. 2007. "Report of the IAU/IAG Working
%         Group on cartographic coordinates and rotational elements: 2006"
%         Celestial Mech. Dyn. Astr. 98:155-180. 
%   [F]   Anderson, J.D., et al. 1987. "The mass, gravity field, and
%         ephemeris of Mercury" Icarus 71:337-349.
%   [G]   Konopliv, A.S., et al. 1999. "Venus gravity: 180th degree and
%         order model" Icarus 139:3-18.
%   [H]   Folkner, W.M. and Williams, J.G. 2008. "Mass parameters and
%         uncertainties in planetary ephemeris DE421." Interoffice Memo.
%         343R-08-004 (internal document), Jet Propulsion Laboratory,
%         Pasadena, CA. 
%   [I]   Jacobson, R.A. 2008. "Ephemerides of the Martian Satellites -
%         MAR080" Interoffice Memo. 343R-08-006 (internal document),
%         Jet Propulsion Laboratory, Pasadena, CA. 
%   [J]   Jacobson, R.A. 2005. "Jovian Satellite ephemeris - JUP230"
%         private communication. 
%   [K]   Jacobson, R.A., et al. 2006. "The gravity field of the Saturnian
%         system from satellite observations and spacecraft tracking data"
%         AJ 132(6):2520-2526. 
%   [L]   Jacobson, R.A. 2007. "The gravity field of the Uranian system and
%         the orbits of the Uranian satellites and rings" BAAS 39(3):453. 
%   [M]   Jacobson, R.A. 2008. "The orbits of the Neptunian satellites and
%         the orientation of the pole of Neptune" BAAS 40(2):296. 
%   [N]   Jacobson, R.A. 2007. "The orbits of the satellites of Pluto -
%         Ephemeris PLU017" private communication.
%   [W1]  http://ssd.jpl.nasa.gov/?planet_phys_par Last retrieved
%         20/03/2013
%   [W2]  http://ssd.jpl.nasa.gov/?sat_phys_par Last retrieved
%         20/03/2013
%   [W3]  http://ssd.jpl.nasa.gov/horizons.cgi Last retrieved
%         20/03/2013
%   [M1]  Bills, B.G. and Ferrari, A.J. 1977. ``A Harmonic Analysis of
%         Lunar Topography'', Icarus 31, 244-259.
%   [M2]  Standish, E. M. 1998. JPL Planetary and Lunar Ephemerides,
%         DE405/LE405.
%   [M3]  Lunar Constants and Models Document, Ralph B. Roncoli, 23 Sept 2005,
%         JPL Technical Document D-32296 
%   [P1]  Mercury Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
%   [P2]  Venus Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
%   [P3]  Earth Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
%   [P4]  Mars Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
%   [P5]  Jupiter Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
%   [P6]  Saturn Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
%   [P7]  Uranus Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html
%   [P8]  Neptune Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
%   [P9]  Pluto Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
%   [P10] Moon Fact sheet: last retrieved 21/11/2024
%         https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
%
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 2006, MATLAB, astroConstants.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 2006, MATLAB, astro_constants.m, Ver. 1.2
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   26/10/2006, Camilla Colombo: Updated.
%   22/10/2007, Camilla Colombo: astroConstants(8) added (Obliquity (angle)
%       of the ecliptic at Epoch 2000).
%   02/10/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   12/11/2010, Camilla Colombo: astroConstants(9) added (J2) Note: the
%       present value of J2 is not consistent with the value of the Earth
%       radius. This value of J2 should be used in conjunction to Earth
%       radius = 6378.1363 km
%   19/03/2013, Camilla Colombo: constants updated to NASA JPL website.
%       References added.
%   20/03/2013, REVISION, Francesca Letizia.
%   22/03/2013, Francesca Letizia: all GM from DE405.
%   21/11/2024, Mathilda Bolis: added constants related to other planets.
%
% -------------------------------------------------------------------------

% 9: J2
% 82: 365.25

out = zeros(1,length(in));
for i=1:length(in)
    switch in(i)
        case 1
            out(i)=6.67259e-20; % From DITAN and Horizon
        case 2
            out(i)=149597870.691; % From DE405
        case 3
            % out(i)=700000; % From DITAN
            out(i)=6.955*10^5; % From Horizon [W3]
        case 4
            % out(i)=0.19891000000000E+31*6.67259e-20; % From DITAN
            out(i)=1.32712440017987E+11; % From DE405 [A]
        case 5
            out(i)=299792.458; % Definition in the SI, Horizon, DE405
        case 6
            out(i)=9.80665; % Definition in Wertz, SMAD
        case 7
            % out(i)=384401; % Definition in Wertz, SMAD
            out(i)=384400; % From Horizon [W3]
        case 8
            % out(i)=23.43928111*pi/180; % Definition in Wertz, SMAD
            out(i)=84381.412/3600*pi/180; % Definition in Horizon
            % obliquity of ecliptic (J2000)    epsilon = 84381.412 (± 0.005) arcsec 
        case 9
            out(i)=0.1082626925638815e-2; % Definition in Wertz, SMAD
        case 11
            % out(i)=0.33020000000000E+24*6.67259e-20; % From DITAN
            %out(i)=0.330104E+24*6.67259e-20;    % From Horizon [F]
            out(i)=2.203208E+4;    % From DE405
        case 12
            % out(i)=0.48685000000000E+25*6.67259e-20; % From DITAN
            %out(i)=4.86732E+24*6.67259e-20;     % From Horizon [G]
            out(i)=3.24858599E+5; % From DE405
        case 13
            % out(i)=0.59736990612667E+25*6.67259e-20; % From DITAN
            % out(i)=5.97219E+24*6.67259e-20;     % From Horizon [H]
            out(i) = 3.98600433e+5; % From DE405
        case 14
            % out(i)=0.64184999247389E+24*6.67259e-20; % From DITAN
            %out(i)=0.641693E+24*6.67259e-20; 	% From Horizon [I]
            out(i) = 4.2828314E+4; %Frome DE405
        case 15
            % out(i)=0.18986000000000E+28*6.67259e-20; % From DITAN
            %out(i)=1898.13E+24*6.67259e-20; 	% From Horizon [J]
            out(i) = 1.26712767863E+08; % From DE405
        case 16
            % out(i)=0.56846000000000E+27*6.67259e-20; % From DITAN
            % out(i)=568.319E+24*6.67259e-20;     % From Horizon [k]
            out(i) = 3.79406260630E+07; % From DE405
        case 17
            % out(i)=0.86832000000000E+26*6.67259e-20; % From DITAN
            % out(i)=86.8103E+24*6.67259e-20;     % From Horizon [L]
            out(i)= 5.79454900700E+06; % From DE405
        case 18
            % out(i)=0.10243000000000E+27*6.67259e-20; % From DITAN
            % out(i)=102.410E+24*6.67259e-20;     % From Horizon [M]
            out(i) = 6.83653406400E+06; % From DE405
        case 19
            % out(i)=0.14120000000000E+23*6.67259e-20; % From DITAN
            %out(i)=.01309E+24*6.67259e-20;     % From Horizon [N]
            out(i) = 9.81601000000E+02; %From DE405
        case 20
            % out(i)=0.73476418263373E+23*6.67259e-20; % From DITAN
             out(i)=4902.801;                 % From Horizon  [M2]
            %out(i)=4902.801076;                % From Horizon  [M3]
        case 21
            % out(i)=0.24400000000000E+04; % From DITAN
            out(i)=2439.7; % From Horizon [D]
        case 22
            % out(i)=0.60518000000000E+04; % From DITAN
            out(i)=6051.8; % From Horizon [D]
        case 23
            % out(i)=0.63781600000000E+04; % From DITAN
            % out(i)=6371.00; % From Horizon [B]
            out(i)=6371.01; % From Horizon [W3]
        case 24
            % out(i)=0.33899200000000E+04; % From DITAN
            % out(i)=3389.50; % From Horizon [D]
            out(i)=3389.9; % From Horizon [W3]            
        case 25
            % out(i)=0.69911000000000E+05; % From DITAN
            out(i)=69911;   % From Horizon [D]
        case 26
            % out(i)=0.58232000000000E+05; % From DITAN
            out(i)=58232;   % From Horizon [D]
        case 27
            % out(i)=0.25362000000000E+05; % From DITAN
            out(i)=25362;   % From Horizon [D]
        case 28
            % out(i)=0.24624000000000E+05; % From DITAN
            % out(i)=24622;   % From Horizon [D]
            out(i)= 24624; % From Horizon [W3]            
        case 29
            % out(i)=0.11510000000000E+04; % From DITAN
            out(i)=1151; 	% From Horizon [C]
        case 30
            % out(i)=0.17380000000000E+04; % From DITAN
            % out(i)=1737.5;  % From Horizon [M1]
            out(i)=1738.0;    % From Horizon  [M3]
        % J2 Gravitational Harmonic coefficent [-]:
        case 31
            out(i) = 50.3*1e-6; % [P1]
        case 32
            out(i) = 4.458*1e-6; % [P2]
        case 33
            out(i) = 0.1082626925638815e-2; % Definition in Wertz, SMAD
        case 34
            out(i) = 1960.45*1e-6; % [P4]
        case 35
            out(i) = 14696.5735*1e-6; % [P5]
        case 36
            out(i) = 16290.573*1e-6; % [P6]
        case 37
            out(i) = 3510.68*1e-6; % [P7]
        case 38
            out(i) = 3408.43*1e-6; % [P8]
        case 39
            out(i) = 9.01*10^-4; % [P9] -> Not present in Nasa fact sheets
        case 40
            out(i) = 202.7*1e-6; % [P10]
        % Planetary oblateness
        case 41
            out(i) = 0.0009; % [P1]
        case 42
            out(i) = 0.00001; % [P2]
        case 43
            out(i) = 0.00335; % [P3]
        case 44
            out(i) = 0.00589; % [P4]
        case 45
            out(i) = 0.06487; % [P5]
        case 46
            out(i) = 0.09796; % [P6]
        case 47
            out(i) = 0.02293; % [P7]
        case 48
            out(i) = 0.01708; % [P8]
        case 49
            out(i) = 0; % [P9]
        case 50
            out(i) = 0.0012; % [P10]
        % Sidereal rotation period [hours]:
        case 51
            out(i) = 1407.6; % [P1]
        case 52
            out(i) = -5832.6; % [P2]
        case 53
            out(i) = 23.9345; % [P3]
        case 54
            out(i) = 24.6229; % [P4]
        case 55
            out(i) = 9.9250; % [P5]
        case 56
            out(i) = 10.656; % [P6]
        case 57
            out(i) = -17.24; % [P7]
        case 58
            out(i) = 16.11; % [P8]
        case 59
            out(i) = -153.2928; % [P9]
        case 60
            out(i) = 655.720; % [P10]
        % Axial tilt [deg]:
        case 61
            out(i) = 0.034; % [P1]
        case 62
            out(i) = 177.36; % [P2]
        case 63
            out(i) = 23.44; % [P3]
        case 64
            out(i) = 25.19; % [P4]
        case 65
            out(i) = 3.13; % [P5]
        case 66
            out(i) = 26.73; % [P6]
        case 67
            out(i) = 97.77; % [P7]
        case 68
            out(i) = 28.32; % [P8]
        case 69
            out(i) = 119.51; % [P9]
        case 70
            out(i) = 6.68; % [P10]
        % Solar irradiance [W/m^2];
        case 71
            out(i) = 9082.7; % [P1]
        case 72
            out(i) = 2601.3; % [P2]
        case 73
            out(i) = 1361.0; % [P3]
        case 74
            out(i) = 586.2; % [P4]
        case 75
            out(i) = 50.26; % [P5]
        case 76
            out(i) = 14.82; % [P6]
        case 77
            out(i) = 3.69; % [P7]
        case 78
            out(i) = 1.508; % [P8]
        case 79
            out(i) = 0.873; % [P9]
        case 80
            out(i) = 1361.0; % [P10]
        % Other constants
        case 81
            out(i)=1367; % From Wertz, SMAD
            % out(i)=1367.6;  % From Horizon  [W3]
        case 82
            out(i)=365.25; % From Horizon
        % Add an identifier and constant here. Prototype:
        % case $identifier$
        %     out(i)=$constant_value$;
        otherwise
            warning('Constant identifier %d is not defined!',in(i));
            out(i)=0;
    end
end
