function [long,lat] = GroundTrack (y,t,w_e,theta_0)
% Compute Longitude and Latitude from the solution of the propagation
% problem expressed in [r,v] form
%
%
% PROTOTYPE:
%   [long,lat] = GroundTrack (y,t,w,theta_0)
%
% INPUT:
%
%   y[Nx6]                  Position and velocity vectors from ode
%                               propagation (rx ry rz vx vy vz) [Km,Km/s]
%   t[N]                    Time (output vector from ode solver)
%   w_e                     angular velocity of selected planet [deg/s]
%   theta_0                 Initial position of the reference meridiad of
%                               the planet at initial time T0 [deg]
%
% OUTPUT:
%
%   lat[Nx1]                latitude vector of orbit propagation [rad]
%   long[Nx1]               longitude vector of orbit propagation [rad]
%
% VERSIONS:
%   2025-11-25        First Version
%   2025-11-28        Added boundary check to remove boundary connecting line
%   2025-12-19      Bug-fix (boundary check section refined for cleaner plot)
%   2026-01-03: Last version
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti

    r_norm = vecnorm(y(:,1:3),2,2);
    delta = y(:,3) ./ r_norm;
    alpha = atan2(y(:,2),y(:,1));


    theta_g = deg2rad(theta_0) + deg2rad(w_e) / 3600 .* (t - t(1));

    in_long = wrapToPi(alpha - theta_g');
    in_lat = wrapToPi(asin(delta)); 

    LON = [];
    LAT = [];

    for ii = 1 : length(in_long) - 1
        LON = [LON; in_long(ii)];
        LAT = [LAT; in_lat(ii)];

        % check if there's a boundary cross
        if abs(in_long(ii+1) - in_long(ii)) > pi
            % check if boundary is on left or right side
            if in_long(ii) > 0
                boundary = pi;
                next_long_adj = in_long(ii+1) + 2*pi;
            else
                boundary = -pi;
                next_long_adj = in_long(ii+1) - 2*pi;
            end

            % interpolation to find approximate LAT in boundary point
            % fraction of path between actual point and next
            frac = (boundary - in_long(ii)) / (next_long_adj - in_long(ii));
            lat_interp = in_lat(ii) + frac * (in_lat(ii+1) - in_lat(ii));

            % add interp point + NaN
            LON = [LON; boundary; NaN; -boundary];
            LAT = [LAT; lat_interp; NaN; lat_interp];
        end
    end

    LON = [LON; in_long(end)];
    LAT = [LAT; in_lat(end)];

    long = LON;
    lat = LAT;
end

% function [long,lat] = GroundTrack (y,t,w_e,theta_0)
% % y       = [r,v] vector, not keplerian element for now
% % t       = time vector
% % w_e     = rotation velocity of planet [deg/h]
% % theta_0 = initial position of orbit [deg]
% 
%     delta  = (y(:,3)./vecnorm(y(:,1:3),2,2));
%     alpha = atan2(y(:,2),y(:,1));
%     theta_g = deg2rad(theta_0) + deg2rad(w_e) / (3600) .* (t-t(1));
%     long = alpha - theta_g';
%     lat = asin(delta);
% 
% 
%     in_long = wrapToPi(long);
%     in_lat = wrapToPi(lat);
% 
%     check = [double(abs(diff(in_long)) > deg2rad(180));0];
%     k = [0; find(check)];
%     LON = [];
%     LAT = [];
%     for ii = 2 : length(k)
%         LON = [LON; in_long(k(ii-1)+1 : k(ii)); NaN];
%         LAT = [LAT; in_lat(k(ii-1)+1:k(ii)); NaN  ];
%     end
%     LON = [LON; in_long(k(end)+1 : end)];
%     LAT = [LAT; in_lat(k(end)+1 : end)];
%     long = LON;
%     lat = LAT;
% 
% end