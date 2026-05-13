function tle2kep(fname)
%tle2kep TLE parsing and conversion to Keplerian elements time histories
%
% PROTOTYPE:
% tle2kep(fname)
%
% DESCRIPTION:
% This function reads a Two-Line Element (TLE) file and extracts the main
% orbital elements as time series. After filtering valid TLE line pairs and
% matching satellite numbers, the routine parses the epoch and line-2
% elements (inclination, RAAN, eccentricity, argument of perigee, mean
% motion). The semi-major axis is then derived from the mean motion using
% Kepler's third law with Earth's gravitational parameter. The resulting
% element histories are cleaned (NaNs removed), sorted by epoch, duplicates
% discarded, and angular elements are unwrapped to obtain continuous trends.
% Finally, the time evolution of the elements is plotted.
%
% INPUT:
% fname[char/string] Path to the input TLE text file containing lines in the
%                    standard "1 …" and "2 …" format [-]
%
% OUTPUT:
% None (the function generates figures of orbital elements versus time)
%
%
% CONTRIBUTORS:
% Andrea Accogli
% Fabrizio Pio Bavetta
% Sophia Trivellin
% Stefano Zotti
%
% VERSIONS:
% 2026-01-03

lines = readlines(fname);
lines = strtrim(lines);
lines(lines=="") = [];        % remove empty

% keep only proper TLE lines
lines = lines(startsWith(lines,"1 ") | startsWith(lines,"2 "));

L1 = lines(startsWith(lines,"1 "));
L2 = lines(startsWith(lines,"2 "));

N = min(numel(L1), numel(L2));
L1 = L1(1:N); 
L2 = L2(1:N);

% optional but recommended: keep only pairs with same satellite number
sat1 = arrayfun(@(s) str2double(extractBetween(s,3,7)), L1); % cols 3-7
sat2 = arrayfun(@(s) str2double(extractBetween(s,3,7)), L2);
ok = (sat1 == sat2) & isfinite(sat1);

L1 = L1(ok); L2 = L2(ok);

muE = astroConstants(13); % km^3/s^2
N = numel(L1);

epochs = NaT(N,1);
inc  = zeros(N,1);
RAAN = zeros(N,1);
e    = zeros(N,1);
argp = zeros(N,1);
nrev = zeros(N,1);

for k = 1:N
    a1 = char(L1(k));
    a2 = char(L2(k));

    % Epoch: YYDDD.DDDDDDDD at cols 19-32
    yy  = sscanf(a1(19:20), '%d');
    doy = sscanf(a1(21:32), '%f');
    if yy < 57, year = 2000+yy; else, year = 1900+yy; end
    epochs(k) = datetime(year,1,1) + days(doy-1);

    % Line 2: standard columns
    inc(k)  = sscanf(a2(9:16),  '%f');           % deg
    RAAN(k) = sscanf(a2(18:25), '%f');           % deg
    e(k)    = sscanf(['0.' a2(27:33)], '%f');    % -
    argp(k) = sscanf(a2(35:42), '%f');           % deg
    nrev(k) = sscanf(a2(53:63), '%f');           % rev/day
end

% Semi-major axis from mean motion
nrad = nrev * 2*pi/86400;
a = (muE ./ nrad.^2).^(1/3); % km

% Remove NaNs / invalid
good = isfinite(a) & isfinite(e) & isfinite(inc) & isfinite(RAAN) & isfinite(argp) & isfinite(nrev);
epochs = epochs(good); a=a(good); e=e(good); inc=inc(good); RAAN=RAAN(good); argp=argp(good);

% Sort by time (SUPER important)
[epochs, idx] = sort(epochs);
a=a(idx); e=e(idx); inc=inc(idx); RAAN=RAAN(idx); argp=argp(idx);

% Remove exact duplicate epochs (keep the last)
[epochs, ia] = unique(epochs,'last');
a=a(ia); e=e(ia); inc=inc(ia); RAAN=RAAN(ia); argp=argp(ia);

% Unwrap angles AFTER sorting
RAAN_u = unwrap(deg2rad(RAAN));
argp_u = unwrap(deg2rad(argp));

figure; 
subplot(3, 2, [1 2]); plot(epochs, a); grid on; ylabel('a [km]'); title('$a$(t) from TLE' , Interpreter='latex');
subplot(3, 2, 3); plot(epochs, e); grid on; ylabel('e [-]');  title('$e$(t) from TLE' , Interpreter='latex');
subplot(3, 2, 4); plot(epochs, inc); grid on; ylabel('i [deg]'); title('$i$(t) from TLE' , Interpreter='latex');
subplot(3, 2, 5); plot(epochs, rad2deg(RAAN_u)); grid on; ylabel('\Omega [deg]'); title('$\Omega$(t) from TLE' , Interpreter='latex');
subplot(3, 2, 6); plot(epochs, rad2deg(argp_u)); grid on; ylabel('\omega [deg]'); title('$\omega(t)$ from TLE' , Interpreter='latex');
[~, name, ~] = fileparts(fname);
sgtitle(sprintf("%s orbital element" , name),'interpreter','none');

end