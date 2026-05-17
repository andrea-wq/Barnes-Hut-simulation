clear
clc
close all

addpath("function\")
M_sun = 1.989e30;
G = 6.67430e-11;
[points,points_mass,points_j2,v_0,L,sizes,N,dt,t_end] = solar_sys();
% [points,points_mass,v_0,points_j2,sizes] = add_moon(points, points_mass , v_0 , points_j2,sizes);
% N=55;
% idx_moon = 11:N;
% idx_planet = 1:10;

% keepIdx = [3, 10, 11];
% % Ensure indices are within current range
% keepIdx = keepIdx(keepIdx >= 1 & keepIdx <= size(points,1));
% 
% points = points(keepIdx, :);
% points_mass = points_mass(keepIdx, :);
% v_0 = v_0(keepIdx, :);
% points_j2 = points_j2(keepIdx, :);
% sizes = sizes(keepIdx);
% N = numel(keepIdx);

%%
clear
clc
close all
addpath("function\")
M_sun = 1.989e30;
[points,points_mass,points_j2,v_0,L,sizes,N,dt,t_end] = galaxy_cluster(4);
v_0 = zeros(4,3);
% points = [points ; 0,0,0];
% points_mass = [points_mass ; 300*M_sun];
% points_j2 = [points_j2;0];
% v_0 = [v_0;0,0,0];
% sizes = [sizes , 150];
% N = N+1;

% points = [points ; -4e12,-4e12,0];
% points_mass = [points_mass ; 10*M_sun];
% v_0 = [v_0 ; 0,0,0];
% sizes = [sizes , 150];
% N=10;

figure;
h = scatter3(points(:,1), points(:,2), points(:,3),3.*sizes, 'filled');
grid on;
axis equal;
hold off

xlim([-L L]);
ylim([-L L]);
zlim([-L L]);

ax = gca;
ax.Clipping = 'off';

% trail_length = 75;
% history = zeros(trail_length, 3, N);
% for i = 1:N
%     history(:,:,i) = repmat(points(i,:), trail_length, 1);
% end
% 
% hold on;
% h_trail = cell(N,1);  
% for i = 1:N
%     h_trail{i} = plot3(points(i,1), points(i,2), points(i,3), '-', 'LineWidth', 0.5);
% end


node = octree_sys([0,0,0],L,points,points_mass,points_j2);

% Draw octree box wireframes in white with semi-transparent fills
figure
scatter3(points(:,1), points(:,2), points(:,3), 0.5*sizes, 'filled');
grid on;
axis equal;
hold off

xlim([-L L]);
ylim([-L L]);
zlim([-L L]);

ax = gca;
ax.Clipping = 'off';
hold on;
axis equal;
view(3);

drawNode(node);
hold off;

%%
n_planets = 10;
n_sub = 100;        % abbassa se troppo lento, non scendere sotto 20
dt_sub = dt / n_sub;
G = 6.67430e-11;

dt = dt/4;

% T_orb = 225 * 1e6 * 3.156e7;   % ~7.1e15 s
% dt    = T_orb / 200;            % ~3.5e13 s ≈ 1.1 Myr — buona risoluzione
% t_end = 20 * T_orb;            % ~4.5 Gyr, una rotazione galattica completa

% Verifica energia iniziale luna-terra
% i_luna = n_planets + 1;  % indice Luna nel vettore
% i_terra = 3;             % Terra
% 
% r = points(i_terra,:) - points(i_luna,:);
% v_rel = v_0(i_terra,:) - v_0(i_luna,:);
% 
% E_kin = 0.5 * points_mass(i_luna) * dot(v_rel, v_rel);
% E_pot = -G * points_mass(i_terra) * points_mass(i_luna) / norm(r);
% E_tot = E_kin + E_pot;
% 
% fprintf('E_cin  = %.4e J\n', E_kin);
% fprintf('E_pot  = %.4e J\n', E_pot);
% fprintf('E_tot  = %.4e J\n', E_tot);
% fprintf('Bound? %s\n', string(E_tot < 0));

t=0;
while t < t_end

    [node,conto] = octree_sys([0,0,0],L,points,points_mass,points_j2);
    node = node_mass(node);

    F = zeros(N,3);
    for i = 1:size(points,1)
        F(i,:) = calc_forza(node,points(i,:),points_mass(i));
    end
    
    a = F ./ points_mass;
    v_0 = v_0 + a * dt/2;          % mezzo step velocità
    points = points + v_0 * dt;    % step posizione

    % ricalcola F con nuove posizioni
    node = octree_sys([0,0,0],L,points,points_mass,points_j2);
    node = node_mass(node);
    F = zeros(N,3);
    for i = 1:N
        F(i,:) = calc_forza(node, points(i,:), points_mass(i));
    end

    a = F ./ points_mass;
    v_0 = v_0 + a * dt/2;         % secondo mezzo step

    % points(end,:) = [0,0,0];
    % v_0(end,:) = [0,0,0];


    % Update scatter positions without drawing force quivers
    if exist('h','var') && isgraphics(h)
        set(h, 'XData', points(:,1), 'YData', points(:,2), 'ZData', points(:,3));
    else
        h = scatter3(points(:,1), points(:,2), points(:,3), 0.5.*sizes, 'filled');
    end

    drawnow limitrate;
    pause(0.01);
    t = t + dt;
end


function [node, count] = octree(center, half, point, point_mass, depth)
    if nargin < 5, depth = 0; end
    count = 1;

    centers = [-1 -1 -1
               -1 -1  1
               -1  1 -1
               -1  1  1
                1 -1 -1
                1 -1  1
                1  1 -1
                1  1  1];

    node.center     = center;
    node.half       = half;
    node.point      = point;
    node.point_mass = point_mass;
    node.child      = {};

    % guardie di stop
    if size(point,1) <= 1 || depth > 25 || half < 1e17
        return;
    end

    % punti coincidenti → tratta come foglia
    if max(max(abs(point - point(1,:)))) < 1e-10 * half
        return;
    end

    % jitter sui punti esattamente sul bordo
    on_boundary = abs(point - center) < 1e-10 * half;
    point(on_boundary) = point(on_boundary) + 1e-9 * half;

    bin = point > center;
    idx = 1 + bin(:,1)*4 + bin(:,2)*2 + bin(:,3);

    node.child = cell(1,8);
    for i = 1:8
        child_pts      = point(idx == i, :);
        child_pts_mass = point_mass(idx == i);

        if isempty(child_pts), continue; end

        % bounding box del figlio con half minimo garantito
        child_center = center + half/2 * centers(i,:);
        child_half   = max(half/2, 1e3);   % mai zero

        [node.child{i}, c] = octree(child_center, child_half, child_pts, child_pts_mass, depth+1);
        count = count + c;
    end
end


function [node, count] = octree_sys(center, half, point, point_mass, point_j2 , depth)
    if nargin < 6, depth = 0; end
    count = 1;

    centers = [-1 -1 -1
               -1 -1  1
               -1  1 -1
               -1  1  1
                1 -1 -1
                1 -1  1
                1  1 -1
                1  1  1];

    node.center     = center;
    node.half       = half;
    node.point      = point;
    node.point_mass = point_mass;
    node.point_j2   = point_j2;
    node.child      = {};

    % guardie di stop
    if size(point,1) <= 1 || depth > 20 || half < 1e3
        return;
    end

    % punti coincidenti → tratta come foglia
    if max(max(abs(point - point(1,:)))) < 1e-10 * half
        return;
    end

    % jitter sui punti esattamente sul bordo
    on_boundary = abs(point - center) < 1e-10 * half;
    point(on_boundary) = point(on_boundary) + 1e-9 * half;

    bin = point > center;
    idx = 1 + bin(:,1)*4 + bin(:,2)*2 + bin(:,3);

    node.child = cell(1,8);
    for i = 1:8
        child_pts      = point(idx == i, :);
        child_pts_mass = point_mass(idx == i);
        child_pts_j2   = point_j2(idx == i);

        if isempty(child_pts), continue; end

        % bounding box del figlio con half minimo garantito
        child_center = center + half/2 * centers(i,:);
        child_half   = max(half/2, 1e3);   % mai zero

        [node.child{i}, c] = octree_sys(child_center, child_half, child_pts, child_pts_mass, child_pts_j2 , depth+1);
        count = count + c;
    end
end


function node = node_mass(node)

    node.mass = sum(node.point_mass);

    x_center = sum(node.point(:,1).*node.point_mass) / node.mass;
    y_center = sum(node.point(:,2).*node.point_mass) / node.mass;
    z_center = sum(node.point(:,3).*node.point_mass) / node.mass;

    node.cm = [x_center,y_center,z_center];

    if isempty(node.child)
        return;
    end
    
    for i =1:8
        if isempty(node.child{i})
            continue;
        end
        node.child{i} = node_mass(node.child{i});
    end

end


function drawNode(n)
    if isempty(n), return; end
    c = n.center;
    h = n.half;
    % 8 corners of the cube
    corners = c + h * [ -1 -1 -1
                        -1 -1  1
                        -1  1 -1
                        -1  1  1
                         1 -1 -1
                         1 -1  1
                         1  1 -1
                         1  1  1 ];
    % faces as indices into corners (each row is a quad)
    F = [1 2 4 3;  % -X face
         5 6 8 7;  % +X face
         1 2 6 5;  % -Y face
         3 4 8 7;  % +Y face
         1 3 7 5;  % -Z face
         2 4 8 6]; % +Z face

    % Plot filled faces with low opacity (10%)
    facecolor = [1 1 1]; % white
    for f = 1:size(F,1)
        verts = corners(F(f,:),:);
        patch('Vertices',corners,'Faces',F(f,:),'FaceColor',facecolor, ...
              'FaceAlpha',0.05,'EdgeColor','none');
    end

    % edges as pairs of corner indices
    E = [1 2;1 3;1 5;
         2 4;2 6;
         3 4;3 7;
         4 8;
         5 6;5 7;
         6 8;
         7 8];
    % plot edges in white
    for k = 1:size(E,1)
        p1 = corners(E(k,1),:);
        p2 = corners(E(k,2),:);
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'w-','LineWidth',0.8);
    end
    % recurse into children
    if ~isempty(n.child)
        for ii = 1:numel(n.child)
            if ~isempty(n.child{ii})
                drawNode(n.child{ii});
            end
        end
    end
end


function f = calc_forza(node,point,point_mass)
    theta = 0.5; % Define the opening angle for the calculation
    G = 6.67430e-11; % Gravitational constant
    eps_soft = 1e5;
    M_sun = 1.9890e30;

    f = [0,0,0];

    r_vec = node.cm - point;
    d2     = dot(r_vec, r_vec) + eps_soft^2;  % distanza ammorbidita al quadrato
    d      = sqrt(d2);

    if isempty(node.child) && size(node.point, 1) == 1
        if norm(node.point - point) < eps_soft
            return;
        end
    end
    
    s     = 2*node.half;

    if isempty(node.child) || s/d < theta
        f = G * point_mass * node.mass * r_vec / d^3; % Calculate gravitational force
        if isempty(node.child)
            ROI = norm(node.point) * (node.point_mass/M_sun)^(2/5);
            if d < ROI
                f = f + node.point_j2.*[r_vec(1)/(d^7) * (6*r_vec(3)^2 - 3/2*(r_vec(1)^2+r_vec(2)^2) )
                                        r_vec(1)/(d^7) * (6*r_vec(3)^2 - 3/2*(r_vec(1)^2+r_vec(2)^2) )
                                        r_vec(1)/(d^7) * (3*r_vec(3)^2 - 9/2*(r_vec(1)^2+r_vec(2)^2) )]';
            end
        end
    else
        for i = 1:8
            if ~isempty(node.child{i})
                f = f + calc_forza(node.child{i}, point, point_mass); % Accumulate forces from child nodes
            end
        end
    end

end


function [points,points_mass,point_j2,v_0,L,sizes,N,dt,t_end] = solar_sys()
    
    N = 10;
    M_sun = 1.989e30;
    radii = [2440, 6052, 6371, 3390, 69911, 58232, 25362, 24622, 2440 , 696000];  % km
    radii_norm = radii / max(radii);   % normalizza 0→1
    
    % comprimi con radice per ridurre il divario
    sizes = (radii_norm .^ 0.2) * 200;   % esponente < 1 comprime i grandi
    sizes(end) = 100;
    
    
    G = 6.67430e-11;
    for i = 1:9
        planet.kep(i,:) = uplanet(0,i);
        [planet.car(i,:),planet.vel(i,:)] = kep2car(planet.kep(i,:),astroConstants(4));
        planet.mass(i) = astroConstants(10+i)*10^9/G;
        point_j2(i)  = astroConstants(30 + i);
    end
    
    points = [planet.car*1000 ; 0,0,0] ;
    points_mass = [planet.mass , M_sun]';
    v_0 = [planet.vel ; 0,0,0]*1000;
    point_j2 = [point_j2, 0]';
    
    for i = 1:9
        dist(i) = norm(points(i,:));
    end
    
    L = max(dist);

    T_mercurio = 88 * 24 * 3600;   % ~7.6e6 s
    dt = T_mercurio / 100;
    t_end = 200*T_mercurio;

end


function [points,points_mass,points_j2,v_0,L,sizes,N,dt,t_end] = galaxy_cluster(N)

    L      = 3e21;
    M_sun  = 1.989e30;
    G      = 6.67430e-11;

    % --- Geometria disco con raggio minimo e massimo ---
    R_min  = L * 0.04;          % raggio interno (buco centrale)
    R_max  = L * 0.45;          % raggio esterno del disco
    
    % distribuzione uniforme in area tra R_min e R_max
    u = rand(N,1);
    r = sqrt(u * (R_max^2 - R_min^2) + R_min^2);

    phi = 2*pi * rand(N,1);
    z   = randn(N,1) * (R_max - R_min) * 0.02;   % spessore 2% del disco

    points = [r.*cos(phi), r.*sin(phi), z];

    % --- Masse: distribuzione di Salpeter ---
    m_min = 0.1 * M_sun;
    m_max = 50  * M_sun;
    u2 = rand(N,1);
    points_mass = (m_min^(-1.35) + u2*(m_max^(-1.35) - m_min^(-1.35))).^(1/-1.35);

    % --- Massa enclosed entro r (disco + alone) ---
    [r_sorted, sort_idx] = sort(r);
    M_enc_sorted = cumsum(points_mass(sort_idx));
    M_enc = zeros(N,1);
    M_enc(sort_idx) = M_enc_sorted;

    M_disk    = sum(points_mass);
    M_halo_enc = M_disk * (r / R_max);      % alone isotermo singolare

    v_circ = sqrt(G * (M_enc + M_halo_enc) ./ r);

    % Clamp di sicurezza
    v_max  = 1000e3;
    n_clip = sum(v_circ > v_max);
    if n_clip > 0
        warning('%d particelle clampate a 1000 km/s', n_clip);
        v_circ = min(v_circ, v_max);
    end

    % --- Velocità tangenziali + dispersione ---
    sigma_v = 0.05 * v_circ;
    v_0 = zeros(N,3);
    v_0(:,1) = -v_circ .* sin(phi) + sigma_v .* randn(N,1);
    v_0(:,2) =  v_circ .* cos(phi) + sigma_v .* randn(N,1);
    v_0(:,3) =                        sigma_v .* randn(N,1) * 0.5;

    % --- Parametri simulazione ---
    r_ref  = (R_min + R_max) / 2;
    M_ref  = M_disk/2 + M_disk*0.5;
    T_ref  = 2*pi * r_ref / sqrt(G * M_ref / r_ref);
    dt     = T_ref / 200;
    t_end  = 10 * T_ref;

    % --- Output ---
    L         = R_max * 2.5;
    sizes     = 5 * ones(N,1)';
    points_j2 = zeros(N,1);

    fprintf('R_min   = %.3e m = %.1f kpc\n', R_min, R_min/3.086e19);
    fprintf('R_max   = %.3e m = %.1f kpc\n', R_max, R_max/3.086e19);
    fprintf('dt      = %.3e s = %.1f Myr\n', dt,    dt   /(1e6*3.156e7));
    fprintf('t_end   = %.3e s = %.1f Gyr\n', t_end, t_end/(1e9*3.156e7));
    fprintf('v_circ  media = %.1f km/s\n',   mean(v_circ)/1e3);
end


function [point,points_mass,v_0,points_j2,sizes] = add_moon(points,points_mass,v_0,points_j2,sizes)
    moon_index = [301,  401,402,  501:505,514:516,  601:608,632,634,609,612,613,614  ,701:705,715,  801:808,814,  901:905];
    moon_masses = [7.342e22, 1.062e16, 1.441e15, 8.932e22, 4.800e22, 1.482e23, 1.076e23, 2.080e18, 4.300e17, 2.000e15, 3.600e16, 3.751e19, 1.080e20, 6.175e20, 1.096e21, 2.308e21, 1.346e23, 5.586e18, 1.806e21, 1.400e13, 3.000e13, 8.292e18, 2.500e15, 7.200e15, 3.600e15, 1.353e21, 1.172e21, 3.527e21, 3.014e21, 6.590e19, 2.900e18, 2.140e22, 3.090e19, 1.900e17, 3.500e17, 2.100e18, 2.120e18, 4.900e18, 4.400e19, 5.000e15, 1.587e21, 4.500e16, 4.800e16, 1.650e16, 7.500e15];
    point = points;
    for i = 1:length(moon_index)
        [point(size(points,1)+i,:),v_0(size(points,1)+i,:)] = horizon_data(0,moon_index(i));
        points_mass(size(points,1)+i) = moon_masses(i);
        points_j2(size(points,1)+i) = 0;
        sizes(size(points,1)+i) = 30;
    end
   pause(1)
end


function [pos,vel] = horizon_data(mjd2000_date , index)

url="https://ssd.jpl.nasa.gov/api/horizons.api";
center='500@ssb';

vec_table='2';

ref_system='ICRF';
ref_plane='ECLIPTIC';
command=sprintf("%d",index);

opts=weboptions('Timeout',60);

app = (mjd20002date(mjd2000_date));
start = sprintf("%d-%d-%d", app(1:3));
end_time = sprintf("%d-%d-%d", app(1) , app(2) , app(3)+1);

data=webread(url,'format','text', ...
                 'COMMAND',command, ...
                 'OBJ_DATA','NO', ...
                 'MAKE_EPHEM','YES',...
                 'EPHEM_TYPE','VECTORS', ...
                 'START_TIME',start, ...
                 'CENTER',center,...
                 'STOP_TIME',end_time, ...
                 'STEP_SIZE',1, ...
                 'QUANTITIES','A',...
                 'OUT_UNITS','KM-S', ...
                 'RANGE_UNITS','KM', ...
                 'VEC_TABLE',vec_table,...
                 'CSV_FORMAT','YES', ...
                 'REF_SYSTEM',ref_system, ...
                 'REF_PLANE',ref_plane,...
                 'ECLIP','B1950',opts);


% Find the ephemeris data block in the returned text and extract the CSV lines
soe_idx = strfind(data, '$$SOE');
eoe_idx = strfind(data, '$$EOE');

if isempty(soe_idx) || isempty(eoe_idx)
    error('Could not find ephemeris start/end markers in data.');
end

% Extract the block between $$SOE and $$EOE (exclusive of markers)
eph_block = data(soe_idx + length('$$SOE') : eoe_idx - 1);

% Split into lines and trim whitespace
lines = regexp(eph_block, '\r\n|\r|\n', 'split');
lines = strtrim(lines);
% Remove empty lines
lines = lines(~cellfun('isempty', lines));

if isempty(lines)
    error('No ephemeris lines found in the data block.');
end

% The CSV from Horizons VECTORS typically has a header line followed by data lines.
% Find the first data line (skip header lines that contain non-numeric text like 'date' or 'UTC')
firstDataIdx = 1;
for k = 1:numel(lines)
    % Heuristic: data lines start with a quote (date) or numeric; header contains 'date' or 'utc'
    if isempty(regexpi(lines{k}, 'date|utc', 'once')) 
        firstDataIdx = k;
        break;
    end
end

% If the first line is a header (contains commas and 'date'), use the next line as first data line
if ~isempty(regexpi(lines{1}, 'date|utc', 'once'))
    if numel(lines) < 2
        error('CSV header present but no data line found.');
    end
    firstLine = lines{2};
else
    firstLine = lines{1};
end

% Split the first CSV data line into fields
fields = strsplit(firstLine, ',');
% Trim possible quotes/spaces
fields = strtrim(strrep(fields, '"', ''));

% Horizons VECTORS CSV format: date, JD, x, y, z, vx, vy, vz, ... (positions in km, velocities in km/s)
% Find numeric fields: positions are typically fields 3-5, velocities 6-8 (1-based)
if numel(fields) < 8
    error('Unexpected CSV format: not enough fields for position and velocity.');
end

pos = str2double(fields(3:5)) * 1000;
vel = str2double(fields(6:8)) * 1000;

if any(isnan(pos)) || any(isnan(vel))
    error('Failed to parse position/velocity from the first ephemeris line.');
end

% pos and vel are 1x3 vectors containing position (km) and velocity (km/s) respectively.

end


function F_m = moon_forces(points, masses, n_planets, idx_moons, G)
    n_m = numel(idx_moons);
    F_m = zeros(n_m, 3);
    eps = 1e5;  % stesso softening piccolo che hai già

    for i = 1:n_m
        mi = idx_moons(i);
        for j = 1:n_planets
            r = points(j,:) - points(mi,:);
            d = sqrt(dot(r,r) + eps^2);
            F_m(i,:) = F_m(i,:) + G * masses(mi) * masses(j) * r / d^3;
        end
    end
end