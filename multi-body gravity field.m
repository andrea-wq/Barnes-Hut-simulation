clear
clc
close all

%#codegen

addpath("function\")
M_sun = 1.989e30;
G = 6.67430e-11;
[points,points_mass,points_j2,v_0,L,sizes,N,dt,t_end] = solar_sys();
% [points,points_mass,v_0,points_j2,sizes] = add_moon(points, points_mass , v_0 , points_j2,sizes);
% N=55;
% 
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

[points,points_mass,points_j2,v_0,L,sizes,N,dt,t_end] = galaxy_cluster(4);
% points = [points ; -4e12,-4e12,0];
% points_mass = [points_mass ; 10*M_sun];
% v_0 = [v_0 ; 0,0,0];
% sizes = [sizes , 150];
% N=10;

figure;
h = scatter3(points(:,1), points(:,2), points(:,3),0.5.*sizes, 'filled');
grid on;
axis equal;
hold off

xlim([-L L]);
ylim([-L L]);
zlim([-L L]);

ax = gca;
ax.Clipping = 'off';

trail_length = 75;
history = zeros(trail_length, 3, N);
for i = 1:N
    history(:,:,i) = repmat(points(i,:), trail_length, 1);
end

hold on;
h_trail = cell(N,1);  
for i = 1:N
    h_trail{i} = plot3(points(i,1), points(i,2), points(i,3), '-', 'LineWidth', 0.5);
end


node = octree_sys([0,0,0],L,points,points_mass,points_j2);

% Draw octree box wireframes in white with semi-transparent fills
figure
scatter3(points(:,1), points(:,2), points(:,3), 0.5.*sizes, 'filled');
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
    v_0 = v_0 + a * dt/2;          % secondo mezzo step

    h.XData = points(:,1);
    h.YData = points(:,2);
    h.ZData = points(:,3);

    history = circshift(history, 1, 1);
    history(1,:,:) = points';
    for i = 1:N
        h_trail{i}.XData = squeeze(history(:,1,i));
        h_trail{i}.YData = squeeze(history(:,2,i));
        h_trail{i}.ZData = squeeze(history(:,3,i));
    end

    % % Draw force vectors for each point using F (force per point)
    % % Remove previous force quiver graphics if they exist
    % if exist('h_force','var') && isgraphics(h_force)
    %     delete(h_force);
    % end
    % % Scale factor for visual clarity (adjust as needed)
    % scale = 1e-12;
    % % Start positions are the point coordinates
    % X = points(:,1);
    % Y = points(:,2);
    % Z = points(:,3);
    % U = F(:,1) * scale;
    % V = F(:,2) * scale;
    % W = F(:,3) * scale;
    % % Draw 3D quiver for forces and store handle in persistent variable
    % h_force = quiver3(X, Y, Z, U, V, W, 0, 'Color', [0 0.6 0], 'LineWidth', 1);

    drawnow limitrate;
    pause(0.01);
    t = t + dt;
end


% TESTs
% node = octree([L/2, L/2, L/2],L/2,points,points_mass);
% 
% node = node_mass(node);
% 
% % Draw octree box wireframes in white with semi-transparent fills
% hold on;
% axis equal;
% view(3);
% 
% drawNode(node);
% hold off;
% 
% 
% 
% 
% for i = 1:length(points)
% 
%     F(i,:) = calc_forza(node,points(i,:),points_mass(i));
% end
% 
% dt = 1000000;
% 
% for i = 1:length(points)
% 
%     v(i,:) = v_0(i,:) + F(i,:)/points_mass(i) *dt;
% 
%     points(i,:) = points(i,:) + v(i,:)*dt;
% 
%     v_0(i,:) = v(i,:);
% end



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
    L = 3e16;
    M_sun = 1.989e30;
    
    R_disk = L * 0.3;
    r      = R_disk * sqrt(rand(N,1));   
    phi    = 2*pi * rand(N,1);
    
    x = r .* cos(phi);
    y = r .* sin(phi);
    z = randn(N,1) * R_disk * 0.05;     % disco sottile ~5% dello spessore radiale
    
    points = [x, y, z] + [L/2, L/2, L/2];  % centra nel box
    points_mass = ((0.3 + 9.9*rand(N,1)) * M_sun);
    
    G     = 6.67430e-11;
    M_tot = sum(points_mass);
    
    v_circ = sqrt(G * M_tot ./ r);      % velocità orbitale kepleriana
    % direzione tangenziale nel piano xy (perpendicolare a r)
    v_0 = zeros(N,3);
    v_0(:,1) = -v_circ .* sin(phi) * 0.3;
    v_0(:,2) =  v_circ .* cos(phi) * 0.3;
    v_0(:,3) = 0;                        % nessuna componente z

    sizes = 20/0.5;

    points = points - [L/2,L/2,L/2];


    v_typ = 10e3;                % velocità tipica, m/s
    t_cross = L / v_typ;         % ~3e12 s ~ 100k anni
    dt    = t_cross / 100;   % ~3e10 s, ~1000 anni  — risolve bene le orbite
    t_end = 50 * t_cross;    % ~3e13 s, ~1M anni    — vedi evoluzione globale

    points_j2 = zeros(size(points,1),size(points,2));

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