clc
close all

N = 5;
L = 3e16;
M_sun = 1.989e30;
sizeng = 0.5e-29;


R_disk = L * 0.3;
r      = R_disk * sqrt(rand(N,1));   
phi    = 2*pi * rand(N,1);

x = r .* cos(phi);
y = r .* sin(phi);
z = randn(N,1) * R_disk * 0.05;     % disco sottile ~5% dello spessore radiale

points = [x, y, z] + [L/2, L/2, L/2];  % centra nel box
points_mass = ((0.3 + 9.9*rand(N,1)) * M_sun)';

G     = 6.67430e-11;
M_tot = sum(points_mass);

v_circ = sqrt(G * M_tot ./ r);      % velocità orbitale kepleriana
% direzione tangenziale nel piano xy (perpendicolare a r)
v_0 = zeros(N,3);
v_0(:,1) = -v_circ .* sin(phi) * 0.3;
v_0(:,2) =  v_circ .* cos(phi) * 0.3;
v_0(:,3) = 0;                        % nessuna componente z

figure;
%%
% Optionally plot the points
h = scatter3(points(:,1), points(:,2), points(:,3), sizeng*points_mass, 'filled');
grid on;
axis equal;
hold off

xlim([0 L]);
ylim([0 L]);
zlim([0 L]);


trail_length = 150;
history = zeros(trail_length, 3, N);
for i = 1:N
    history(:,:,i) = repmat(points(i,:), trail_length, 1);
end

hold on;
h_trail = cell(N,1);  
for i = 1:N
    h_trail{i} = plot3(points(i,1), points(i,2), points(i,3), '-', 'LineWidth', 0.5);
end

%% Meh draw
node = octree([L/2, L/2, L/2],L/2,points,points_mass);

% Draw octree box wireframes in white with semi-transparent fills
hold on;
axis equal;
view(3);

drawNode(node);
hold off;

%%

v_typ = 10e3;                % velocità tipica, m/s
t_cross = L / v_typ;         % ~3e12 s ~ 100k anni
dt    = t_cross / 100;   % ~3e10 s, ~1000 anni  — risolve bene le orbite
t_end = 10 * t_cross;    % ~3e13 s, ~1M anni    — vedi evoluzione globale

t=0;
while t < t_end

    [node,conto] = octree([L/2, L/2, L/2],L/2,points,points_mass);
    node = node_mass(node);

    F = zeros(N,3);
    for i = 1:size(points,1)

        F(i,:) = calc_forza(node,points(i,:),points_mass(i));
    end
    
    v = v_0 + F./points_mass' * dt;
    points = points + v*dt;
    v_0 = v;

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

    drawnow limitrate;
    pause(0.01);
    t = t + dt;
end


%% TESTs
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

%%
function [node,count] = octree(center,half,point,point_mass)
count = 1;
centers = [-1 -1 -1
           -1 -1 1
            -1 1 -1
            -1 1 1
            1 -1 -1
            1 -1 1
            1 1 -1
            1 1 1];

    node.center = center;
    node.half = half;
    node.point = point;
    node.point_mass = point_mass;
    node.child = {};

if size(point,1) <= 1  % end recursion
    return;
end

    node.child = cell(1,8);
    
    bin = point > center;
    idx = 1 + bin(:,1)*4 + bin(:,2)*2 + bin(:,3);

    for i = 1:8
        
        child_pts = point(idx==i,:);
        child_pts_mass = point_mass(idx==i);

        if isempty(child_pts)
            continue
        end
        child_center = center + half/2 * centers(i,:);
        [node.child{i},c] = octree(child_center,half/2,child_pts,child_pts_mass);
        count = c+count;
    end

end



function node = node_mass(node)

    node.mass = sum(node.point_mass);

    x_center = sum(node.point(:,1).*node.point_mass') / node.mass;
    y_center = sum(node.point(:,2).*node.point_mass') / node.mass;
    z_center = sum(node.point(:,3).*node.point_mass') / node.mass;

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



% recursive function to draw boxes
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
    eps_soft = 1e13;

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
    else
        for i = 1:8
            if ~isempty(node.child{i})
                f = f + calc_forza(node.child{i}, point, point_mass); % Accumulate forces from child nodes
            end
        end
    end

end