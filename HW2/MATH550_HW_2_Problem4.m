clc; clear; close all;

%% 1. DOWNLOAD LONGS PEAK TOPOGRAPHY DATA
fprintf('Downloading Longs Peak elevation data from USGS...\n');

% Longs Peak Center in UTM Zone 13N (Meters)
utm_easting_center  = 447600; % Easting (m)
utm_northing_center = 4456300; % Northing (m)

lat_min = 40.2505; lat_max = 40.2595;
lon_min = -105.6209; lon_max = -105.6093;

temp_tif = [tempname, '.tif'];
try
    websave(temp_tif, sprintf(['https://elevation.nationalmap.gov/arcgis/rest/services/' ...
        '3DEPElevation/ImageServer/exportImage?bbox=%.4f,%.4f,%.4f,%.4f' ...
        '&bboxSR=4326&size=200,200&format=tiff&f=image'], lon_min, lat_min, lon_max, lat_max));
    % Read raster directly. Imread returns Z(row, col) = Z(y, x).
    Z_elevation = double(imread(temp_tif));
    delete(temp_tif);
catch
    fprintf('USGS API unreachable. Generating synthetic Longs Peak profile...\n');
    [Xgrid, Ygrid] = meshgrid(linspace(-1,1,200));
    Z_elevation = 4346 * exp(-2*(Xgrid.^2 + Ygrid.^2)) + 50*cos(8*Xgrid).*sin(8*Ygrid);
end

%% 2. EXTRACT STRICTLY CLOSED CONTOURS
Z_smooth = imgaussfilt(Z_elevation, 1.5);
z_min = min(Z_smooth(:));
z_max = max(Z_smooth(:));

num_contours = 4; % Changed to 3 contours
[x_raw, y_raw, x_scaled, y_scaled] = deal(cell(num_contours, 1));
margin = 3.0; 

% STEP A: Find lowest closed contour
valid_base_found = false;
test_frac = 0.20;
while test_frac < 0.95 && ~valid_base_found
    z_test = z_min + test_frac*(z_max - z_min);
    C = contourc(1:size(Z_smooth,2), 1:size(Z_smooth,1), Z_smooth, [z_test z_test]);
    
    idx = 1; longest_len = 0; best_contour = [];
    while idx < size(C,2)
        num_pts = C(2, idx);
        segment = C(:, idx+1 : idx+num_pts);
        
        if min(segment(1,:)) > margin && max(segment(1,:)) < size(Z_smooth,2) - margin && ...
           min(segment(2,:)) > margin && max(segment(2,:)) < size(Z_smooth,1) - margin
            if num_pts > longest_len
                longest_len = num_pts;
                best_contour = segment;
            end
        end
        idx = idx + num_pts + 1;
    end
    
    if ~isempty(best_contour)
        valid_base_found = true;
        base_frac = test_frac;
    else
        test_frac = test_frac + 0.02; 
    end
end

% STEP B: Distribute elevations evenly
top_frac = min(0.95, base_frac + 0.4); 
if (top_frac - base_frac) < 0.15, top_frac = base_frac + 0.15; end
fractions = linspace(base_frac, top_frac, num_contours);
target_elevations_m = zeros(num_contours, 1);

for c = 1:num_contours
    z_test = z_min + fractions(c)*(z_max - z_min);
    target_elevations_m(c) = z_test; % True physical elevation in meters
    
    C = contourc(1:size(Z_smooth,2), 1:size(Z_smooth,1), Z_smooth, [z_test z_test]);
    idx = 1; longest_len = 0; best_contour = [];
    while idx < size(C,2)
        num_pts = C(2, idx);
        segment = C(:, idx+1 : idx+num_pts);
        
        if min(segment(1,:)) > margin && max(segment(1,:)) < size(Z_smooth,2) - margin && ...
           min(segment(2,:)) > margin && max(segment(2,:)) < size(Z_smooth,1) - margin
            if num_pts > longest_len
                longest_len = num_pts;
                best_contour = segment;
            end
        end
        idx = idx + num_pts + 1;
    end
    
    x_c = best_contour(1,:)'; y_c = best_contour(2,:)';
    if norm([x_c(1)-x_c(end), y_c(1)-y_c(end)]) > 1e-3
        x_c(end+1) = x_c(1); y_c(end+1) = y_c(1);
    end
    x_raw{c} = x_c; y_raw{c} = y_c;
end

scaled_elevs = target_elevations_m; 

% STEP C: Scale domain globally across ALL contours
all_x = cat(1, x_raw{:}); all_y = cat(1, y_raw{:});
center_x = (min(all_x) + max(all_x)) / 2;
center_y = (min(all_y) + max(all_y)) / 2;

max_radius = max(max(abs(all_x - center_x)), max(abs(all_y - center_y)));
scale_factor = (pi * 0.75) / max_radius;

for c = 1:num_contours
    x_scaled{c} = pi + (x_raw{c} - center_x) * scale_factor;
    y_scaled{c} = pi + (y_raw{c} - center_y) * scale_factor;
end

%% 3 & 4. FOURIER CONVERSION & IB PDE SETUP
font = 24; lw = 4.5; lwlrg = 6; lwsml = 3;
window_size_for_export = [0, 50, 1920, 975];
a_s = 1.2; 
[xdomain, ydomain] = deal([0, 2*pi]);

N = 200; 
h = (xdomain(2) - xdomain(1)) / (N - 1);
x = linspace(xdomain(1), xdomain(2), N);
y = linspace(ydomain(1), ydomain(2), N);
[X, Y] = meshgrid(x, y);

% Map PDE grid [0, 2*pi] to 1km x 1km UTM coordinates
x_utm = utm_easting_center  + (x - pi) * (1000 / (2*pi));
y_utm = utm_northing_center + (y - pi) * (1000 / (2*pi));
[X_utm, Y_utm] = meshgrid(x_utm, y_utm);

theta = linspace(0, 2*pi, max(400, floor(10*pi/h)));
num_harmonics = 12;

xib = []; yib = []; f2 = [];
xib_plot = cell(num_contours, 1); yib_plot = cell(num_contours, 1);
xib_utm  = cell(num_contours, 1); yib_utm  = cell(num_contours, 1);

for c = 1:num_contours
    [xc, yc, a_cf, b_cf, c_cf, d_cf] = contour_to_fourier(x_scaled{c}, y_scaled{c}, num_harmonics);
    [x_full, y_full] = generate_fourier_epicycle(theta, xc, yc, a_cf, b_cf, c_cf, d_cf);
    
    xib_plot{c} = x_full(1:end-1)'; yib_plot{c} = y_full(1:end-1)';
    xib = [xib; xib_plot{c}]; yib = [yib; yib_plot{c}];
    
    xib_utm{c} = utm_easting_center  + (xib_plot{c} - pi) * (1000 / (2*pi));
    yib_utm{c} = utm_northing_center + (yib_plot{c} - pi) * (1000 / (2*pi));
    
    f2 = [f2; ones(length(xib_plot{c}), 1) * scaled_elevs(c)];
end

%% 5. FAST VECTORIZED MATRIX ASSEMBLY
fprintf('Assembling Matrices...\n');
tic;

row_L = zeros(N^2 * 5, 1); col_L = zeros(N^2 * 5, 1); val_L = zeros(N^2 * 5, 1);
idx_L = 1;
for j = 1:N
    for i = 1:N
        g = global_index(j, i, N);
        if i == 1 || j == 1 || i == N || j == N
            row_L(idx_L) = g; col_L(idx_L) = g; 
            val_L(idx_L) = 1 / h^2; 
            idx_L = idx_L + 1;
        else
            row_L(idx_L:idx_L+4) = g;
            col_L(idx_L:idx_L+4) = [g, global_index(j, i+1, N), global_index(j, i-1, N), ...
                                       global_index(j+1, i, N), global_index(j-1, i, N)];
            val_L(idx_L:idx_L+4) = [-4, 1, 1, 1, 1] / h^2; 
            idx_L = idx_L + 5;
        end
    end
end
Lap = sparse(row_L(1:idx_L-1), col_L(1:idx_L-1), val_L(1:idx_L-1), N^2, N^2);

row_S = []; col_S = []; val_S = [];
r_cut = 4 * a_s * h;
const_factor = 1 / (2 * pi * (a_s*h)^2) * h^2;
X_vec = X(:); Y_vec = Y(:);

for k = 1:length(xib)
    idx = find(abs(X_vec - xib(k)) <= r_cut & abs(Y_vec - yib(k)) <= r_cut);
    if ~isempty(idx)
        dx = xib(k) - X_vec(idx); dy = yib(k) - Y_vec(idx);
        dist_sq = (dx.^2 + dy.^2) / (a_s*h)^2;
        vals = const_factor * exp(-0.5 * dist_sq);
        
        row_S = [row_S; idx];
        col_S = [col_S; repmat(k, length(idx), 1)];
        val_S = [val_S; vals];
    end
end
Sib = sparse(row_S, col_S, val_S, N^2, length(xib));
Jib = sparse(col_S, row_S, val_S, length(xib), N^2);

% Set base boundary condition on domain exterior
RHS_Lap = zeros(N^2, 1);
base_elevation_m = min(target_elevations_m) - 200; 
for j = 1:N
    for i = 1:N
        if i == 1 || j == 1 || i == N || j == N
            RHS_Lap(global_index(j, i, N)) = base_elevation_m / h^2;
        end
    end
end

reg_param = 1e-3; 
LHSMat = [Lap, -Sib; Jib, -reg_param * speye(length(xib))];
RHSVec = [RHS_Lap; f2];
toc;

%% 6. DIRECT LINEAR SOLVE
fprintf('Solving sparse matrix system...\n');
tic;
ResVec = LHSMat \ RHSVec; 
u = reshape(ResVec(1:N^2), N, N);
toc;

%% 7. VISUALIZATION
figure(1); set(gcf, 'Position', window_size_for_export);

% Plot 1: 3D Surface Solution with UTM and Elevation
subplot(1, 2, 1);
surf(X_utm, Y_utm, u, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
hold on; 
colormap(slanCM('plasma')); % Applied plasma colormap

for c = 1:num_contours
    plot3(xib_utm{c}, yib_utm{c}, ones(size(xib_utm{c})) * scaled_elevs(c), 'r-', 'LineWidth', lwlrg/2);
end

xlabel('UTM Easting [m]', 'Interpreter', 'latex'); 
ylabel('UTM Northing [m]', 'Interpreter', 'latex');
zlabel('Elevation [m]', 'Interpreter', 'latex');
title(sprintf('PDE Topography Solution ($N=%d$)', N), 'Interpreter', 'latex');
fontsize(font, "points");
set(gca, 'TickLabelInterpreter', 'latex');
ax1 = gca; ax1.LineWidth = lwsml*.5;
view(-37.5, 30); % Restored 3D isometric view

% Plot 2: Top-down view of contour fit comparison
subplot(1, 2, 2);
hold on;
for c = 1:num_contours
    x_raw_utm = utm_easting_center + (x_scaled{c} - pi) * (1000 / (2*pi));
    y_raw_utm = utm_northing_center + (y_scaled{c} - pi) * (1000 / (2*pi));
    
    p1 = plot(x_raw_utm, y_raw_utm, 'k--', 'LineWidth', lwsml*.75); 
    p2 = plot(xib_utm{c}, yib_utm{c}, 'r-', 'LineWidth', lw*.75);
end

plot(utm_easting_center, utm_northing_center, 'bx', 'MarkerSize', 15, 'LineWidth', 3);

grid on; axis equal;
xlabel('UTM Easting [m]', 'Interpreter', 'latex'); 
ylabel('UTM Northing [m]', 'Interpreter', 'latex');
title('Contours vs. Fourier Epicycles', 'Interpreter', 'latex');
legend([p1, p2], {'Raw USGS Contours', sprintf('Fourier (%d Harmonics)', num_harmonics)}, 'Interpreter', 'latex', 'Location', 'northeast');
fontsize(font, "points");
set(gca, 'TickLabelInterpreter', 'latex');
ax2 = gca; ax2.LineWidth = lwsml*.5;

%% 8. HELPER FUNCTIONS
function g = global_index(j, i, N)
    g = N*(j-1) + i;
end

function [xc, yc, a_coeffs, b_coeffs, c_coeffs, d_coeffs] = contour_to_fourier(x_contour, y_contour, num_harmonics)
    x_contour = x_contour(:); y_contour = y_contour(:);
    M = length(x_contour);
    t = linspace(0, 2*pi, M+1)'; t(end) = [];
    
    xc = mean(x_contour); yc = mean(y_contour);
    a_coeffs = zeros(1, num_harmonics); b_coeffs = zeros(1, num_harmonics);
    c_coeffs = zeros(1, num_harmonics); d_coeffs = zeros(1, num_harmonics);
    
    for k = 1:num_harmonics
        a_coeffs(k) = (2/M) * sum(x_contour .* cos(k*t));
        b_coeffs(k) = (2/M) * sum(x_contour .* sin(k*t));
        c_coeffs(k) = (2/M) * sum(y_contour .* cos(k*t));
        d_coeffs(k) = (2/M) * sum(y_contour .* sin(k*t));
    end
end

function [x_out, y_out] = generate_fourier_epicycle(theta, xc, yc, a, b, c, d)
    x_out = xc * ones(size(theta)); y_out = yc * ones(size(theta));
    max_k = max([length(a), length(b), length(c), length(d)]);
    a(end+1:max_k) = 0; b(end+1:max_k) = 0; 
    c(end+1:max_k) = 0; d(end+1:max_k) = 0;
    
    for k = 1:max_k
        x_out = x_out + a(k)*cos(k*theta) + b(k)*sin(k*theta);
        y_out = y_out + c(k)*cos(k*theta) + d(k)*sin(k*theta);
    end
end

exportgraphics(gcf,"HW2_4.jpg")