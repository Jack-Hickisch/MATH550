clc; clear; close all;

%% 1. DOWNLOAD LONGS PEAK TOPOGRAPHY DATA
% Longs Peak Center in UTM Zone 13N (Meters)
utm_easting_center  = 447600; 
utm_northing_center = 4456300; 
domain_length_m     = 1000; % 1km x 1km domain for physical PDE grid

lat_min = 40.2505; lat_max = 40.2595;
lon_min = -105.6209; lon_max = -105.6093;

temp_tif = [tempname, '.tif'];
try
    websave(temp_tif, sprintf(['https://elevation.nationalmap.gov/arcgis/rest/services/' ...
        '3DEPElevation/ImageServer/exportImage?bbox=%.4f,%.4f,%.4f,%.4f' ...
        '&bboxSR=4326&size=200,200&format=tiff&f=image'], lon_min, lat_min, lon_max, lat_max));
    Z_elevation = double(imread(temp_tif));
    delete(temp_tif);
catch
    [Xgrid, Ygrid] = meshgrid(linspace(-1,1,200));
    Z_elevation = 4346 * exp(-2*(Xgrid.^2 + Ygrid.^2)) + 50*cos(8*Xgrid).*sin(8*Ygrid);
end

%% 2. EXTRACT STRICTLY CLOSED CONTOURS
Z_smooth = imgaussfilt(Z_elevation, 1.5);
z_min = min(Z_smooth(:)); z_max = max(Z_smooth(:));

num_contours = 3;
[x_raw, y_raw, x_scaled, y_scaled] = deal(cell(num_contours, 1));
margin = 3.0; 

valid_base_found = false; test_frac = 0.20;
while test_frac < 0.95 && ~valid_base_found
    z_test = z_min + test_frac*(z_max - z_min);
    C = contourc(1:size(Z_smooth,2), 1:size(Z_smooth,1), Z_smooth, [z_test z_test]);
    
    idx = 1; longest_len = 0; best_contour = [];
    while idx < size(C,2)
        num_pts = C(2, idx); segment = C(:, idx+1 : idx+num_pts);
        if min(segment(1,:)) > margin && max(segment(1,:)) < size(Z_smooth,2) - margin && ...
           min(segment(2,:)) > margin && max(segment(2,:)) < size(Z_smooth,1) - margin
            if num_pts > longest_len, longest_len = num_pts; best_contour = segment; end
        end
        idx = idx + num_pts + 1;
    end
    if ~isempty(best_contour), valid_base_found = true; base_frac = test_frac;
    else, test_frac = test_frac + 0.02; end
end

top_frac = min(0.95, base_frac + 0.4); 
if (top_frac - base_frac) < 0.15, top_frac = base_frac + 0.15; end
fractions = linspace(base_frac, top_frac, num_contours);
target_elevations_m = zeros(num_contours, 1);

for c = 1:num_contours
    z_test = z_min + fractions(c)*(z_max - z_min);
    target_elevations_m(c) = z_test; 
    
    C = contourc(1:size(Z_smooth,2), 1:size(Z_smooth,1), Z_smooth, [z_test z_test]);
    idx = 1; longest_len = 0; best_contour = [];
    while idx < size(C,2)
        num_pts = C(2, idx); segment = C(:, idx+1 : idx+num_pts);
        if min(segment(1,:)) > margin && max(segment(1,:)) < size(Z_smooth,2) - margin && ...
           min(segment(2,:)) > margin && max(segment(2,:)) < size(Z_smooth,1) - margin
            if num_pts > longest_len, longest_len = num_pts; best_contour = segment; end
        end
        idx = idx + num_pts + 1;
    end
    
    x_c = best_contour(1,:)'; y_c = best_contour(2,:)';
    if norm([x_c(1)-x_c(end), y_c(1)-y_c(end)]) > 1e-3, x_c(end+1) = x_c(1); y_c(end+1) = y_c(1); end
    x_raw{c} = x_c; y_raw{c} = y_c;
end

all_x = cat(1, x_raw{:}); all_y = cat(1, y_raw{:});
center_x = (min(all_x) + max(all_x)) / 2; center_y = (min(all_y) + max(all_y)) / 2;
max_radius = max(max(abs(all_x - center_x)), max(abs(all_y - center_y)));
scale_factor = (pi * 0.75) / max_radius;

for c = 1:num_contours
    x_scaled{c} = pi + (x_raw{c} - center_x) * scale_factor;
    y_scaled{c} = pi + (y_raw{c} - center_y) * scale_factor;
end

%% 3 & 4. FOURIER CONVERSION & PHYSICAL UTM PDE SETUP
N = 80; 
x_utm = linspace(utm_easting_center - domain_length_m/2, utm_easting_center + domain_length_m/2, N);
y_utm = linspace(utm_northing_center - domain_length_m/2, utm_northing_center + domain_length_m/2, N);
h_utm = x_utm(2) - x_utm(1); 
[X_utm, Y_utm] = meshgrid(x_utm, y_utm);

num_harmonics = 12;
theta = linspace(0, 2*pi, max(400, floor(10*pi/(h_utm/(domain_length_m/(2*pi))))));

xib = []; yib = []; f2 = [];
xib_utm_cell = cell(num_contours, 1); yib_utm_cell = cell(num_contours, 1);

for c = 1:num_contours
    [xc, yc, a_cf, b_cf, c_cf, d_cf] = contour_to_fourier(x_scaled{c}, y_scaled{c}, num_harmonics);
    [x_full, y_full] = generate_fourier_epicycle(theta, xc, yc, a_cf, b_cf, c_cf, d_cf);
    
    xib_utm_c = utm_easting_center  + (x_full(1:end-1)' - pi) * (domain_length_m / (2*pi));
    yib_utm_c = utm_northing_center + (y_full(1:end-1)' - pi) * (domain_length_m / (2*pi));
    
    xib_utm_cell{c} = xib_utm_c; yib_utm_cell{c} = yib_utm_c;
    xib = [xib; xib_utm_c]; yib = [yib; yib_utm_c];
    f2  = [f2; ones(length(xib_utm_c), 1) * target_elevations_m(c)];
end

%% 5. IB OPERATORS (PHYSICAL SCALING)
fprintf('Assembling Physical IB Matrices...\n');
a_s = 1.2;
r_cut = 4 * a_s * h_utm;
const_factor = 1 / (2 * pi * (a_s*h_utm)^2) * h_utm^2;
X_vec = X_utm(:); Y_vec = Y_utm(:);

row_S = []; col_S = []; val_S = [];
for k = 1:length(xib)
    idx = find(abs(X_vec - xib(k)) <= r_cut & abs(Y_vec - yib(k)) <= r_cut);
    if ~isempty(idx)
        dx = xib(k) - X_vec(idx); dy = yib(k) - Y_vec(idx);
        dist_sq = (dx.^2 + dy.^2) / (a_s*h_utm)^2;
        vals = const_factor * exp(-0.5 * dist_sq);
        row_S = [row_S; idx]; col_S = [col_S; repmat(k, length(idx), 1)]; val_S = [val_S; vals];
    end
end
Sib = sparse(row_S, col_S, val_S, N^2, length(xib));
Jib = sparse(col_S, row_S, val_S, length(xib), N^2);

base_elevation_m = min(target_elevations_m) - 200; 
RHS_PDE = zeros(N^2, 1);
for j = 1:N
    for i = 1:N
        if i == 1 || j == 1 || i == N || j == N
            RHS_PDE(N*(j-1)+i) = base_elevation_m; 
        end
    end
end

% MODIFICATION 1: Soften the Immersed Boundary Constraints to prevent ringing
reg_param = 5.0; 
RHSVec = [RHS_PDE; f2];

%% 6. MINIMAL SURFACE PDE SOLVER (STABILIZED PICARD ITERATION)
fprintf('Solving Stabilized Minimal Surface PDE...\n');

% MODIFICATION 2: Increase Tikhonov regularization for stronger baseline smoothness
tikhonov_alpha = 1.5; 
max_iters = 25;
tol = 0.5; 

% INITIAL GUESS
row_L = zeros(N^2 * 5, 1); col_L = zeros(N^2 * 5, 1); val_L = zeros(N^2 * 5, 1);
idx_L = 1;
for j = 1:N
    for i = 1:N
        g = N*(j-1) + i;
        if i == 1 || j == 1 || i == N || j == N
            row_L(idx_L) = g; col_L(idx_L) = g; val_L(idx_L) = 1; idx_L = idx_L + 1;
        else
            row_L(idx_L:idx_L+4) = g; col_L(idx_L:idx_L+4) = [g, g+1, g-1, g+N, g-N];
            val_L(idx_L:idx_L+4) = [-4, 1, 1, 1, 1] / (h_utm^2); idx_L = idx_L + 5;
        end
    end
end
Lap = sparse(row_L(1:idx_L-1), col_L(1:idx_L-1), val_L(1:idx_L-1), N^2, N^2);
ResVec = [Lap, -Sib; Jib, -reg_param * speye(length(xib))] \ RHSVec; 
u = reshape(ResVec(1:N^2), N, N);

tic;
for iter = 1:max_iters
    [ux, uy] = gradient(u, h_utm); 
    D = 1 ./ sqrt(1 + ux.^2 + uy.^2); 
    
    idx_L = 1;
    for j = 1:N
        for i = 1:N
            g = N*(j-1) + i;
            if i == 1 || j == 1 || i == N || j == N
                row_L(idx_L) = g; col_L(idx_L) = g; val_L(idx_L) = 1; idx_L = idx_L + 1;
            else
                Dx_r = (D(j,i) + D(j,i+1))/2 + tikhonov_alpha; 
                Dx_l = (D(j,i) + D(j,i-1))/2 + tikhonov_alpha;
                Dy_t = (D(j,i) + D(j+1,i))/2 + tikhonov_alpha; 
                Dy_b = (D(j,i) + D(j-1,i))/2 + tikhonov_alpha;
                
                row_L(idx_L:idx_L+4) = g;
                col_L(idx_L:idx_L+4) = [g, g+1, g-1, g+N, g-N];
                val_L(idx_L:idx_L+4) = [-(Dx_r+Dx_l+Dy_t+Dy_b), Dx_r, Dx_l, Dy_t, Dy_b] / (h_utm^2);
                idx_L = idx_L + 5;
            end
        end
    end
    L_D = sparse(row_L(1:idx_L-1), col_L(1:idx_L-1), val_L(1:idx_L-1), N^2, N^2);
    
    LHSMat = [L_D, -Sib; Jib, -reg_param * speye(length(xib))];
    ResVec = LHSMat \ RHSVec;
    u_new = reshape(ResVec(1:N^2), N, N);
    
    diff = max(abs(u_new(:) - u(:)));
    fprintf('  Iter %d: Max diff = %.2f m\n', iter, diff);
    
    % MODIFICATION 3: Heavy under-relaxation to suppress high-frequency oscillations
    u = 0.2 * u_new + 0.8 * u; 
    
    if diff < tol, fprintf('  Converged!\n'); break; end
end
toc;

%% 7. VISUALIZATION
window_size_for_export = [0, 50, 1920, 975];
font = 24; lw = 4.5; lwlrg = 6; lwsml = 3;
figure(1); set(gcf, 'Position', window_size_for_export);

subplot(1, 2, 1);
surf(X_utm, Y_utm, u, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
hold on; colormap('plasma');

for c = 1:num_contours
    plot3(xib_utm_cell{c}, yib_utm_cell{c}, ones(size(xib_utm_cell{c})) * target_elevations_m(c), 'r-', 'LineWidth', lwlrg/2);
end

xlabel('UTM Easting [m]', 'Interpreter', 'latex'); 
ylabel('UTM Northing [m]', 'Interpreter', 'latex');
zlabel('Elevation [m]', 'Interpreter', 'latex');
title(sprintf('Stabilized Minimal Surface ($N=%d$)', N), 'Interpreter', 'latex');
fontsize(font, "points");
set(gca, 'TickLabelInterpreter', 'latex');
ax1 = gca; ax1.LineWidth = lwsml*.5;
view(-37.5, 30); 

subplot(1, 2, 2); hold on;
for c = 1:num_contours
    x_raw_utm = utm_easting_center + (x_scaled{c} - pi) * (domain_length_m / (2*pi));
    y_raw_utm = utm_northing_center + (y_scaled{c} - pi) * (domain_length_m / (2*pi));
    
    p1 = plot(x_raw_utm, y_raw_utm, 'k--', 'LineWidth', lwsml*.75); 
    p2 = plot(xib_utm_cell{c}, yib_utm_cell{c}, 'r-', 'LineWidth', lw*.75);
end

plot(utm_easting_center, utm_northing_center, 'bx', 'MarkerSize', 15, 'LineWidth', 3);
grid on; axis equal;
xlabel('UTM Easting [m]', 'Interpreter', 'latex'); 
ylabel('UTM Northing [m]', 'Interpreter', 'latex');
title('USGS Contours vs. Fourier Epicycles', 'Interpreter', 'latex');
legend([p1, p2], {'Raw USGS Contours', sprintf('Fourier (%d Harmonics)', num_harmonics)}, 'Interpreter', 'latex', 'Location', 'northeast');
fontsize(font, "points");
ax2 = gca; ax2.LineWidth = lwsml*.5;

%% 8. HELPER FUNCTIONS
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
    a(end+1:max_k) = 0; b(end+1:max_k) = 0; c(end+1:max_k) = 0; d(end+1:max_k) = 0;
    
    for k = 1:max_k
        x_out = x_out + a(k)*cos(k*theta) + b(k)*sin(k*theta);
        y_out = y_out + c(k)*cos(k*theta) + d(k)*sin(k*theta);
    end
end