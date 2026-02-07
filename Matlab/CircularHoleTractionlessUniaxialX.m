%% Master Visualization Script: Infinite Plate with Circular Hole
% Load Case: Uniaxial Tension in X-Direction (Sigma_xx = Sig_Inf)
clear; clc; close all;
export_flag = true; 

%% 1. Parameters & Grid Generation
R = 1.0;                  % Hole Radius
Sig_Inf = 100;            % Far-field Stress (MPa)
E = 210e3;                % Young's Modulus (MPa)
nu = 0.3;                 % Poisson's Ratio

% Grid Setup (Polar is best for this geometry, then convert to Cartesian)
r_vec = linspace(R, 5*R, 80);
th_vec = linspace(0, 2*pi, 180);
[r, th] = meshgrid(r_vec, th_vec);

x = r .* cos(th);
y = r .* sin(th);

%% 2. Analytical Solution (Kirsch Equations for X-Tension)
% Using standard Kirsch solution for Tension in X to match your request 
% (Stress concentration at +/- 90 degrees)

% Stress Components (Polar)
s_r = (Sig_Inf/2) * (1 - R^2./r.^2) + ...
      (Sig_Inf/2) * (1 - 4*(R^2./r.^2) + 3*(R^4./r.^4)) .* cos(2*th);

s_th = (Sig_Inf/2) * (1 + R^2./r.^2) - ...
       (Sig_Inf/2) * (1 + 3*(R^4./r.^4)) .* cos(2*th);

t_rth = -(Sig_Inf/2) * (1 + 2*(R^2./r.^2) - 3*(R^4./r.^4)) .* sin(2*th);

% Convert Stresses to Cartesian for Contours
s_x = s_r.*cos(th).^2 + s_th.*sin(th).^2 - 2*t_rth.*sin(th).*cos(th);
s_y = s_r.*sin(th).^2 + s_th.*cos(th).^2 + 2*t_rth.*sin(th).*cos(th);
t_xy = (s_r - s_th).*sin(th).*cos(th) + t_rth.*(cos(th).^2 - sin(th).^2);

% Von Mises Stress (Plane Stress)
s_vm = sqrt(s_x.^2 + s_y.^2 - s_x.*s_y + 3*t_xy.^2);

% Principal Stress Difference (For Photoelasticity)
% Sigma1 - Sigma2 = sqrt((sx-sy)^2 + 4*txy^2)
diff_principal = sqrt((s_x - s_y).^2 + 4*t_xy.^2);

%% 3. Displacement Calculation
% Plane Stress conditions
G = E / (2*(1+nu));
k = (3-nu) / (1+nu); % Kolosov constant for Plane Stress

% Displacement (Radial and Tangential)
u_r = (Sig_Inf*R / (8*G)) * ( (r/R).*(k-1) + (2*R./r) + ...
      2*((r/R).*(k+1) + R./r - R^3./r.^3).*cos(2*th) );
      
v_th = -(Sig_Inf*R / (8*G)) * ( 2*((r/R).*(k+1) + R./r + R^3./r.^3).*sin(2*th) );

% Convert Displacements to Cartesian
u_x = u_r.*cos(th) - v_th.*sin(th);
u_y = u_r.*sin(th) + v_th.*cos(th);

%% FIGURE 1: Stress Field Dashboard (Contours)
figure('Color','w', 'Name', 'Stress Field Analysis', 'Position', [50 50 1000 800]);
t = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Stress Field Dashboard (Uniaxial Tension X-Direction)', 'FontSize', 16);

% --- Subplot 1: Sigma_xx ---
nexttile; 
contourf(x, y, s_x, 20, 'LineStyle','none');
axis equal; 
colormap(gca, parula); 
colorbar; 
title('\sigma_{xx} (Longitudinal)', 'Color','k');
viscircles([0 0], R, 'Color','k', 'LineWidth',1); 
grid on;

% --- Subplot 2: Sigma_yy ---
nexttile; 
contourf(x, y, s_y, 20, 'LineStyle','none');
axis equal; 
colormap(gca, parula); 
colorbar; 
title('\sigma_{yy} (Transverse)');
viscircles([0 0], R, 'Color','k', 'LineWidth',1); 
grid on;

% --- Subplot 3: Shear Tau_xy ---
nexttile; 
contourf(x, y, t_xy, 20, 'LineStyle','none');
axis equal; 
colormap(gca, jet); 
colorbar; 
title('\tau_{xy} (Shear)');
viscircles([0 0], R, 'Color','k', 'LineWidth',1); 
grid on;

% --- Subplot 4: Von Mises ---
nexttile; 
contourf(x, y, s_vm, 20, 'LineStyle','none');
axis equal; 
colormap(gca, hot); 
colorbar; 
title('Von Mises Equivalent Stress');
viscircles([0 0], R, 'Color','k', 'LineWidth',1); 
grid on;

% One-off export with white background (ignores screen colors)
if(export_flag)
     exportgraphics(gcf, 'KirschX1.png', 'BackgroundColor', 'w');
end
exportgraphics(gcf, 'KirschX1.png', 'BackgroundColor', 'w');

%% FIGURE 2: Engineering Analysis (Polar & Radial)
figure('Color','w', 'Name', 'Engineering Validation', 'Position', [100 100 900 500]);
t2 = tiledlayout(1, 2);

% --- Left: Polar Plot of Hoop Stress at Boundary ---
nexttile;
% Extract boundary values (r = R is the first index of r_vec)
s_th_rim = s_th(:,1); 
polarplot(th_vec, s_th_rim, 'LineWidth', 2, 'Color', 'r');
title({'Boundary Hoop Stress \sigma_{\theta\theta}', '(r=R)'});
pax = gca; 
pax.ThetaZeroLocation = 'top'; % Orientation
text(pi/2, max(s_th_rim), '  \leftarrow K_t = 3', 'Color', 'r', 'FontSize',12);

% --- Right: Radial Decay at Key Angles ---
nexttile;
hold on;
% Find indices for specific angles: 0 (0 deg), 90 (pi/2)
[~, idx_0] = min(abs(th_vec - 0));
[~, idx_90] = min(abs(th_vec - pi/2));

plot(r_vec/R, s_th(idx_90, :)/Sig_Inf, 'r-', 'LineWidth', 2, 'DisplayName', '\theta=90^\circ (Max Stress)');
plot(r_vec/R, s_th(idx_0, :)/Sig_Inf, 'b--', 'LineWidth', 1.5, 'DisplayName', '\theta=0^\circ (Compressive)');
yline(1, 'k:', 'DisplayName', 'Far Field \sigma_{\infty}');

xlabel('Radial Distance (r/R)');
ylabel('Normalized Stress (\sigma_{\theta\theta} / \sigma_{\infty})');
title('Stress Decay Analysis');
legend; grid on; box on;
xlim([1 5]);

% One-off export with white background (ignores screen colors)
if(export_flag)
     exportgraphics(gcf, 'KirschX2.png', 'BackgroundColor', 'w');
end

%% FIGURE 3: Displacement & Deformation (Deformed Mesh)
figure('Color','w', 'Name', 'Deformation Visualization', 'Position', [150 150 800 600]);

% Create a coarser grid for mesh visualization so lines are visible
r_c = linspace(R, 4*R, 15);
th_c = linspace(0, 2*pi, 30);
[rc, thc] = meshgrid(r_c, th_c);
xc = rc .* cos(thc);
yc = rc .* sin(thc);

% Re-calc displacement for coarse grid
ur_c = (Sig_Inf*R / (8*G)) * ( (rc/R).*(k-1) + (2*R./rc) + ...
       2*((rc/R).*(k+1) + R./rc - R^3./rc.^3).*cos(2*thc) );
vth_c = -(Sig_Inf*R / (8*G)) * ( 2*((rc/R).*(k+1) + R./rc + R^3./rc.^3).*sin(2*thc) );
ux_c = ur_c.*cos(thc) - vth_c.*sin(thc);
uy_c = ur_c.*sin(thc) + vth_c.*cos(thc);

% Scale Factor for Visibility
scale = 0.15 * R / max(max(abs(ux_c))); % Scale max disp to 15% of Radius

% Plot Original Mesh (Gray dashed)
% Plotting radial lines
plot(xc, yc, 'Color', [0.7 0.7 0.7], 'LineStyle', ':'); hold on;
plot(xc', yc', 'Color', [0.7 0.7 0.7], 'LineStyle', ':');

% Plot Deformed Mesh (Blue solid)
X_def = xc + ux_c * scale;
Y_def = yc + uy_c * scale;
plot(X_def, Y_def, 'b', 'LineWidth', 1.2);
plot(X_def', Y_def', 'b', 'LineWidth', 1.2);

viscircles([0 0], R, 'Color', 'k', 'LineStyle', '--');
title(sprintf('Deformed Mesh (Scale Factor: %.0f)', scale));
axis equal; axis off;

% One-off export with white background (ignores screen colors)
if(export_flag)
     exportgraphics(gcf, 'KirschX3.png', 'BackgroundColor', 'w');
end         


%% FIGURE 4: Photoelasticity (Isochromatics)
figure('Color','k', 'Name', 'Simulated Photoelasticity', 'Position', [200 200 700 600]);

% Isochromatics depend on the Principal Stress Difference
% We use 'hsv' or 'jet' to simulate interference fringes
contourf(x, y, diff_principal, 30, 'LineStyle', 'none');
colormap(gca, 'jet'); 
colorbar('Color', 'w');
viscircles([0 0], R, 'Color','w', 'LineWidth', 2);
axis equal; axis off;
title('Simulated Isochromatics (\sigma_1 - \sigma_2)', 'Color', 'w');

% One-off export with white background (ignores screen colors)
if(export_flag)
     exportgraphics(gcf, 'KirschX4.png', 'BackgroundColor', 'w');
end   
