%% 1. Setup Parameters and Symbolic Variables
clear; clc;
syms z zeta sigma r theta real
syms R sig_inf positive

% We assume uniaxial tension at infinity: σ_xx = 0, σ_yy =  σ_inf, τ_xy = 0.

%% 2. Mapping and Derivative
% Mapping: z = ω(ζ) = R*ζ. This maps the unit circle (|ζ|=1) to a circle of radius R.

% For a circular hole, the mapping is simple scaling.
omega = R * zeta;
omega_dot = diff(omega, zeta); % ω'(ζ) = R

%% 3. Boundary Condition and Potential Form
% Boundary Condition on the hole (|ζ|=1):
% φ(σ) + ω(σ)/conj(ω'(σ)) * conj(φ'(σ)) + conj(ψ(σ)) = 0
% For a circle: ω(σ)/conj(ω'(σ)) = (R*σ)/R = σ.
% BC becomes: φ(σ) + σ * conj(φ'(σ)) + conj(ψ(σ)) = 0.



% We assume general Laurent series for exterior problem:
% φ(ζ) = A*R*ζ + B*R/ζ
% ψ(ζ) = C*R*ζ + D*R/ζ + E*R/ζ^3
syms A B C D E

% Infinity Conditions (Uniaxial Tension):
% A = σ_inf / 4
% C = σ_inf / 2 (Wait, C logic derived below)

A = sig_inf / 4;
C = sig_inf / 2;

phi = A * R * zeta + B * R / zeta;
phi_prime_zeta = diff(phi, zeta);

%% 4. Solving for Coefficients (Symbolic Match)
% On boundary ζ = σ (where σ*conj(σ) = 1)
% conj(phi_prime_zeta) becomes A*R - B*R*σ^2
% conj(psi) becomes C*R/σ + D*R*σ + E*R*σ^3
term1 = A*R*sigma + B*R/sigma;
term2 = sigma * (A*R - B*R*sigma^2);
term3 = C*R/sigma + D*R*sigma + E*R*sigma^3;

BC_eqn = term1 + term2 + term3; % Must be zero for all sigma
% Collect terms: (2A + D)*R*sigma + (B + C)*R/sigma + (E - B)*R*sigma^3 = 0

% Solve the system:
D_val = -2*A;         % From coeff sigma
B_val = -C;           % From coeff 1/sigma
E_val = B_val;        % From coeff sigma^3

%% Final Symbolic Potentials in z-domain (ζ = z/R)
phi_z = A*z + (B_val*R^2)/z;
psi_z = C*z + (D_val*R^2)/z + (E_val*R^4)/z^3;

%% 5. Stress Field Derivation
% σ_xx + σ_yy = 4 * Re[phi'(z)]
% σ_yy - σ_xx + 2i*tau_xy = 2 * [conj(z)*phi''(z) + psi'(z)]

phi_p = diff(phi_z, z);
phi_pp = diff(phi_p, z);
psi_p = diff(psi_z, z);

%% Substitute z = r*exp(i*theta)
z_polar = r * exp(1i*theta);
z_conj = r * exp(-1i*theta);

S_sum = 4 * real(subs(phi_p, z, z_polar));
S_diff = 2 * (z_conj * subs(phi_pp, z, z_polar) + subs(psi_p, z, z_polar));
%% cartesian stress in polars (r,theta)
clc
sig_xx = (S_sum - real(S_diff))/2
sig_yy = (S_sum + real(S_diff))/2
tau_xy = imag(S_diff)/2

% disp(sig_xx); disp(sig_yy); disp(tau_xy)

%% Far field verification
limit(sig_xx,r,inf)
limit(sig_yy,r,inf)

%% verification at the rim
subs(sig_xx,[r,theta],[R,0])
subs(sig_yy,[r,theta],[R,0])
subs(tau_xy,[r,theta],[R,0])

%%
subs(sig_xx,[r,theta],[R,pi/2])
subs(sig_yy,[r,theta],[R,pi/2])
subs(tau_xy,[r,theta],[R,pi/2])
%% 
fprintf('Symbolic σ_yy at (r=R, theta=0): %s\n', char(simplify(subs(sig_yy, [r, theta], [R, 0]))));
% Result should be 3*sig_inf (Stress Concentration Factor)

%% 10. Parameters and Grid Setup
clear; clc; close all;

% Physical Parameters
R = 1.0;            % Radius of hole
sig_inf = 10.0;     % Far-field stress (Tension in Y)
E_mod = 210e3;      % Young's Modulus (e.g., Steel in MPa)
nu = 0.3;           % Poisson's ratio
mu = E_mod / (2*(1+nu));        % Shear Modulus
kappa = 3 - 4*nu;               % Kolosov constant (Plane Strain)

% Create Grid (Polar for easy masking, then convert to Cartesian)
r_min = R; 
r_max = 5*R;
nr = 60; ntheta = 90;
r_vec = linspace(r_min, r_max, nr);
theta_vec = linspace(0, 2*pi, ntheta);
[r_grid, theta_grid] = meshgrid(r_vec, theta_vec);

% Convert to Complex Plane z
X = r_grid .* cos(theta_grid);
Y = r_grid .* sin(theta_grid);
z = X + 1i*Y;

%% 2. Calculate Potentials (Using your coefficients)
% Coefficients for Y-Tension (as verified previously)
A = sig_inf / 4;
C = sig_inf / 2; 
D = -2*A;
B = -C;
E = B;

% Potentials phi(z) and psi(z)
phi = A*z + (B*R^2)./z;
psi = C*z + (D*R^2)./z + (E*R^4)./z.^3;

% Derivatives needed for stress/disp
% phi'(z) = A - B*R^2/z^2
phi_prime = A - (B*R^2)./(z.^2);
% phi''(z) = 2*B*R^2/z^3
phi_double_prime = (2*B*R^2)./(z.^3);
% psi'(z) = C - D*R^2/z^2 - 3*E*R^4/z^4
psi_prime = C - (D*R^2)./(z.^2) - (3*E*R^4)./(z.^4);

%% 3. Compute Stress Field
% Formulas:
% sig_xx + sig_yy = 4 * Re[phi'(z)]
% sig_yy - sig_xx + 2i*tau_xy = 2 * [conj(z)*phi''(z) + psi'(z)]

S_sum = 4 * real(phi_prime);
S_diff_term = 2 * (conj(z).*phi_double_prime + psi_prime);

sig_xx = (S_sum - real(S_diff_term)) / 2;
sig_yy = (S_sum + real(S_diff_term)) / 2;
tau_xy = imag(S_diff_term) / 2;

% Calculate Hoop Stress (Sigma_theta) for the contour plot
% Transformation: sig_theta = sig_x*sin^2 + sig_y*cos^2 - 2*tau*sin*cos
sig_theta = sig_xx.*sin(theta_grid).^2 + ...
            sig_yy.*cos(theta_grid).^2 - ...
            2*tau_xy.*sin(theta_grid).*cos(theta_grid);

%% 4. Compute Displacement Field
% Formula: 2*mu*(u + iv) = kappa*phi(z) - z*conj(phi'(z)) - conj(psi(z))
disp_complex = (kappa*phi - z.*conj(phi_prime) - conj(psi)) / (2*mu);
U = real(disp_complex);
V = imag(disp_complex);

%% 5. Visualization
figure('Color', 'w', 'Position', [100 100 900 700]);

% --- Plot: Hoop Stress Contour + Displacement Quiver ---
% We use 20 contour levels for smooth gradients
[C_plot, h] = contourf(X, Y, sig_theta, 20, 'LineStyle', 'none'); 
hold on;
colormap(jet); 
c = colorbar;
c.Label.String = 'Hoop Stress \sigma_{\theta\theta} (MPa)';
clim([-sig_inf, 3*sig_inf]); % Set limits to highlight concentration

% Add Displacement Vectors (Quiver)
% We subsample the grid for quiver so arrows aren't too dense
sample_step = 4; 
q = quiver(X(1:sample_step:end, 1:sample_step:end), ...
           Y(1:sample_step:end, 1:sample_step:end), ...
           U(1:sample_step:end, 1:sample_step:end), ...
           V(1:sample_step:end, 1:sample_step:end), ...
           'k', 'LineWidth', 1.2);

% Draw the hole boundary circle
viscircles([0 0], R, 'Color', 'k', 'LineWidth', 2);

% Formatting
axis equal;
title(['Stress Concentration (Y-Tension) & Displacement Field', newline, ...
       'Max \sigma_{\theta\theta} \approx 3\sigma_{inf} at \theta=0, \pi']);
xlabel('X (mm)');
ylabel('Y (mm)');
grid on;
set(gca, 'FontSize', 12);
hold off;
