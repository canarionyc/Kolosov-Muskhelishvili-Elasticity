%% 1. Setup Parameters and Symbolic Variables
clear; clc;
syms z % the complex variable  in the physical plane z=x+iy 
syms zeta % the variable in the transformed plane, where z = ω(ζ) 
syms sig_inf theta real
syms R r positive

% We assume uniaxial tension at infinity: σ_xx = 0, σ_yy =  σ_inf, τ_xy = 0.

%% 2. Mapping and Derivative
% Mapping: z = ω(ζ) = R*ζ. This maps the unit circle (|ζ|=1) to a circle of radius R.

% For a circular hole, the mapping is simple scaling.
omega = R * zeta;
omega_dot = diff(omega, zeta); % ω'(ζ) = R

%% 3. Potential Form
% We assume general Laurent series for exterior problem:
% φ(ζ) = A*ζ + B/ζ
% ψ(ζ) = C*ζ + D/ζ + E/ζ^3
syms A B C D E real

% Infinity Conditions (Uniaxial Tension):
% A = σ_inf / 4
% C = σ_inf / 2 (Wait, C logic derived below)

% A = sig_inf / 4;
% C = sig_inf / 2;

Phi = A * zeta + B / zeta
Phi_prime_zeta = diff(Phi, zeta)

Psi = C*zeta + D/zeta + E/zeta^3

%% 3. Boundary Condition   
% Boundary Condition on the hole (|ζ|=1):
% φ(σ) + ω(σ)/conj(ω'(σ)) * conj(φ'(σ)) + conj(ψ(σ)) = 0
% For a circle: ω(σ)/conj(ω'(σ)) = (R*σ)/R = σ.
% BC becomes: φ(σ) + σ * conj(φ'(σ)) + conj(ψ(σ)) = 0.

%% 4. Solving for Coefficients (Symbolic Match)
% On boundary ζ = σ (where σ*conj(σ) = 1)
% conj(phi_prime_zeta) becomes A*R - B*R*σ^2
% conj(psi) becomes C*R/σ + D*R*σ + E*R*σ^3
syms sigma
term1 = A*sigma + B/sigma;
term2 = sigma * (A - B*sigma^2);
term3 = C/sigma + D*sigma + E*sigma^3;

BC_eqn = term1 + term2 + term3; % Must be zero for all sigma
%%
% simplify(BC_eqn)
expr=expand(BC_eqn*sigma)
expr=collect(expr,sigma)
eqns=coeffs(expr, sigma)

[Bsol, Dsol, Esol] = solve(eqns,[B,D,E])
%% Collect terms: (2A + D)*sigma + (B + C)/sigma + (E - B)*sigma^3 = 0
% collect(BC_eqn, sigma)
% coeffs(BC_eqn, sigma)
%% Solve the system:
% D_val = -2*A;         % From coeff sigma
% B_val = -C;           % From coeff 1/sigma
% E_val = B_val;        % From coeff sigma^3

%% Final Symbolic Potentials in z-domain (ζ = z/R)

% phi_z = A*z + (Bsol*R^2)/z
% psi_z = C*z + (Dsol*R^2)/z + (Esol*R^4)/z^3
Phi=subs(Phi,[B],[Bsol])
Psi=subs(Psi,[D,E],[Dsol,Esol])

Phi_polar = @(r,theta) Phi(r .* exp(1i*theta))
Psi_polar = @(r,theta) Psi(r .* exp(1i*theta))

%%
phi=subs(Phi,zeta, z/R)
psi=subs(Psi,zeta, z/R)
%% 5. Stress Field Derivation
% σ_xx + σ_yy = 4 * Re[phi'(z)]
% σ_yy - σ_xx + 2i*tau_xy = 2 * [conj(z)*phi''(z) + psi'(z)]

phi_p = diff(phi, z)
phi_pp = diff(phi_p, z)
psi_p = diff(psi, z)

%% Substitute z = r*exp(i*theta)
z_polar = r * exp(1i*theta);

S_sum = 4 * real(subs(phi_p, z, z_polar))
S_diff = 2 * (conj(z_polar) * subs(phi_pp, z, z_polar) + ...
	subs(psi_p, z, z_polar))
%% cartesian stress in polars (r,theta)
% clc
sig_xx = (S_sum - real(S_diff))/2
sig_yy = (S_sum + real(S_diff))/2
tau_xy = imag(S_diff)/2

% disp(sig_xx); disp(sig_yy); disp(tau_xy)

%% Far field verification
limit(sig_xx,r,inf)
limit(sig_yy,r,inf)

far_eqns=[limit(sig_xx,r,inf)==0, 
	limit(sig_yy,r,inf)==sig_inf]

[Asol,Csol]=solve(far_eqns,[A,C],"Real",true)
%%
phi=subs(phi,[A,C],[Asol,Csol])
psi=subs(psi,[A,C],[Asol,Csol])

phi_p = diff(phi, z);
phi_pp = diff(phi_p, z);
psi_p = diff(psi, z);

S_sum = 4 * real(phi_p)
S_diff = 2 * (conj(z) * phi_pp + psi_p)

syms sig_xx sig_yy x y real
% eqns=[
% 	sig_xx+sig_yy==subs(S_sum,z, x+i*y)
% 	sig_yy-sig_xx==real(subs(S_diff,z, x+i*y))
% 	]
% soln=solve(eqns, [sig_xx,sig_yy])

%% polar components
% R = so2(theta,"theta") # requires  Navigation Toolbox

% syms sig_xx sig_yy tau_xy real
sigma = [sig_xx, tau_xy;
         tau_xy, sig_yy]

theta                 % scalar angle
R = [cos(theta), -sin(theta);
     sin(theta),  cos(theta)]

sigma_rot = R * sigma * R.'   % rotated stress tensor
sig_xx_p = simplify(sigma_rot(1,1))
sig_yy_p = simplify(sigma_rot(2,2))
sig_xy_p = simplify(sigma_rot(1,2))


%% verification at the rim
subs(sig_xx,[r,theta],[R,0])
subs(sig_yy,[r,theta],[R,0])
subs(tau_xy,[r,theta],[R,0])

%% verification at the rim
subs(sig_xx,[r,theta],[R,0])
subs(sig_yy,[r,theta],[R,0])
subs(tau_xy,[r,theta],[R,0])

%%
clc
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
Phi = A*z + (B*R^2)./z;
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
disp_complex = (kappa*Phi - z.*conj(phi_prime) - conj(psi)) / (2*mu);
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
