%% 1. Setup Symbolic Variables
clear; clc;
syms zeta R m Gamma GammaP real
syms s positive % A real parameter for evaluating near the boundary

% Constants for your specific test case:
% a=2, b=1 => R=1.5, m=0.333
% P1=0, P2=100, Lambda=0 => Gamma=25, GammaP=-50
R_val = 1.5; m_val = 1/3;
G_val = 25; GP_val = -50;

%% 2. Define Mapping and Potentials
omega = R * (zeta + m/zeta)
omega_dot = diff(omega, zeta)

% Primary Potential phi(zeta)
phi = Gamma * R * zeta + GammaP * R / zeta
phi_dot = diff(phi, zeta)
% disp(phi_dot)

% Secondary Potential psi(zeta) derived for Traction-Free Hole
% This exact term ensures the BC: phi(s) + omega(s)/conj(omega'(s)) * conj(phi'(s)) + conj(psi(s)) = 0
interaction = (zeta * (1 + m*zeta^2) / (zeta^2 - m)) * phi_dot
psi = conj(GammaP) * R * zeta - (Gamma + m * conj(GammaP)) * R / zeta - interaction
psi_dot = diff(psi, zeta)

%% 3. Define Stress Components
% Phi(z) = phi'(zeta) / omega'(zeta)
Phi = phi_dot / omega_dot
Phi_prime_zeta = diff(Phi, zeta)

% Psi(z) = psi'(zeta) / omega'(zeta)
Psi_z = psi_dot / omega_dot

%% Stress Formulas
% sigma_x + sigma_y = 4 * Re[Phi]
% sigma_y - sigma_x + 2i*tau_xy = 2 * [ conj(omega)/omega' * Phi' + Psi ]
S_sum = 4 * Phi
S_diff = 2 * ( (conj(omega)/omega_dot) * Phi_prime_zeta + Psi_z )

%% 4. Boundary Verification (at zeta = 1)
% At the boundary (zeta=1), the normal is in X direction. 
% Therefore, Tx = sigma_x and Ty = tau_xy. Both should be zero.

% We substitute numerical values and set zeta = 1.0001 to avoid the coordinate singularity
test_map = {R, m, Gamma, GammaP, zeta}
test_vals = {R_val, m_val, G_val, GP_val, 1.0001}

SUM = double(subs(S_sum, test_map, test_vals))
DIFF = double(subs(S_diff, test_map, test_vals))

sigX = (real(SUM) - real(DIFF)) / 2
sigY = (real(SUM) + real(DIFF)) / 2
tauXY = imag(DIFF) / 2
%%
fprintf('--- Verification at Ellipse Tip (Theta=0, z=a) ---\n');
fprintf('Sigma X (Traction Tx): %10.6f (Target: 0)\n', sigX);
fprintf('Sigma Y (Hoop Stress): %10.6f (Target: 500 for a=2,b=1)\n', sigY);
fprintf('Tau XY (Traction Ty):  %10.6f (Target: 0)\n', tauXY);