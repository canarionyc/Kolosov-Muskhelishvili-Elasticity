%% 1. Setup Parameters and Symbolic Variables
clear; clc;
syms z zeta s R m Gamma GammaP kappa G positive
syms x y real

% Geometry (2:1 Ellipse)
a = 2; b = 1;
R_val = (a + b)/2;
m_val = (a - b)/(a + b);

% Loading (Uniaxial Tension 100MPa along Y)
% P1 = 0 (x), P2 = 100 (y)
P1 = 0; P2 = 100; Lambda = 0;
Gamma_val = (P1 + P2) / 4;
GammaP_val = -(P1 - P2) / 2 * exp(-2i * Lambda);

%% 2. Define Potentials Symbolically
% Conformal Mapping: z = omega(zeta)
omega = R * (zeta + m/zeta);
omega_dot = diff(omega, zeta);

% Primary Potential phi(zeta)
phi = Gamma * R * zeta + GammaP * R / zeta;
phi_dot = diff(phi, zeta);

% Secondary Potential psi(zeta) derived from Traction-Free BC
% The term (zeta * (1 + m*zeta^2) / (zeta^2 - m)) is the simplification 
% of [conj(omega(1/sigma)) / omega'(sigma)]
interaction_term = (zeta * (1 + m*zeta^2) / (zeta^2 - m)) * phi_dot;
psi = conj(GammaP_val) * R * zeta - (Gamma + m * conj(GammaP_val)) * R / zeta - interaction_term;
psi_dot = diff(psi, zeta);

%% 3. Define Stress Formulas
% Phi(zeta) = phi'(zeta) / omega'(zeta)
Phi = phi_dot / omega_dot;
Phi_dot_zeta = diff(Phi, zeta);

% Psi(zeta) = psi'(zeta) / omega'(zeta)
Psi = psi_dot / omega_dot;

% Stress sum: sx + sy = 4 * Re(Phi)
% Stress diff: sy - sx + 2i*txy = 2 * [ conj(omega)/omega' * Phi' + Psi ]
% Note: In symbolic, we use 'conj' carefully or assume zeta is on boundary
sx_plus_sy = 4 * real(Phi);
complex_disp = 2 * ( (conj(omega)/omega_dot) * Phi_dot_zeta + Psi );

%% 4. Numerical Verification at Boundary
% Map back to numerical values
% subs_map = [R, m, Gamma, GammaP] -> [R_val, m_val, Gamma_val, GammaP_val];

% Let's test at the tip of the ellipse (theta = 0)
% At the boundary, rho = 1, zeta = exp(i*0) = 1
zeta_test = 1.000001; % Slightly outside to avoid singularity in omega_dot if m=1

Phi_eval = subs(Phi, [zeta, R, m, Gamma, GammaP], [zeta_test, R_val, m_val, Gamma_val, GammaP_val]);
Psi_eval = subs(Psi, [zeta, R, m, Gamma, GammaP], [zeta_test, R_val, m_val, Gamma_val, GammaP_val]);
omega_eval = subs(omega, [zeta, R], [zeta_test, R_val]);
omega_dot_eval = subs(omega_dot, [zeta, R, m], [zeta_test, R_val, m_val]);
Phi_dot_eval = subs(Phi_dot_zeta, [zeta, R, m, Gamma, GammaP], [zeta_test, R_val, m_val, Gamma_val, GammaP_val]);

S_sum = 4 * real(Phi_eval);
S_diff = 2 * ( (conj(omega_eval)/omega_dot_eval) * Phi_dot_eval + Psi_eval );

sigX = double((S_sum - real(S_diff))/2);
sigY = double((S_sum + real(S_diff))/2);
tauXY = double(imag(S_diff)/2);

fprintf('--- Results at Ellipse Tip (Theta=0) ---\n');
fprintf('Sigma X: %.4f\n', sigX);
fprintf('Sigma Y: %.4f (Stress Concentration)\n', sigY);
fprintf('Tau XY:  %.4f\n', tauXY);

%% 5. Boundary Traction Check
% Normal vector at theta=0 for this ellipse is [1, 0]
% Traction T = [sigX, tauXY; tauXY, sigY] * [1; 0] = [sigX; tauXY]
fprintf('\nTraction Tx: %.4f (Should be 0)\n', sigX);
fprintf('Traction Ty: %.4f (Should be 0)\n', tauXY);