%% Symbolic Verification of Traction-Free Boundary
clear; syms zeta R m Gamma GammaP real

% 1. Define Symbolic Potentials
omega = R * (zeta + m/zeta);
omega_dot = diff(omega, zeta);

% φ(ζ) for Traction Free
phi = Gamma * R * zeta + GammaP * R / zeta;
phi_dot = diff(phi, zeta);

% ψ(ζ) derived via Muskhelishvili's BC (the "Interaction")
% On the boundary (zeta=sigma), conj(omega) = R(1/sigma + m*sigma)
% This leads to the standard interaction term:
interaction = (zeta * (1 + m*zeta^2) / (zeta^2 - m)) * phi_dot;
psi = conj(GammaP) * R * zeta - (Gamma + m * conj(GammaP)) * R / zeta - interaction;
disp(psi)
%% 2. Calculate Stresses Symbolically
% Phi(z) = phi'(zeta) / omega'(zeta)
Phi = phi_dot / omega_dot;
disp(Phi)

% Psi(z) = psi'(zeta) / omega'(zeta)
Psi = diff(psi, zeta) / omega_dot;

% Stress formulas
% sigma_sum = 4 * Re(Phi)
sigma_sum = 4 * real(Phi)
disp(sigma_sum)
% sigma_diff = 2 * ( conj(omega)/omega' * Phi' + Psi )

%% 3. Numerical Verification (The "Map")
R_v = 1.5; m_v = 0.333; % Example values
G_v = 25; GP_v = -50;   % P1=0, P2=100 -> Gamma=25, GammaP=-50

% We test at zeta = 1 (the boundary)
% Note: Use a limit if the expression is singular, or a value like 1.0001
test_zeta = 1.0001; 

S_sum = subs(4 * Phi, {R, m, Gamma, GammaP, zeta}, {R_v, m_v, G_v, GP_v, test_zeta});

% Note: real() and imag() in MATLAB don't always play nice with symbolic 
% objects until they are fully numerical.
S_sum_val = double(S_sum);

fprintf('Sum of Stresses at boundary: %.4f\n', real(S_sum_val));