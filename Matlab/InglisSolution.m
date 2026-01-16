%% Parameters
a = 2.0; b = 1.0;
R = (a + b)/2; m = (a - b)/(a + b);
P1 = 0; P2 = 100; Lambda = 0;

%% Constants
Gamma = (P1 + P2) / 4;
GammaP = -(P1 - P2) / 2 * exp(-2i * Lambda);

%% Evaluation Point (Far Field)
z = 100 + 0i;
zeta = (z + sqrt(z^2 - 4*R^2*m))/(2*R);

%% Potentials
dPhi_dZeta = Gamma*R - GammaP*R/zeta^2;
omegaP = R*(1 - m/zeta^2);
Phi = dPhi_dZeta / omegaP;

%% Stress Traces
sigma_sum = 4 * real(Phi);
fprintf('Far Field sigmaX + sigmaY: %.4f (Expected: 100)\n', sigma_sum);