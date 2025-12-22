#%% setup
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import os


def plot_epitrochoid(b_val=1, m_val=0.3, n_val=3):
    """
    Plots the epitrochoid shape for given parameters.
    """
    # 1. Create a range of angles for the unit circle in the zeta-plane
    theta = np.linspace(0, 2 * np.pi, 400)
    zeta_vals = np.exp(1j * theta)

    # 2. Define the mapping function using numpy
    # z = b * (zeta + m * zeta**n)
    z_vals = b_val * (zeta_vals + m_val * zeta_vals**n_val)

    # 3. Get the real and imaginary parts for plotting
    x_vals = np.real(z_vals)
    y_vals = np.imag(z_vals)

    # 4. Create the plot
    plt.figure(figsize=(8, 8))
    title = f'Epitrochoid (n={n_val}, m={m_val:.2f})'
    if m_val == 1/n_val:
        title = 'Triangle with Sharp Corners (Deltoid)'
    elif m_val < 1/n_val:
        title = 'Rounded Triangle'

    plt.plot(x_vals, y_vals, label=f'n={n_val}, m={m_val:.2f}')
    plt.title(title)
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True, linestyle='--')
    plt.axis('equal')
    plt.legend()
    plt.show()


def calc_rigidity_epitrochoid():
	# 1. Setup Conformal Map Variables
	zeta = sp.symbols('zeta')  # Complex variable in circle plane
	sigma = sp.symbols('sigma')  # Boundary variable (on unit circle)
	n = sp.symbols('n', integer=True)  # n=3 is triangle-like
	b, m = sp.symbols('b m', real=True)

	# 2. Define Mapping Function (Eq 134.1a in text)
	# z = omega(zeta) = b * (zeta + m * zeta**n)
	omega = b * (zeta + m * zeta ** n)

	# On the boundary, sigma_bar = 1/sigma
	omega_sigma = omega.subs(zeta, sigma)
	omega_bar_sigma = b * (1 / sigma + m * (1 / sigma) ** n)
	d_omega_sigma = sp.diff(omega_sigma, sigma)

	# 3. Integrate using Residue Theorem
	# The integral is a closed contour loop around sigma=0
	integrand = omega_sigma ** 2 * omega_bar_sigma * d_omega_sigma

	# We find the Residue at sigma=0 (coeff of 1/sigma)
	# This replaces the integration step
	residue = sp.residue(integrand, sigma, 0)

	# The integral result is 2*pi*i * Residue
	integral_val = 2 * sp.pi * sp.I * residue

	# 4. Calculate D (Simplified term I + Integral)
	# (Note: calculating Polar Moment of Inertia 'I' via mapping is also possible)
	print(f"Integral Term Value: {sp.simplify(integral_val)}")

if __name__ == '__main__':
	print(os.getcwd())
	calc_rigidity_epitrochoid()
	
	# Plot a triangle with sharp corners (a deltoid)
	plot_epitrochoid(b_val=1, m_val=1/3, n_val=3)

	# Plot a triangle with rounded corners
	plot_epitrochoid(b_val=1, m_val=0.2, n_val=3)
