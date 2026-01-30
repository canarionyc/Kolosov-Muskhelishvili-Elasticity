# %% setup
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import os

import numpy as np
import matplotlib.pyplot as plt


def map_triangle():
	theta = np.linspace(0, 2 * np.pi, 500)
	rho = 1.0  # Boundary of unit circle
	zeta = rho * np.exp(1j * theta)

	# Parameters for Triangle
	# z = b * (zeta + m * zeta^n)
	n = -2  # Wait, for a triangle we often use 1/zeta terms for EXTERIOR maps.
	# But for the EPITROCHOID (Interior map), we use positive powers.

	# Let's stick to the text's Epitrochoid definition (Eq 134.1a):
	# z = b(zeta + m * zeta^n)

	# Try n=4 (aiming for 3 corners? or n=3 for 2?)
	# Actually, visual inspection:
	# 1 + m*e^(i*3*theta) -> This has period 2pi/3 (3 lobes).
	# This comes from n=4 in the expansion? No, terms are usually 1 and n.
	# Actually, simply: z = zeta + m * zeta^(-2) ??

	# Let's use the PROVEN mapping for a curvilinear triangle often used in Savin:
	# z = A * (zeta + m * zeta^(-2))  <-- This maps EXTERIOR of circle to EXTERIOR of triangle
	# But Muskhelishvili §134a uses POSITIVE powers for bounded regions (Epitrochoid).

	# Let's try the Epitrochoid formula specifically:
	# z = rho * e^(it) + m * rho^n * e^(int)

	# Set n = -2 (Hypotrochoid) vs n = 4 (Epitrochoid)

	# TRIANGLE APPROXIMATION (Epitrochoid variant)
	# We want 3 lobes. We need cos(3*theta).
	# That comes from interference between exp(i*theta) and exp(-i*2*theta).

	m_triangle = 1 / 3.5  # Approx 0.28

	# Mapping for bounded region (Epitrochoid / curvilinear polygon):
	# z = zeta + (m * zeta^4) ? No, that gives 3 lobes but they point "in".
	# Let's try:
	z_tri = zeta + (1.0 / 4.0) * zeta ** (-2)  # Commonly cited for triangle exterior.

	# Let's stick to the BOOK's Epitrochoid (Section 134a):
	# z = b * (zeta + m * zeta^n)
	# This actually creates shapes with n-1 loops if m is large, or n-1 "corners".
	# Let's plot n=4, m=0.2

	z_epitrochoid = zeta + 0.22 * zeta ** 4

	plt.figure(figsize=(6, 6))
	plt.plot(z_epitrochoid.real, z_epitrochoid.imag, 'b-', linewidth=2)
	plt.title("Epitrochoid with n=4 (Approx Triangle?)")
	plt.axis('equal')
	plt.grid(True)
	plt.show()


map_triangle()