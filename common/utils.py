"""
Utility functions for the Muskhelishvili Theory of Vibrations project.
This module provides reusable functions for calculations and visualizations.
"""
from typing import Tuple, Union, Callable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from . import ureg, Quantity

#%% Complex Variable Methods

def complex_potential_phi(z: complex, a: float, sigma_0: float, omega: float, t: float) -> complex:
    """
    Calculate the phi complex potential for a crack of length 2a.

    Args:
        z: Complex coordinate
        a: Half-length of the crack
        sigma_0: Amplitude of the applied load
        omega: Circular frequency of vibration
        t: Time

    Returns:
        Complex value of phi(z)
    """
    return sigma_0 * np.exp(1j * omega * t) * z / np.sqrt(z**2 - a**2)

def complex_potential_psi(z: complex, a: float, sigma_0: float, omega: float, t: float) -> complex:
    """
    Calculate the psi complex potential for a crack of length 2a.

    Args:
        z: Complex coordinate
        a: Half-length of the crack
        sigma_0: Amplitude of the applied load
        omega: Circular frequency of vibration
        t: Time

    Returns:
        Complex value of psi(z)
    """
    return -sigma_0 * np.exp(1j * omega * t) * a**2 / (z**2 - a**2)**(3/2)

#%% Visualization Functions

def plot_complex_function(func: Callable,
                         x_range: Tuple[float, float],
                         y_range: Tuple[float, float],
                         nx: int = 100,
                         ny: int = 100,
                         cmap: str = 'viridis',
                         title: str = '',
                         figsize: Tuple[float, float] = (10, 8)) -> Tuple[Figure, Axes]:
    """
    Plot a complex function as a heatmap.

    Args:
        func: Complex function to plot
        x_range: Range of x values as (min, max)
        y_range: Range of y values as (min, max)
        nx: Number of points in x direction
        ny: Number of points in y direction
        cmap: Colormap to use
        title: Plot title
        figsize: Figure size

    Returns:
        Matplotlib figure and axes
    """
    x = np.linspace(x_range[0], x_range[1], nx)
    y = np.linspace(y_range[0], y_range[1], ny)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j*Y

    # Evaluate function and get magnitude
    values = func(Z)
    magnitude = np.abs(values)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.pcolormesh(X, Y, magnitude, cmap=cmap, shading='auto')
    fig.colorbar(im, ax=ax, label='Magnitude')

    ax.set_title(title)
    ax.set_xlabel('Re(z)')
    ax.set_ylabel('Im(z)')
    ax.set_aspect('equal')

    return fig, ax
