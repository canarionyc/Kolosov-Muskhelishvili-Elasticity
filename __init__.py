"""
Provides a single unit registry for consistent unit handling throughout the package.
"""
import matplotlib.pyplot as plt

import pint
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Create a single unit registry for the entire package
ureg = pint.UnitRegistry()
ureg.setup_matplotlib(enable=True)

# Define Quantity as a shorthand for creating quantities
Quantity = ureg.Quantity

# Export the registry and Quantity for use in other modules
__all__ = ['ureg', 'Quantity']