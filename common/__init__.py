"""
Common initialization module for the Muskhelishvili Theory of Vibrations project.
This module provides shared resources like the unit registry and common imports.
"""
import matplotlib.pyplot as plt

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pint

# Create a single unit registry for the entire project
ureg = pint.UnitRegistry()
ureg.setup_matplotlib(enable=True)

# Define Quantity as a shorthand for creating quantities
Quantity = ureg.Quantity

# Export commonly used components
__all__ = ['ureg', 'Quantity', 'np', 'sp', 'plt']
