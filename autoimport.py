"""
Auto-import module for Muskhelishvili Theory of Vibrations project.

Usage:
    Add `%run autoimport.py` at the top of each notebook.

This will import all common resources without cluttering your notebooks with import statements.
"""
# Import common module contents into global namespace
from common import ureg, Quantity, np, sp, plt
from common.utils import complex_potential_phi, complex_potential_psi, plot_complex_function

# Print a short confirmation message
print("✅ Common resources loaded: unit registry, NumPy, SymPy, Matplotlib")
print("✅ Utility functions loaded: complex_potential_phi, complex_potential_psi, plot_complex_function")
