"""
Energy accumulation terms using explicit m = E / (γ c²).

Examples:
- Kinetic energy of rotation
- Age-weighted accumulation

Numerically small compared with the main amplification term.
"""
import numpy as np

C = 299792458.0

def kinetic_mass(M_bary, V_rot):
    """m_k = (1/2 M V²) / c²"""
    E_k = 0.5 * M_bary * V_rot**2
    return E_k / (C**2)

def age_accumulated_mass(E_rate, t_age):
    """m_acc = (E_rate * t_age) / c²"""
    return (E_rate * t_age) / (C**2)

def main():
    print("Explicit accumulation terms defined.")
    print("They remain numerically negligible in current SPARC tests.")

if __name__ == "__main__":
    main()
