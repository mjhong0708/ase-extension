from ase_extension import _ext

from .bias import BiasPotential


class LogFermiSphericalWallPotential(BiasPotential):
    """Apply logfermi potential for confined molecular dynamics.
    Confines the system to be inside a sphere by applying wall potential.

    Method referenced from https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#confining-in-a-cavity
    """

    def __init__(self, radius=5.0, temperature=300, beta=6):
        self.radius = radius
        self.temperature = temperature
        self.beta = beta

    def _get_bias_energy_and_force(self, atoms):
        E, E_grad = _ext.log_fermi_spherical_potential(atoms.positions, self.radius, self.temperature, self.beta)
        return E, -E_grad
