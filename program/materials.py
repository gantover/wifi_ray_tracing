from data import *


class Material:
    def __init__(self, eps_r, sig, color):
        self.eps_r = eps_r  # relative permittivity
        self.sig = sig  # conductivity
        self.eps = eps_r * EPS0
        self.t_eps = eps_r * EPS0 -1.0j * (sig / OMEGA)
        self.Z = np.sqrt(MU0 / self.t_eps)
        self.gamma = 1.0j * OMEGA * np.sqrt(MU0 * self.t_eps)
        self.color = color


brick = Material(3.95, 0.073, "firebrick")
concrete = Material(6.4954, 1.43, "red")
division = Material(2.7, 0.05346, "forestgreen")
glass = Material(6.3919, 0.00107, "aqua")
metal = Material(1.0, 1.0 * 1e7, "silver")


"""Sample exercise data"""
# concrete = Material(4.8, 0.018)
