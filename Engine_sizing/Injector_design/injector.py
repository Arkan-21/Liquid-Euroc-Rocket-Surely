import numpy as np
import cea

class Injector_Sizing:
    def __init__(self, fuel, ox, Pc, MR, expansion_ratio, injection_velocity, injector_diameter, m_dot, m_dot_film):
        self.fuel = fuel
        self.ox = ox
        self.Pc = Pc
        self.MR = MR
        self.expansion_ratio = expansion_ratio
        self.injection_velocity = injection_velocity
        self.injector_diameter = injector_diameter
        self.m_dot = m_dot
        self.m_dot_film = m_dot_film