import numpy as np
import cea
from matplotlib import pyplot as plt

class injector_sizing:
    def __init__(self, ox_mass, fuel_mass, ox_cd, fuel_cd):
        self.ox_mass = ox_mass # kg/s
        self.fuel_mass = fuel_mass # kg/s
        self.ox_cd = ox_cd
        self.fuel_cd = fuel_cd
