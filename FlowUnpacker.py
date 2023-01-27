from netCDF4 import Dataset as CDFData
from math import atan, sqrt
import matplotlib.pyplot as plt
import numpy as np

class Database:

    def __init__(self, filepath):
        rootgroup = CDFData(filepath, "r", format="NETCDF4")
        self.ERRX, self.ERRY = rootgroup.variables["ERRX"][()], rootgroup.variables["ERRY"][()]
        self.lat, self.lon = rootgroup.variables["lat"][()], rootgroup.variables["lon"][()]
        self.VX, self.VY = rootgroup.variables["VX"][()], rootgroup.variables["VY"][()]
        self.x, self.y = rootgroup.variables["x"][()], rootgroup.variables["y"][()]

    def compute_angle(self):
        self.angle = np.arctan(self.VY / self.VX)

    def compute_speed(self):
        self.speed = sqrt(self.VX**2 + self.VY**2)

    def compute_error(self):
        self.error = sqrt(self.ERRX**2 + self.ERRY**2)

    def compute_angle_error(self):
        try:
            self.angle_error = self.error / (2*self.speed)
        except:
            raise Exception(f"Both error and speed need to be computed first to calculate angle error")

    def compute_all(self):
        self.compute_angle()
        self.compute_speed()
        self.compute_error()
        self.compute_angle_error()

    def close(self):
        self.rootgroup.close()
