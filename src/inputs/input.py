import json
import math
import os

class InputParameters:
    def __init__(self, file_path):
        with open(file_path, 'r') as f:
            config_json = json.load(f)

        # Physical constants
        self.hbar = config_json["physical_constants"]["hbar"]
        self.kB = config_json["physical_constants"]["kB"]

        # Open chain parameters
        self.mass = config_json["open_chain_parameters"]["mass"]
        self.P = config_json["open_chain_parameters"]["P"]
        self.P_plus_1 = self.P + 1

        # NHC parameters
        self.ExDoF = config_json["NHC_parameters"]["ExDoF"]
        self.tau = config_json["NHC_parameters"]["tau"]
        self.dN = config_json["NHC_parameters"]["dN"]

        # Langevin parameters
        self.gamma = config_json["LG_parameters"]["gamma"]

        # Simulation parameters
        self.T = config_json["simulation_parameters"]["T"]
        self.beta = 1.0 / (self.kB * self.T)
        self.omegaP = math.sqrt(self.P) / (self.beta * self.hbar)
        self.N = config_json["simulation_parameters"]["N"]
        self.x_min = config_json["simulation_parameters"]["x_min"]
        self.x_max = config_json["simulation_parameters"]["x_max"]
        self.dt = config_json["simulation_parameters"]["dt"]
        self.steps = config_json["simulation_parameters"]["steps"]
        self.save_interval = config_json["simulation_parameters"]["save_interval"]
        self.potential_type = config_json["simulation_parameters"]["potential_type"]

        # Derived result directory
        self.result_dir = f"../../results/{self.potential_type}_beta={self.beta:.2f}/"

        # Measurement parameters
        self.s_min = config_json["measurement_parameters"]["s_min"]
        self.s_max = config_json["measurement_parameters"]["s_max"]
        self.p_min = config_json["measurement_parameters"]["p_min"]
        self.p_max = config_json["measurement_parameters"]["p_max"]
        self.q_min = config_json["measurement_parameters"]["q_min"]
        self.q_max = config_json["measurement_parameters"]["q_max"]
        self.nbins = config_json["measurement_parameters"]["nbins"]

        # Save raw dictionary for optional use
        self.params = config_json

    def get_param(self, *keys):
        data = self.params
        for key in keys:
            data = data[key]
        return data