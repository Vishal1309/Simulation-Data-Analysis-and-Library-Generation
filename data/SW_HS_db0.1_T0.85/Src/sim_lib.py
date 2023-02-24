import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt
# from lnpi_gen import lnpiGenerator
# from moments_gen import momentsCalculator
from prop_gen import propGenerator
# from error_propagation import Complex, arrays_to_complex

class sim_lib:
    def __init__(self):
        return

    def run_functions(self, generate_lnpi = False, generate_moments = False, generate_properties = False, n = 1, first_iteration = 1, last_iteration = 1, x_start = 0, x_end = 100, x_id = 0, y_id = 1, strategy1 = True, strategy2 = True, isPix = False, data_file_name = "gn_m.dat", property_params = []):
        print("ENTERED RUN FUNCTIONS")
        print("booleans", generate_lnpi, generate_moments, generate_properties)
        # if generate_lnpi == True:
        #     generator = lnpiGenerator()
        #     generator.generate_lnpi()
        # if generate_moments == True:
        #     calculator = momentsCalculator()
        #     calculator.generate_moments(n, first_iteration, last_iteration, x_start, x_end, x_id, y_id, strategy1, strategy2, isPix, data_file_name)
        # if generate_properties == True:
        #     generator_prop = propGenerator(first_iteration, last_iteration)
        #     generator_prop.generate_properties()

