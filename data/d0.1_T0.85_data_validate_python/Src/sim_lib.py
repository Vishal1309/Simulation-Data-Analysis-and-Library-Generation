import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt
from lnpi_gen import lnpiGenerator
from moments_gen import momentsCalculator
from prop_gen import propGenerator
from clustering import clustering
from radial_distribution_function import radial_distribution_function
# from error_propagation import Complex, arrays_to_complex

class sim_lib:
    def __init__(self):
        pass

    def generate_lnpi(self):
        generator = lnpiGenerator()
        generator.generate_lnpi()

    def generate_moments(self, n = 1, first_iteration = 1, last_iteration = 1, x_start = 0, x_end = 100, x_id = 0, y_id = 1, strategy1 = True, strategy2 = True, isPix = False, data_file_name = "gn_m.dat"):
        calculator = momentsCalculator()
        calculator.generate_moments(n, first_iteration, last_iteration, x_start, x_end, x_id, y_id, strategy1, strategy2, isPix, data_file_name)
    
    def generate_properties(self, first_iteration, last_iteration, skewness_and_kurtosis = True, compressibility=True):
        generator_prop = propGenerator(first_iteration, last_iteration)
        generator_prop.generate_properties(skewness_and_kurtosis, compressibility)

    def generate_clustering(self, process_no, file_start, file_end, r_cut):
        clustering_generator = clustering()
        clustering_generator.assign_cluster(process_no, file_start, file_end, r_cut)

    def get_rdf(self, molecule_type, particle_id, rho, nhis, filename, box):
        rdf_generator = radial_distribution_function()
        rdf_generator.calculate(molecule_type, particle_id, rho, nhis, filename, box)

    # def run_functions(self, generate_lnpi = False, generate_moments = True, generate_properties = True, n = 1, first_iteration = 1, last_iteration = 1, x_start = 0, x_end = 100, x_id = 0, y_id = 1, strategy1 = True, strategy2 = True, isPix = False, data_file_name = "gn_m.dat", property_params = []):
    #     if generate_lnpi:
    #         generator = lnpiGenerator()
    #         generator.generate_lnpi()
    #     if generate_moments:
    #         calculator = momentsCalculator()
    #         calculator.generate_moments(n, first_iteration, last_iteration, x_start, x_end, x_id, y_id, strategy1, strategy2, isPix, data_file_name)
    #     if generate_properties:
    #         generator_prop = propGenerator(first_iteration, last_iteration)
    #         generator_prop.generate_properties(True, True)

