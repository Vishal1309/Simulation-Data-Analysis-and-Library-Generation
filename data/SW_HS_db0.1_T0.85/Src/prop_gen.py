import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt
from stats import statHelp

# from error_propagation import Complex, arrays_to_complex

class propGenerator:
    def __init__(self, first_iteration, last_iteration):
        self.statHelp = statHelp()
        self.first_iteration = first_iteration
        self.last_iteration = last_iteration
        pass


    def gen_data(self):
        for i in range(self.first_iteration, self.last_iteration + 1):
            mFilePath = f"../R{i}/moment.dat"
            cmFilePath = f"../R{i}/central_moment.dat"
            mdf = self.statHelp.dat_to_dataframe(mFilePath)
            cmdf = self.statHelp.dat_to_dataframe(cmFilePath)
            print(mdf.head())
            print(cmdf.head())

    def generate_properties(self, skewness = True, kurtosis = True, compressibility = True):
        self.gen_data()
        if skewness:
            self.gen_skewness()
        if kurtosis:
            self.gen_kurtosis()
        if compressibility:
            self.gen_compressibility()

