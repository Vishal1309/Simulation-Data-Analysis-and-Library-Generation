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
        self.mdfs = None
        self.cmdfs = None
        pass


    def gen_data(self):
        mdfs = []
        cmdfs = []
        for i in range(self.first_iteration, self.last_iteration + 1):
            mFilePath = f"../R{i}/moments.dat"
            cmFilePath = f"../R{i}/central_moments.dat"
            mdfs.append(self.statHelp.dat_to_dataframe(mFilePath))
            cmdfs.append(self.statHelp.dat_to_dataframe(cmFilePath))    
        self.mdfs = mdfs
        self.cmdfs = cmdfs

    def generate_properties(self, skewness_and_kurtosis = True, compressibility = True):
        self.gen_data()
        if skewness_and_kurtosis == True:
            self.gen_skewness_and_kurtosis()
        if compressibility:
            self.gen_compressibility()
    
    def gen_skewness_and_kurtosis(self):
        skewness_list = list()
        kurtosis_list = list()
        for i in range(0, self.last_iteration - self.first_iteration + 1):
            second_order_cm = self.cmdfs[i].iloc[2,0]
            std_dev = np.sqrt(second_order_cm)
            skewness = self.cmdfs[i].iloc[3,0] / np.power(std_dev, 3)
            kurtosis = self.cmdfs[i].iloc[4,0] / np.power(std_dev, 4)
            skewness_list.append(skewness)
            kurtosis_list.append(kurtosis)
        avg_skewness = np.average(skewness_list)
        atd_skewness = np.average(skewness_list)
        avg_kurtosis = np.average(kurtosis_list)
        std_kurtosis = np.average(kurtosis_list)
        print("skewness_list: ", skewness_list)
        print("kurtosis_list: ", kurtosis_list)

    def gen_compressibility(self):
        beta = 1/0.85
        box = 10
        compressibility_list = list()
        for i in range(0, self.last_iteration - self.first_iteration + 1):
            compressibility = beta * (np.power(box, 3)) * (self.mdfs[i].iloc[2,0] - np.power(self.mdfs[i].iloc[1, 0], 2)) / (np.power(self.mdfs[i].iloc[1, 0], 2))
            compressibility_list.append(compressibility)
        avg_compressiblity = np.average(compressibility_list)
        std_compressibility = np.average(compressibility_list)
        print("compressibility_list: ", compressibility_list)
