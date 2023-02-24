import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt

class statHelp:
    def __init__(self):
        pass

    def custom_power(self, a, n):
        if(n == 0):
            return 1
        if(n == 1):
            return a
        if(n % 2 == 0):
            return self.custom_power(a, n // 2) * self.custom_power(a, n // 2)
        else:
            return a * self.custom_power(a, n // 2) * self.custom_power(a, n // 2)

    def normalize(self, Pix):
        sum = np.sum(Pix)
        norm_Pix = np.divide(Pix, sum)
        return norm_Pix

    def exponent(self, lnPix):
        lnPix = np.array(lnPix)
        mxlnPix = np.max(lnPix)
        lnPix = np.subtract(lnPix, mxlnPix)
        Pix = np.exp(lnPix)
        return Pix

    def dat_to_dataframe(self, dat_file):
        if dat_file.split(".")[1] == "csv":
            df = pd.read_csv(dat_file, dtype=np.float64, header=None)
            return df
        with open(dat_file, 'r') as input_file:
            lines = input_file.readlines()
            newLines = []
            for line in lines:
                newLine = line.strip(' ').split()
                newLines.append(newLine)

        csv_file = dat_file.replace('dat', 'csv')
        with open(csv_file, 'w') as output_file:
            file_writer = csv.writer(output_file)
            file_writer.writerows(newLines)

        df = pd.read_csv(csv_file, dtype=np.float64, header=None)
        return df

#
def dat_to_dataframe(self, dat_file):
        with open(dat_file, 'r') as input_file:
            lines = input_file.readlines()
            newLines = []
            count = 0
            for line in lines:
                count+=1
                if(count >= 3):
                    newLine = line.strip(' ').split()
                    newLines.append(newLine)

        csv_file = dat_file.replace('XYZ', 'csv')
        with open(csv_file, 'w') as output_file:
            file_writer = csv.writer(output_file)
            file_writer.writerows(newLines)

        df = pd.read_csv(csv_file, dtype=np.float64, header=None)
        return df