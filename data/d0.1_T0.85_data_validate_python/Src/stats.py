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

    def dat_to_dataframe(self, dat_file, hasExtraInfo = False, needBoxSize = False):
        if dat_file.split(".")[1] == "csv":
            df = pd.read_csv(dat_file, dtype=np.float64, header=None)
            return df
        count = 0
        with open(dat_file, 'r') as input_file:
            lines = input_file.readlines()
            newLines = []
            ##### second line last element should be stored somehow
            for line in lines:
                count += 1
                newLine = line.strip(' ').split()
                if hasExtraInfo:
                    if count == 2 and needBoxSize:
                       return newLine[-1]
                    if(count > 8):
                        newLines.append(newLine)
                else:
                    newLines.append(newLine)

        csv_file = dat_file.replace('dat', 'csv')
        with open(csv_file, 'w') as output_file:
            file_writer = csv.writer(output_file)
            file_writer.writerows(newLines)

        df = pd.read_csv(csv_file, dtype=np.float64, header=None)
        return df

