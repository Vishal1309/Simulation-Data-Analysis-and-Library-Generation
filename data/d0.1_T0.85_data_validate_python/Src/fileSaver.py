import pandas as pd
import csv
import numpy as np
import os
import cmath

class fileSaver:
    def __init__(self):
        pass

    def save_file(self, folder_name, file_name, file_format, file_identifiers, file_data):
        fileName = f"../{folder_name}/{file_name}"
        for identifier in file_identifiers:
            fileName += f"_{identifier}"
        fileName += file_format
        np.savetxt(fileName, file_data, delimiter=" ")