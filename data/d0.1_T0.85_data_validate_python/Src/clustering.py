import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt
from stats import statHelp
from fileSaver import fileSaver

class clustering:
    def __init__(self):
        self.statHelp = statHelp()
        self.fileSaver = fileSaver()
        self.L = 10
        pass
    def assign_cluster(self, process_no, file_start, file_end, r_cut):  
        frequencies = []
        for i in range(file_start, file_end + 1):
            filePath = f"XYZ_{process_no}_{i}"
            df = self.statHelp.dat_to_dataframe(filePath)
            x = list(df[:][1])
            y = list(df[:][2])
            z = list(df[:][3])
            frequencies = self.cluster(x, y, z, r_cut, process_no, i, frequencies)
        plt.bar(np.arange(len(frequencies) + 1), frequencies)
        plt.savefig(f'freq_{process_no}_{file_start}_{file_end}.png')

    def cluster(self, x, y, z, r_cut, process_no, file_no, frequencies):
        n = len(x)
        cluster_id = 0
        cluster_assignments = {}
        if len(frequencies) == 0:
            frequencies = np.zeros(n)
        assigned = np.full(n, False, dtype=bool)
        for i in range(n):
            if assigned[i]:
                continue
            currPoint = [x[i], y[i], z[i]]
            assigned[i] = True
            cluster_id += 1
            cluster_assignments[cluster_id] = list([i])
            for j in range(i + 1, n):
                nextPoint = [x[j], y[j], z[j]]
                dis = self.distance(currPoint, nextPoint)
                if dis <= r_cut:
                    cluster_assignments[cluster_id].append(j)
                    assigned[j] = True
        for point_list in cluster_assignments.values():
            currsize = len(point_list)
            frequencies[currsize] += 1
        mat = []
        mat2 = []
        for c_id in range(1, cluster_id + 1):
            for point in cluster_assignments[c_id]:
                row = [c_id, x[point], y[point], z[point]]
                row2 = [c_id, point + 1]
                mat.append(row)
                mat2.append(row2)
        # df = pd.DataFrame(mat, columns = ['pid', 'X', 'Y', 'Z'])
        # df = pd.DataFrame(mat2, columns = ['pid', 'cid'])
        self.fileSaver.save_file('clustered', 'c_XYZ', '.dat', [process_no, file_no], mat)
        self.fileSaver.save_file('clustered', 'c_pid', '.dat', [process_no, file_no], mat2)
        return frequencies

    def distance(self, p1, p2):
        a = np.abs(p1[0] - p2[0])
        b = np.abs(p1[1] - p2[1])
        c = np.abs(p1[2] - p2[2])
        disx, disy, disz = a, b, c
        if a > self.L / 2:
            disx = self.L - a
        if b > self.L / 2:
            disy = self.L - b
        if c > self.L / 2:
            disz = self.L - c
        dis = disx**2 + disy**2 + disz**2
        dis = np.sqrt(dis)
        return dis