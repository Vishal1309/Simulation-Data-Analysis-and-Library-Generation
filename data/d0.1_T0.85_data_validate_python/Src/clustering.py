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
        self.point_assignments = None
        self.L = 10
        pass

    def get_cm(self, x2d, y2d, z2d):
        pass

    def get_relevant_XYZ(self, df, molecule_type, particle_id):
        x = list()
        y = list()
        z = list()
        for row in df.iterrows():
            if row[1][1] == molecule_type and row[1][2] == particle_id:
                x.append(row[1][3])
                y.append(row[1][4])
                z.append(row[1][5])
        return x, y, z
    
    def append_zeros(self, n, k):
        curr = str(n)
        curr = curr[::-1]
        while len(curr) != k:
            curr += '0'
        curr = curr[::-1]
        return curr
    
    def get_cluster_assignments(self):
        return self.point_assignments

    def assign_cluster(self, process_no_begin, process_no_end, file_start, file_end, molecule_type, particle_id, r_cut):  
        frequencies = []   
        for pr in range(int(process_no_begin), int(process_no_end) + 1):
            for i in range(int(file_start), int(file_end) + 1):
                curr_pr = self.append_zeros(pr, 2)
                curr_file_id = self.append_zeros(i, 4)
                filePath = f"../config_d0.7_T0.45/xyz_{curr_pr}_{curr_file_id}.dat"
                df = self.statHelp.dat_to_dataframe(filePath, hasExtraInfo=True)
                print(df.head())
                x, y, z = self.get_relevant_XYZ(df, molecule_type, particle_id)
                ####### center of mass - MODIFY LATER
                # x, y, z = self.get_cm(x2d, y2d, z2d)
                frequencies = self.cluster(x, y, z, r_cut, pr, i, frequencies)
        mat = np.array([np.arange(1, len(frequencies) + 1), frequencies]).T
        self.fileSaver.save_file('clusters/cluster_frequencies', 'c', '.dat', [process_no_begin, process_no_end, file_start, file_end], mat)
        ########### print in a dat file N vs freq

    def cluster(self, x, y, z, r_cut, process_no, file_no, frequencies):
        n = len(x)
        last_cluster_id = 0
        curr_cluster_id = 0
        cluster_assignments = {}
        point_assignments = {}
        assigned_cluster = {}
        if len(frequencies) == 0:
            frequencies = np.zeros(n)
        assigned = np.full(n, False, dtype=bool)
        for i in range(n):
            currPoint = [x[i], y[i], z[i]]
            if assigned[i]:
                curr_cluster_id = assigned_cluster[i]
            else:
                assigned[i] = True
                last_cluster_id += 1
                curr_cluster_id = last_cluster_id
                cluster_assignments[curr_cluster_id] = list([i])
                point_assignments[curr_cluster_id] = list([x[i], y[i], z[i]])
            for j in range(i + 1, n):
                nextPoint = [x[j], y[j], z[j]]
                dis = self.distance(currPoint, nextPoint)
                if dis <= r_cut:
                    cluster_assignments[curr_cluster_id].append(j)
                    point_assignments[curr_cluster_id].append([x[j], y[j], z[j]])
                    assigned[j] = True
        for point_list in cluster_assignments.values():
            currsize = len(point_list)
            frequencies[currsize] += 1
        mat = []
        mat2 = []
        for c_id in range(1, last_cluster_id + 1):
            for point in cluster_assignments[c_id]:
                row = [c_id, x[point], y[point], z[point]]
                row2 = [c_id, point + 1]
                mat.append(row)
                mat2.append(row2)
        self.point_assignments = point_assignments
        # df = pd.DataFrame(mat, columns = ['pid', 'X', 'Y', 'Z'])
        # df = pd.DataFrame(mat2, columns = ['pid', 'cid'])
        self.fileSaver.save_file('clusters/cluster_assignments', 'c_XYZ', '.dat', [process_no, file_no], mat)
        self.fileSaver.save_file('clusters/cluster_assignments', 'c_pid', '.dat', [process_no, file_no], mat2)
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