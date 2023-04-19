import os
import cmath
import numpy as np
import matplotlib.pyplot as plt
from stats import statHelp
from fileSaver import fileSaver
from clustering import clustering
from radial_distribution_function import radial_distribution_function

class specific_rdf:
    def __init__(self):
        self.cluster_assignments = None
        self.fileSaver = fileSaver()
        pass

    def separate_points_across_clusters(self, file_idx, cluster_id, molecule_in_cluster, molecule_outside_cluster, inSameCluster = False):
        x1, y1, z1 = list(), list(), list()
        x2, y2, z2 = list(), list(), list()
        for key in self.cluster_assignments[file_idx].keys():
            if key == cluster_id or inSameCluster:
                for coordinate in self.cluster_assignments[file_idx][key]:
                    if coordinate[3] == molecule_in_cluster:
                        x1.append(coordinate[0])
                        y1.append(coordinate[1])
                        z1.append(coordinate[2])
                    if inSameCluster and coordinate[3] == molecule_outside_cluster:
                        x2.append(coordinate[0])
                        y2.append(coordinate[1])
                        z2.append(coordinate[2])
            else:
                for coordinate in self.cluster_assignments[file_idx][key]:
                    if coordinate[3] == molecule_outside_cluster:
                        x2.append(coordinate[0])
                        y2.append(coordinate[1])
                        z2.append(coordinate[2])
        return x1, y1, z1, x2, y2, z2
                
    def get_secific_rdf(self, process_no_begin, process_no_end, file_start, file_end, molecule_in_cluster, molecule_outside_cluster, particle_id1, particle_id2, r_cut, rho, nhis, box = 10): 
        cluster_generator = clustering()
        cluster_generator.assign_cluster(process_no_begin, process_no_end, file_start, file_end, r_cut)
        self.cluster_assignments = cluster_generator.get_cluster_assignments()
        denominator = 0
        g_avg = np.zeros(nhis)
        for pr in range(int(process_no_begin), int(process_no_end) + 1):
            for i in range(int(file_start), int(file_end) + 1):
                rdf_calculator = radial_distribution_function()
                denominator += 1
                temp_g = np.zeros(nhis)
                n_clusters = 0
                for cluster_id in self.cluster_assignments[denominator]:
                    n_clusters += 1
                    x1, y1, z1, x2, y2, z2 = self.separate_points_across_clusters(denominator, cluster_id, molecule_in_cluster, molecule_outside_cluster, True)
                    g, xplt = rdf_calculator.calculate_processed_data(x1, y1, z1, x2, y2, z2, molecule_in_cluster, molecule_outside_cluster, particle_id1, particle_id2, rho, nhis, box)
                    for i in range(nhis):
                        temp_g += g[i]
                temp_g = np.divide(temp_g, n_clusters)
                for i in range(nhis):
                    g_avg += temp_g[i]
        g_avg = np.divide(g, denominator)
        mat = np.array([xplt, g]).T
        self.fileSaver.save_file('rdf', 'rdf', '.dat', [process_no_begin, process_no_end, file_start, file_end, molecule_in_cluster, molecule_outside_cluster, particle_id1, particle_id2, rho, nhis, box], mat)
        return g

