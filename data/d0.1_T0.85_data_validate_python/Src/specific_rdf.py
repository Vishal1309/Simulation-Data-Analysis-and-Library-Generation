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

    def separate_points_across_clusters(self, cluster_id_list):
        x1, y1, z1 = list(), list(), list()
        x2, y2, z2 = list(), list(), list()
        for key in self.cluster_assignments.keys():
            if key in cluster_id_list:
                for coordinate in self.cluster_assignments[key]:
                    x1.append(coordinate[0])
                    y1.append(coordinate[1])
                    z1.append(coordinate[2])
            else:
                for coordinate in self.cluster_assignments[key]:
                    x2.append(coordinate[0])
                    y2.append(coordinate[1])
                    z2.append(coordinate[2])
        return x1, y1, z1, x2, y2, z2
                


    def get_secific_rdf(self, process_no_begin, process_no_end, file_start, file_end, molecule_type, particle_id, r_cut, rho, nhis, box = 10): 
        ############################# PROCESS BEGIN and PROCESS END MUST BE THE SAME, FILE START AND FILE END MUST ALSO BE THE SAME
        cluster_generator = clustering()
        cluster_generator.assign_cluster(process_no_begin, process_no_end, file_start, file_end, molecule_type, particle_id, r_cut)
        self.cluster_assignments = cluster_generator.get_cluster_assignments()
        x1, y1, z1, x2, y2, z2 = self.separate_points_across_clusters([1])
        rdf_calculator = radial_distribution_function()
        g = rdf_calculator.calculate_processed_data(x1, y1, z1, x2, y2, z2, molecule_type, molecule_type, particle_id, particle_id, rho, nhis, box)
        mat = np.array([self.xplt, self.g]).T
        self.fileSaver.save_file('rdf', 'rdf', '.dat', [process_no_begin, process_no_end, file_start, file_end, molecule_type, molecule_type, particle_id, particle_id, rho, nhis, box], mat)
        return g

