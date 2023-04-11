import os
import cmath
import numpy as np
import matplotlib.pyplot as plt
from stats import statHelp
from fileSaver import fileSaver

class radial_distribution_function:
    def __init__(self):
        self.statHelp = statHelp()
        self.fileSaver = fileSaver()
        self.g = None
        self.xplt = None
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
    
    def periodic_boundary_condition(self, box, xi, xj, yi, yj, zi, zj):
        dx = xi - xj
        dy = yi - yj
        dz = zi - zj
        dx = dx - box * round(dx / box)
        dy = dy - box * round(dy / box)
        dz = dz - box * round(dz / box)
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        return r
    
    def append_zeros(self, n, k):
        curr = str(n)
        curr = curr[::-1]
        while len(curr) != k:
            curr += '0'
        curr = curr[::-1]
        return curr
    
    def calculate_processed_data(self, x1, y1, z1, x2, y2, z2, molecule_type1, molecule_type2, particle_id1, particle_id2, rho, nhis, box = 10):
        npart = len(x1)
        mpart = len(x2)
        dr = box / (2 * nhis)
        g = np.zeros(nhis)
        for i in range(1, npart):
            next = 1
            if molecule_type1 == molecule_type2 and particle_id1 == particle_id2:
                next = i + 1
            for j in range(next, mpart):
                r = self.periodic_boundary_condition(box, x1[i], x2[j], y1[i], y2[j], z1[i], z2[j])
                if r < box / 2:
                    idx = int(r / dr)
                    g[idx] = g[idx] + 2
        
        self.xplt = list()
        ngr = 1 ###########################################################################3
        for i in range(1, nhis):
            r = dr * (i + 0.5)
            self.xplt.append(r)
            vb = ( np.power(i + 1, 3) - np.power(i, 3) ) * np.power(dr, 3)
            nid = (4 / 3) * (np.pi) * (vb) * (rho)
            g[i] = g[i] / (ngr * npart * nid)
        return g, self.xplt

    def calculate(self, process_no_begin, process_no_end, file_start, file_end, molecule_type1, molecule_type2, particle_id1, particle_id2, rho, nhis, box = 10): ######### WHATTTT ISSS NGR???????
        ########## any atom type and any molecules
        ########## separate file for combining clustering and rdf - molecules from 1 cluster 
        denominator = 0
        self.g = np.zeros(nhis)
        for pr in range(int(process_no_begin), int(process_no_end) + 1):
            for i in range(int(file_start), int(file_end) + 1):
                curr_pr = self.append_zeros(pr, 2)
                curr_file_id = self.append_zeros(i, 4)
                filePath = f"XYZ_{curr_pr}_{curr_file_id}"
                df = self.statHelp.dat_to_dataframe(filePath, hasExtraInfo=True)
                x1, y1, z1 = self.get_relevant_XYZ(df, molecule_type1, particle_id1)
                x2, y2, z2 = self.get_relevant_XYZ(df, molecule_type2, particle_id2)
                denominator += 1
                g, self.xplt = self.calculate_processed_data(self, x1, y1, z1, x2, y2, z2, molecule_type1, molecule_type2, particle_id1, particle_id2, rho, nhis, denominator, box = 10)
                for i in range(nhis):
                    self.g[i] += g[i]
        np.divide(self.g, denominator)
        mat = np.array([self.xplt, self.g]).T
        self.fileSaver.save_file('rdf', 'rdf', '.dat', [process_no_begin, process_no_end, file_start, file_end, molecule_type1, molecule_type2, particle_id1, particle_id2, rho, nhis, box], mat)
        # self.plot(box, nhis)
        return self.g

    
    # def plot(self, box, nhis): ############### will the width be dr ???
    #     dr = box / (2 * nhis)
    #     fig = plt.figure()
    #     plt.bar(self.xplt, self.g, dr)
    #     plt.savefig("barplot.png") ############## what should be the name here?????
    #     ############# print in a file

                
                