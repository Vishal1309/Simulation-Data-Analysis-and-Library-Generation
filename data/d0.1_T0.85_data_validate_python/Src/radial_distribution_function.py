import os
import cmath
import numpy as np
import matplotlib.pyplot as plt
from stats import statHelp
from fileSaver import fileSaver

class radial_distribution_function:
    def __init__(self):
        self.statHelp = statHelp()
        self.g = None
        self.xplt = None
        pass

    def get_relevant_XYZ(self, df, molecule_type, particle_id):
        x = list()
        y = list()
        z = list()
        for row in df.iterrows():
            if row[1] == molecule_type and row[2] == particle_id:
                x.append(row[3])
                y.append(row[4])
                z.append(row[5])
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

    def calculate(self, molecule_type, particle_id, rho, nhis, filename, box = 10): ######### WHATTTT ISSS NGR???????
        filePath = f"../{filename}" # PLEASE SPECIFY FILE NAME
        df = self.statHelp.dat_to_dataframe(filePath, hasExtraInfo=True)
        x, y, z = self.get_relevant_XYZ(df, molecule_type, particle_id)
        npart = len(x)
        dr = box / (2 * nhis)
        self.g = np.zeros(nhis)
        for i in range(1, npart - 1):
            for j in range(i + 1, npart):
                r = self.periodic_boundary_condition(box, x[i], x[j], y[i], y[j], z[i], z[j])
                if r < box / 2:
                    idx = int(r / dr)
                    self.g[idx] = self.g[idx] + 2
        self.xplt = list()
        ngr = 1 ###########################################################################3
        for i in range(1, nhis):
            r = dr * (i + 0.5)
            self.xplt.append(r)
            vb = ( np.power(i + 1, 3) - np.power(i, 3) ) * np.power(dr, 3)
            nid = (4 / 3) * (np.pi) * (vb) * (rho)
            self.g[i] = self.g[i] / (ngr * npart * nid)
        self.plot(box, nhis)
        return self.g
    
    def plot(self, box, nhis): ############### will the width be dr ???
        dr = box / (2 * nhis)
        fig = plt.figure()
        plt.bar(self.xplt, self.g, dr)
        plt.savefig("barplot.png") ############## what should be the name here?????

                
                