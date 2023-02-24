import pandas as pd
import csv
import numpy as np
import os
import cmath
import matplotlib.pyplot as plt
from stats import statHelp
from fileSaver import fileSaver
# from error_propagation import Complex, arrays_to_complex


class momentsCalculator:
    def __init__(self):
        self.statHelp = statHelp()
        self.fileSaver = fileSaver()
        self.first_iteration = None
        self.last_iteration = None
        self.x_start = None
        self.x_end = None
        self.isPix = None
        self.x_id = None
        self.y_id = None
        self.n = None
        self.data_file_name = None
        self.moments_for_all_steps = {}
        self.avg_ith_order_moment = None
        self.avg_ith_order_central_moment = None
        self.std_itd_order_moment = None
        self.std_itd_order_central_moment = None
        pass

    def generate_moments(self, n, first_iteration, last_iteration, x_start, x_end, x_id = 0, y_id = 1, strategy1 = True, strategy2 = True, isPix = False, data_file_name = "gn_m.dat"):
        self.n = n 
        self.first_iteration = first_iteration
        self.last_iteration = last_iteration
        self.x_start = x_start
        self.x_end = x_end
        self.x_id = x_id
        self.y_id = y_id
        self.isPix = isPix
        self.data_file_name = data_file_name
        if strategy1:
            self.compute_moments()
        if strategy2:
            self.compute_moments_with_error()

    def ith_order_moment(self, order, norm_Pix, x, first_order_moment=None):
        moment = 0
        central_moment = 0
        sum_norm_Pix = 0
        begin_id = np.where(x == self.x_start)[0][0]
        end_id = np.where(x == self.x_end)[0][0]
        print(f"entered {order}th order moment calculation")
        print(begin_id, end_id)
        for xi in range(begin_id, end_id + 1):
            moment += self.statHelp.custom_power(xi, order) * norm_Pix[xi]
            if(first_order_moment):
                central_moment += self.statHelp.custom_power(xi - first_order_moment, order) * norm_Pix[xi]
            sum_norm_Pix += norm_Pix[xi]
        print(moment, sum_norm_Pix)
        moment /= sum_norm_Pix
        if(first_order_moment):
            central_moment /= sum_norm_Pix
        print(moment, central_moment)
        print("---------------------------------------------------------------")
        return moment, central_moment

    
    def get_complex_y(self, df):
        N = df.shape[0]
        all_y = np.zeros((N, self.last_iteration - self.first_iteration + 1))
        for i in range(self.first_iteration, self.last_iteration + 1):
            filePath = self.get_step_file_path(i)
            df = self.statHelp.dat_to_dataframe(filePath)
            norm_Pix = self.gen_norm_pix(df)
            for j in range(0, N):
                all_y[j][i - self.first_iteration] = norm_Pix[j]
        avg_y = np.zeros(N)
        std_y = np.zeros(N)
        for i in range(0, N):
            avg_y[i] = np.average(all_y[i])
            std_y[i] = np.std(all_y[i])
        y_with_error = [complex(avg_y[i], std_y[i]) for i in range(0, N)]
        return y_with_error
    
    def compute_moments_with_error(self):
        filePath = self.get_step_file_path(1)
        df = self.statHelp.dat_to_dataframe(filePath)
        y_with_error = self.get_complex_y(df)
        all_moments_with_error = {
            "moments": [],
            "central_moments": []
        }
        all_moments_with_error["moments"].append(complex(1, 0))
        all_moments_with_error["central_moments"].append(complex(1, 0))
        first_order_moment, _ = self.ith_order_moment(
            1, y_with_error, df[self.x_id])
        for order in range(1, self.n + 1):
            curr_order_moment, curr_order_central_moment = self.ith_order_moment(order, y_with_error, df[self.x_id], first_order_moment)
            all_moments_with_error["moments"].append(curr_order_moment)
            all_moments_with_error["central_moments"].append(curr_order_central_moment)
        split_moments = [(v.real, v.imag) for v in all_moments_with_error["moments"]]
        split_central_moments = [(v.real, v.imag) for v in all_moments_with_error["central_moments"]]
        self.fileSaver.save_file("avg_std_2", "complex_m", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], split_moments)
        self.fileSaver.save_file("avg_std_2", "complex_cm", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], split_central_moments)
        
        # VALIDATION USING MANUAL COMPUTATION
        return all_moments_with_error

    def get_step_file_path(self, i):
        return f"../R{i}/{self.data_file_name}"
    
    def gen_norm_pix(self, df): 
        if not self.isPix:
            lnPix = df[self.y_id]
            Pix = self.statHelp.exponent(lnPix)
        else:
            Pix = df[self.y_id]
        norm_Pix = self.statHelp.normalize(Pix)
        return norm_Pix

    def gen_avg_std(self, mat_moments, mat_central_moments):
        avg_ith_order_moment = [1] * (self.n + 1)
        avg_ith_order_central_moment = [1] * (self.n + 1)
        std_ith_order_moment = [0] * (self.n + 1)
        std_ith_order_central_moment = [0] * (self.n + 1)
        for order in range(1, self.n + 1):
            avg_ith_order_moment[order] = np.average(mat_moments[order])
            avg_ith_order_central_moment[order] = np.average(mat_central_moments[order])
            std_ith_order_moment[order] = np.std(mat_moments[order])
            std_ith_order_central_moment[order] = np.std(mat_central_moments[order])
        self.fileSaver.save_file("avg_std_1", "avg_m", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], avg_ith_order_moment)
        self.fileSaver.save_file("avg_std_1", "avg_cm", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], avg_ith_order_moment)
        self.fileSaver.save_file("avg_std_1", "std_m", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], std_ith_order_moment)
        self.fileSaver.save_file("avg_std_1", "std_cm", ".dat", [self.first_iteration, self.last_iteration, self.x_start, self.x_end], std_ith_order_central_moment)
        self.avg_ith_order_moment = avg_ith_order_moment
        self.avg_ith_order_central_moment = avg_ith_order_central_moment
        self.std_itd_order_moment = std_ith_order_moment
        self.avg_ith_order_central_moment = std_ith_order_central_moment

    def compute_moments(self):
        moments_for_all_steps = {}
        mat_moments = np.zeros((self.n + 1, self.last_iteration - self.first_iteration + 1))
        mat_central_moments = np.zeros((self.n + 1, self.last_iteration - self.first_iteration + 1))
        for i in range(self.first_iteration, self.last_iteration + 1):
            file_path = self.get_step_file_path(i)
            df = self.statHelp.dat_to_dataframe(file_path)
            norm_Pix = self.gen_norm_pix(df)
            all_moments = {
                "moments": [],
                "central_moments": []
            }
            all_moments["moments"].append(1)
            all_moments["central_moments"].append(1)
            first_order_moment, _ = self.ith_order_moment(
                1, norm_Pix, df[self.x_id])
            for order in range(1, self.n + 1):
                curr_order_moment, curr_order_central_moment = self.ith_order_moment(order, norm_Pix, df[self.x_id], first_order_moment)
                mat_moments[order][i - self.first_iteration] = curr_order_moment
                mat_central_moments[order][i - self.first_iteration] = curr_order_central_moment
                all_moments["moments"].append(curr_order_moment)
                all_moments["central_moments"].append(curr_order_central_moment)
            moments_for_all_steps[i] = all_moments
            np.savetxt(f"../R{i}/moments.dat", all_moments["moments"], delimiter=" ")
            np.savetxt(f"../R{i}/central_moments.dat", all_moments["central_moments"], delimiter=" ")
        self.moments_for_all_steps = moments_for_all_steps
        self.gen_avg_std(mat_moments, mat_central_moments)
        
        # GET THE X, LNPIX FILE SEPARATELY FROM GLOBAL FILE AND THEN PERFORM OPERATION ON THE MIDFILE
        # CHECKING THE CORRECTNESS OF OUR ALGO
        # MAKE A SEPARATE FUNCTION TO PRINT/ SAVE EVERYTHING
        
        return moments_for_all_steps

    

    def plot(self, first_iteration, last_iteration, x_id, y_id, data_file_name):
        self.data_file_name = data_file_name
        for i in range(first_iteration, last_iteration + 1):
            filePath = self.get_step_file_path(i)
            df = self.statHelp.dat_to_dataframe(filePath)
            plt.plot(df[x_id], df[y_id])
            plt.savefig('../Plots/plot_R' + str(i) + '.png')
            plt.clf()

#### TEST1
# obj = momentsCalculator()
# obj.generate_moments(4, 1, 1, 345, 750, 0, 2, True, True, False, "gnn.dat")




# {
#     3:  {
#         'moments': [1, 9.961757286839747, 111.4177089220305, 1371.0508940213213, 18318.79145394204],
#         'central_moments': [1, 4.601108376401581e-16, 12.181100680125729, 18.444383360351804, 483.0664603267238]
#     },
#     4:  {
#         'moments': [1, 9.962653609694904, 111.43502216456254, 1371.3223375733196, 18322.856184798075],
#         'central_moments': [1, 5.009173166243004e-16, 12.180555217795602, 18.432506941791065, 483.013157086012]
#     },
#     5:  {
#         'moments': [1, 9.963358313591712, 111.45089118289772, 1371.6310842459445, 18328.808303324822],
#         'central_moments': [1, -9.037931130659682e-15, 12.182382297880714, 18.451041127731948, 483.2725587454138]
#     },
#     6: {
#         'moments': [1, 9.961261809100714, 111.4041695010152, 1370.775873482867, 18313.768105787593],
#         'central_moments': [1, -7.992677086823902e-15, 12.177432671566923, 18.444584449473346, 482.93602356252404]
#     }
# }
