from MainCCA import main
from FunctionsCCA import *
from FractalDimensionCCA import fractalDimension
import datetime
import time
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def mainScript(lat_size, particles, total):
    Cluster.classReset()
    field_names = ["date_and_time", "elapsed_time_sec", "lattice_size", "num_particles", "fractal_dim", "percolates",
                   "z_1", "z_2", "z_3", "z_4", "z_5", "z_6"]
    # Create results directory to store a .csv containing results from all clusters.
    if not os.path.isdir('Results'):
        os.mkdir('Results')

    # Creates .csv file for specific cluster lattice size and amount of particles, will store date run, elapsed time,
    # lattice size, number of particles, fractal dimension and if the cluster percolates the system.
    if not os.path.isfile('Results/ResultsOffLSinSolap.csv'):
        with open('Results/ResultsOffLSinSolap.csv', mode='w') as file:
            writer = csv.DictWriter(file, fieldnames=field_names)
            writer.writeheader()

    # If a previous file existed this will count the number of lines to only add the amount of required rows.
    df = pd.read_csv('Results/ResultsOffLSinSolap.csv'.format(lat_size,particles))
    try:
        num_lines = len(df[df['num_particles'] == particles])
    except:
        num_lines = 0

    # Runs the main loop and fractal dimension code n times and writes it to .csv to keep results.
    for i in range(num_lines, total):
        with open('Results/ResultsOffLSinSolap.csv', mode='+a', newline='') as file:
            dict_writer = csv.DictWriter(file, fieldnames=field_names)

            date = datetime.datetime.now() # Gets current date and time
            start = time.time() # Starts internal clock
            x, y, z, percolation = main(lat_size,particles) # Runs main code and creates a cluster utilizing CCA
            *_, z_count = np.unique(np.array(z), return_counts = True)
            for _ in range(4 - len(z_count)):
                z_count = np.append(z_count, 0)
            print(z_count)

            center_x, center_y = centerOfMass(lat_size,particles,x,y) # Calculates center of mass
            fractal_dimension = fractalDimension(lat_size, particles, center_x, center_y) # Calculates fractal dimension
            plt.close('all')

            end = time.time() # End time for internal clock

            # Writes all values to data base
            print('Writing')
            dict_of_elements = {field_names[0]: str(date), field_names[1]: (end - start),
                                field_names[2]: lat_size, field_names[3]: particles,
                                field_names[4]: fractal_dimension, field_names[5]: percolation,
                                field_names[6]: z_count[0], field_names[7]: z_count[1],
                                field_names[8]: z_count[2], field_names[9]: z_count[3],
                                field_names[10]: z_count[4], field_names[11]: z_count[5]}
            dict_writer.writerow(dict_of_elements)
            print(dict_of_elements)

            np.savetxt("D:\ASE II Resultados/Sin Solapamiento/" + str(i) + "cluster" + str(lat_size) + ', particle' + str(particles) + '.csv',
                       np.column_stack([x, y, z]), delimiter=',')


            Cluster.classReset()





mainScript(500, 12000, 100)
