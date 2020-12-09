"""
========================================================================================================================
Code: Fractal Dimension for a 2D Cluster that is off-cetner
File: FractalDimensionCCA.py
Author: Diego Rosenberg

Description: This code will calculate the fractal dimension for a 2D cluster that has its center of
mass away from the center of the grid.
========================================================================================================================
"""

# -----------------------------------------------------------------------------
# Imported libraries
# -----------------------------------------------------------------------------
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, List
from LinearReg import linearReg
from FunctionsDLCA import centerOfMass, surroundingSquare


# -----------------------------------------------------------------------------
# Main Code
# -----------------------------------------------------------------------------

# Function that counts the average number of particles per box in different divisions, calculates the size of these boxes,
# and calculates the amount of boxes that contain particles for a given division size.
def boxCount(L: int, particles: int, x: List[float], y: List[float]) -> Tuple[list, list, list]:
    div = 4 + (4 * (particles <= 8000)) # Sets start division
    break_condition = 2 + (2 * (L > 100)) # Sets ending place for amount of divisions
    box_count = [] # List containing the amount of boxes with particles
    avg_counts = [] # List containing the amount of particles average per box
    box_size = [] # List containing the size of the current boxes
    min_x, max_x, min_y, max_y, x, y = surroundingSquare(L, x, y)
    while True: # Emulation of do while loop to calculate individual vales to put into previous lists
        boxcount = 0
        division_x = (max_x - min_x) / div
        division_y = (max_y - min_y) / div
        division = (division_x + division_y) / 2
        for i in range(1,div+1):
            for j in range(1,div+1):
                lower_x = ((i - 1) * division) + min_x
                upper_x = (i * division) + min_x
                lower_y = ((j-1) * division) + min_y
                upper_y = (j * division) + min_y
                for xi, yi in zip(x, y):
                    if lower_x <= xi <= upper_x and lower_y <= yi <= upper_y:
                        boxcount += 1
                        break # If a box contains 1 particle the code can move on to the next

        # Appends found values to the lists
        box_count.append(boxcount)
        avg_counts.append(particles/boxcount)
        box_size.append(division)

        print("{:.3f} divisions with minimum {}".format((division_x + division_y) / 2, break_condition))
        if division < break_condition: # Ends do while loop
            break

        div *= 2

    return box_count, avg_counts, box_size

# Function that finds the relative radius of all particles with respect to the center of mass,
# and will calculate the amount of particles inside certain radius sizes.
def radius(L: int, x: List[float], y: List[float], cx: float, cy: float) -> Tuple[list, list]:
    particle_radius = [] # List that will contain the radius of all particles in the cluster
    halfL = L/2
    for xi, yi in zip(x, y):
        dx = abs(xi - cx)
        dy = abs(yi - cy)

        if dx > halfL:
            dx = L - dx
        if dy > halfL:
            dy = L - dy

        particle_radius.append(np.sqrt(dx**2 + dy**2))

    # Returned lists that calculate different radii and the amount of particles inside that radius.
    radius = [halfL/32, halfL/16, halfL/8, halfL/4, halfL/2]
    count = [len([1 for r in particle_radius if r <= halfL/32]),
             len([1 for r in particle_radius if r <= halfL/16]), len([1 for r in particle_radius if r <= halfL/8]),
             len([1 for r in particle_radius if r <= halfL/4]), len([1 for r in particle_radius if r <= halfL/2])]

    return radius, count


# Main function that will calculate fractal dimension with 3 techniques:
#   1. The slope of the linear regression between logarithm of the box size vs logarithm the average particles per box.
#   2. The absolute value of the slope of the linear regression of the logarithm box size vs logarithm number of boxes.
#   3. The value of the slope of the linear of the linear regression between the logarithm radius vs logarithm of
#       amount of particles.
def fractalDimension(lat_size: int, particles: int, center_x = None, center_y = None) -> float:
    # If entered lat size and particles is not found as a file it will run the main code.
    if not os.path.isfile("cluster"+str(lat_size)+", particle"+str(particles)+".csv"):
        from MainDLCA import main
        main(lat_size, particles)

    # Reads created csv file and remembers both columns as x and y
    cluster = np.loadtxt("cluster"+str(lat_size)+", particle"+str(particles)+".csv", delimiter=',')
    x, y, *_ = np.hsplit(cluster, cluster.shape[1])
    x = x.flatten()
    y = y.flatten()

    # Run box count function and radius function to have all variables
    box_count, avg_count, box_size = boxCount(lat_size, particles, x, y)
    if (center_x is None) and (center_y is None):
        center_x, center_y = centerOfMass(lat_size, particles, x, y)
    # rads, rad_count = radius(lat_size, x, y, center_x, center_y)

    # Gets logarithm of all the interested lists
    log_box_count = np.log(box_count)
    log_avg_count = np.log(avg_count)
    log_box_size = np.log(box_size)
    # log_radius = np.log(rads)
    # log_radius_count = np.log(rad_count)

    # PLOTS BOX_SIZE VS AVG. PARTICLES PER BOX
    m_box_pcount, b_box_pcount = linearReg(log_box_size, log_avg_count)
    line_x_box_pcount = np.linspace(0,max(log_box_size),num=100)
    line_y_box_pcount = m_box_pcount * line_x_box_pcount + b_box_pcount

    plt.figure()
    plt.scatter(log_box_size, log_avg_count, c='black', marker='.', label='Values')
    plt.plot(line_x_box_pcount, line_y_box_pcount, c = 'royalblue', ls = "--", label='Linear Fit')
    plt.title("Plot of Box Size vs Average Particles per Box (In log scale)")
    plt.xlabel("Log Box Size")
    plt.ylabel("Log Average Particles per Box")
    plt.legend()
    plt.text(0.7 * np.max(log_box_size), 0.2 * np.max(log_avg_count),
             "Ecuation:\ny = {:.3f}x + {:.3f}".format(m_box_pcount, b_box_pcount))

    # PLOT BOX_SIZE VS NUMBER OF BOXES
    m_boxes, b_boxes = linearReg(log_box_size, log_box_count)
    line_x_boxes = np.linspace(0,max(log_box_size),num=100)
    line_y_boxes = m_boxes * line_x_boxes + b_boxes

    plt.figure()
    plt.scatter(log_box_size, log_box_count, c="black", marker='.', label = "Values")
    plt.plot(line_x_boxes, line_y_boxes, c = "royalblue", ls = "--", label = "Linear Fit")
    plt.title("Plot of Box Size vs Number of Boxes (In log scale)")
    plt.xlabel("Log Box Size")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.text(0.7 * np.max(log_box_size), 0.8 * np.max(log_box_count),
             "Ecuation:\ny = {:.3f}x + {:.3f}".format(m_boxes, b_boxes))

    # PLOT RADIUS SIZE VS AMOUNT OF PARTICLES IN RAIDUS
    # m_radius, b_radius = linearReg(log_radius, log_radius_count)
    # line_x_radius = np.linspace(0, max(log_radius), num=100)
    # line_y_radius = m_radius * line_x_radius + b_radius
    #
    # plt.figure()
    # plt.scatter(log_radius, log_radius_count, c='black', marker='.', label='Values')
    # plt.plot(line_x_radius, line_y_radius, 'g--', label="Aproximation")
    # plt.title("Plot of Radius vs Mass (In log scale)")
    # plt.xlabel("Log Radius")
    # plt.ylabel("Log Mass")
    # plt.legend()
    # plt.text(0.5 * max(log_box_size), 0.20 * max(log_radius_count),
    #          "Ecuation:\ny = {:.3f}x + {:.3f}".format(m_radius, b_radius))

    fractal_dimension = m_box_pcount
    return fractal_dimension


if __name__ == '__main__':
    start = time.time()
    fractalDimension(500, 1000)
    end = time.time()
    print(end - start)
    plt.show()
