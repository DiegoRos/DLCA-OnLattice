"""
========================================================================================================================
Code: CCA Cluster generator in 2D
File: mainCCA.py
Author: Diego Rosenberg

Description: This code generates a CCA cluster through a random collocation of particles and
a random movement of these particles to generate larger clusters. This code is separated into
three separate files containing various functions.
    - mainCCA.py: Contains the main code and is responsible for generating the clusters
    - FunctionsCCA.py: Contains the functions necessary to move the individual clusters and join
        the clusters together.
    - ClusterClass.py: Contains the programming description of the Clusters that are formed
        this .py file will contain the ring size for the particles along with the linked lists
        that allow clusters to join each other.
========================================================================================================================
"""

# -----------------------------------------------------------------------------
# Imported libraries
# -----------------------------------------------------------------------------
import os
import shutil
import time
import matplotlib.pyplot as plt
import imageio
from typing import Tuple, List
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from ClusterClass import Cluster
from FunctionsCCA import *  # Imports all functions from FunctionsCCA.py


# -----------------------------------------------------------------------------
# Main Code
# -----------------------------------------------------------------------------
# Function that initializes the grid and places all particles in it.
def initialize(L: int, particles: int) -> List[Cluster]:
    # Sets class variables that will be utilized during the
    Cluster.setMassDictionary(particles)
    Cluster.grid = np.full((L, L), -1)
    Cluster.particle_grid = np.full((L, L), -1)

    particle_list = []
    for i in range(particles):
        # while True: if: break is creted to emulate a do while loop present in other languages
        while True:  # Places all particles in a empty space inside of the lattice
            x = random.randint(0, L - 1)  # Selects a number from 0 to L-1 (inclusive)
            y = random.randint(0, L - 1)  # Selects a number from 0 to L-1 (inclusive)
            if Cluster.grid[x, y] == -1:
                break

        Cluster.grid[x, y] = Cluster.number_of_clusters  # Places particle in the grid
        Cluster.particle_grid[x, y] = i

        # Creates an instance of the Cluster class and creates a "particle"
        particle = Cluster(x, y, i, 1)

        # Appends the instance of the particle to a list containing all instances (particles)
        particle_list.append(particle)

        Cluster.setA()

        # Checks adjacent spots for particles
        checkCluster(L, particle, particle_list)

    return particle_list


# Function that plots the complete (or incomplete cluster)
def plot(L, particles, num_steps, animation = False):
    if not animation:
        # Plots final cluster and saves a .csv file and a png for future use
        plt.figure(figsize=(8,8))
        plt.title("CCA Cluster L = {}, Particles = {}".format(L, particles))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.imshow(Cluster.grid, cmap = "Greys", origin='lower')
        plt.savefig("L{}P{}.png".format(L, particles), dpi=200)
        if os.path.isdir('images'):
            plt.savefig("images/cluster.png", dpi=200)

    # Plots all intervals after initializing (where a step is performed) the cluster to create a the images necessary
    # for creating an animation
    if animation:
        label = str(num_steps)
        plt.title("CCA Cluster L = {}, Particles = {}".format(L, particles))
        plt.imshow(Cluster.grid, cmap='Greys')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig("images/cluster{}.png".format(label), dpi=200)
        plt.close()


# Main function, runs entire code and generates a COMPLETE CCA cluster.
def main(lat_size: int, particles: int, animation = False, images_per_frame = 10) -> Tuple[List, List, bool]:
    if animation:
        if not os.path.isdir("images"):
            os.mkdir("images")
        else:
            shutil.rmtree('images') # Eliminates existing path with (possible) previous images, will not work if folder open
            os.mkdir('images')

    particle_list = initialize(lat_size, particles)

    num_steps = 0
    while Cluster.number_of_clusters > 1:
        particle_list, st = step(lat_size, particle_list)
        if st:
            if animation and num_steps % images_per_frame == 0:
                plot(lat_size, particles, num_steps, animation = animation)
            num_steps += 1

        if num_steps % 100 == 0:
            print('Total number of existing Clusters:', Cluster.number_of_clusters)

    print('Total number of existing Clusters:', Cluster.number_of_clusters)

    # Plots final result
    plot(lat_size, particles, num_steps)

    percolation = clusterPercolates(lat_size)

    # Creates animation and deletes all previous images as well.
    if animation:
        with imageio.get_writer('images/movie.gif', mode='I') as writer:
            for i in range(0, num_steps, images_per_frame):
                filename = "images/cluster" + str(i) + ".png"
                image = imageio.imread(filename)
                writer.append_data(image)
                os.remove(filename)
            image = imageio.imread("images/cluster.png")
            for i in range(15):
                writer.append_data(image)
            os.remove("images/cluster.png")

    return [particle.x for particle in particle_list], [particle.y for particle in particle_list], percolation


if __name__ == "__main__":
    start = time.time()
    # random.seed(783732)
    main(500, 10000)
    # main(50,200)
    # main(6, 15)
    finish = time.time()
    # Displays in console the time taken to complete the cluster
    print(finish - start)
    plt.show()  # Shows final cluster
