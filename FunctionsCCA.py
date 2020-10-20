"""
========================================================================================================================
Code: CCA Cluster generator (OFF-LATTICE) in 2D
File: FunctionsCCA.py
Author: Diego Rosenberg

Description: This code generates an on-lattice DLCA cluster through a random collocation of particles and
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
import numpy as np
import random
from typing import Tuple, List
from ClusterClass import Cluster


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

# Function that will perform a step in a direction of a Cluster of a selected cluster.
def step(L, particle_list):
    sclust = random.randint(0, Cluster.number_of_clusters - 1)
    # Probability of that cluster moving
    prob_num = Cluster.mass_list[particle_list[Cluster.firstp[sclust]].mass]
    prob_den = Cluster.A / particle_list[Cluster.firstp[sclust]].mass
    prob = prob_num / prob_den
    z = random.random()
    if z > prob:
        # This if is utilized to limit the movement of the larger clusters in the set, if the random probability z is
        # smaller than the probability of a standardized constant A*cluster_mass then the particle will take a step
        # if not the particle chosen for movement will be discarded and the process will continue
        return particle_list, False

    dx = (1, 0, -1, 0)
    dy = (0, 1, 0, -1)
    dir = random.randint(0,3)
    particle = Cluster.firstp[sclust]
    while particle != -1:
        Cluster.grid[particle_list[particle].x, particle_list[particle].y] = -1  # Resets value of previous position to -1
        particle_list[particle].x += dx[dir]
        particle_list[particle].y += dy[dir]
        # All of these ifs are created so that if a particle is to overstep, it wraps around to the other end of the Cluster.
        if particle_list[particle].x == L:
            particle_list[particle].x = 0
        elif particle_list[particle].y == L:
            particle_list[particle].y = 0
        elif particle_list[particle].x == -1:
            particle_list[particle].x = L-1
        elif particle_list[particle].y == -1:
            particle_list[particle].y = L-1

        particle = Cluster.nextp[particle]

    particle = Cluster.firstp[sclust]
    while particle != -1:
        # Relabels new position to the name of the cluster
        Cluster.grid[particle_list[particle].x, particle_list[particle].y] = sclust
        particle = Cluster.nextp[particle]

    particle = Cluster.firstp[sclust]
    while particle != -1:
        # Checks every particle in cluster for neighbor from other cluster
        checkCluster(L, particle_list[particle], particle_list)

        particle = Cluster.nextp[particle]

    return particle_list, True



# Checks particle to nearby particles (including wrapped borders)
def checkCluster(L, particle, particle_list):
    dx = (1, 0, -1, 0)
    dy = (0, 1, 0, -1)
    for i, j in zip(dx, dy):
        px, py = particle.x + i, particle.y + j
        if px == L:
            px = 0
        elif py == L:
            py = 0
        elif px == -1:
            px = L-1
        elif py == -1:
            py = L-1
        if Cluster.grid[px, py] != -1 and Cluster.grid[px, py] != Cluster.grid[particle.x, particle.y]:
            # particle_list[particle.number].addSideParticle()
            # particle_list[Cluster.grid[px, py]].addSideParticle()
            joinClusters(Cluster.grid[px, py], Cluster.grid[particle.x, particle.y], particle_list)


# Function in charge of merging the two clusters that have collided in the simulation.
# Function in charge of merging the two clusters that have collided in the simulation.
def joinClusters(c1, c2, particle_list):
    # lc_label = larger cluster label and, sc_label = smaller cluster label
    if particle_list[Cluster.firstp[c1]].mass > particle_list[Cluster.firstp[c2]].mass:
        lc_label = c1
        sc_label = c2
    else:
        lc_label = c2
        sc_label = c1

    Cluster.changeDenominator(particle_list[Cluster.firstp[lc_label]].mass, particle_list[Cluster.firstp[sc_label]].mass)
    Cluster.changeMassList(particle_list[Cluster.firstp[lc_label]].mass, particle_list[Cluster.firstp[sc_label]].mass)
    # Manipulation of linked lists
    Cluster.nextp[Cluster.lastp[lc_label]] = Cluster.firstp[sc_label]
    Cluster.lastp[lc_label] = Cluster.lastp[sc_label]
    particle_list[Cluster.firstp[lc_label]].mass += particle_list[Cluster.firstp[sc_label]].mass
    current = Cluster.firstp[sc_label]
    while current != -1:
        Cluster.grid[particle_list[current].x, particle_list[current].y] = lc_label
        particle_list[current].index = lc_label
        current = Cluster.nextp[current]

    Cluster.number_of_clusters -= 1

    # Resets label of the grid to the smaller of the two labels to have the last cluster existing be a 0 cluster.
    if sc_label != Cluster.number_of_clusters:
        current = Cluster.firstp[Cluster.number_of_clusters]
        while current != -1:
            Cluster.grid[particle_list[current].x, particle_list[current].y] = sc_label
            particle_list[current].index = sc_label
            current = Cluster.nextp[current]

        particle_list[Cluster.firstp[sc_label]].mass = particle_list[Cluster.number_of_clusters].mass
        Cluster.firstp[sc_label] = Cluster.firstp[Cluster.number_of_clusters]
        Cluster.lastp[sc_label] = Cluster.lastp[Cluster.number_of_clusters]


    Cluster.setA()



# Gets the center of mass of a finished system
# Formula found in article: Calculating Center of Mass in an Unbounded 2D Environment
# Authors: Linge Bai and David E. Breen
# DOI: 10.1080/2151237X.2008.10129266
def centerOfMass(L: int, particles: int, *args) -> tuple:
    if len(args) == 1:
        arr = args[0]
        x, y = np.hsplit(arr, 2)
        x = x.flatten()
        y = y.flatten()
    elif len(args) == 2:
        x = args[0]
        y = args[1]
    else:
        raise ValueError("Incorrect array entered for x and y.")

    epsilon_x, zeta_x, epsilon_y, zeta_y = 0, 0, 0, 0
    for (xi, yi) in zip(x, y):
        epsilon_x += np.cos((xi / L) * 2 * np.pi)
        zeta_x += np.sin((xi / L) * 2 * np.pi)
        epsilon_y += np.cos((yi / L) * 2 * np.pi)
        zeta_y += np.sin((yi / L) * 2 * np.pi)
    epsilon_x = epsilon_x / particles
    zeta_x = zeta_x / particles
    epsilon_y = epsilon_y / particles
    zeta_y = zeta_y / particles

    theta_x = np.arctan2(-zeta_x, -epsilon_x) + np.pi
    theta_y = np.arctan2(-zeta_y, -epsilon_y) + np.pi

    center_x = (theta_x / (2 * np.pi)) * L
    center_y = (theta_y / (2 * np.pi)) * L

    return center_x, center_y


# Checks completed cluster to see if cluster has percolated lattice
# REVIEW AND CHANGE
"""
def clusterPercolates(L: int, particle_list: List[Cluster]) -> bool:
    percolation = False
    x_percolation = True
    y_percolation = True

    i = 0
    while True:
        if x_percolation:
            break_condition_x = False  # Ends one instance of wile and for loop to move on to the next column
            for j in range(Cluster.cells):  # This exterior for will check in every column space for a row
                current_kx_left = (j * Cluster.cells) + i
                for num, particle in enumerate(particle_list):  # Checks every particle for the given column
                    if particle.k == current_kx_left:
                        break_condition_x = True  # Breaks loop
                        break
                if break_condition_x:
                    break

            if not break_condition_x:
                x_percolation = False

        if y_percolation:
            break_condition_y = False  # Ends one instance of wile and for loop to move on to the next row
            for j in range(Cluster.cells):  # This exterior for will check in every row space for a column
                current_ky_bottom = (i * Cluster.cells) + j
                for num, particle in enumerate(particle_list):
                    if particle.k == current_ky_bottom:  # Breaks loop
                        break_condition_y = True
                        break
                if break_condition_y:
                    break

            if not break_condition_y:
                y_percolation = False

        if x_percolation is False and y_percolation is False:
            break

        if i > L:
            percolation = True
            break

        i += 1

    return percolation
"""

if __name__ == '__main__':
    pass
