"""
========================================================================================================================
Code: DLCA Cluster generator (OFF-LATTICE) in 2D
File: ClusterClass.py
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

# Cluster Class: In charge of saving all of the particles placed on the grid with all of their respective characteristics
# this includes their position, their mass, and has class variables firstp, nextp, and lastp which store the
# relationship between all of the particles.
class Cluster:
    # Class variables, will control the linked lists and the total number of clusters.
    number_of_clusters = 0  # Total number of clusters present in the Lattice
    firstp = []  # Indicates the first particle of the cluster
    nextp = []  # Linked list that will store the consecutive particles of the clusters
    lastp = []  # Final particle in the cluster
    grid = None # Grid variable will be utilized place and move all particles and will contain all clusters (indexes)
    particle_grid = None # Grid that contains all particles by number
    mass_list = dict()
    denominator = 0
    A = 0  # Constant of normalization for mass distribuition
    MAX_Z = 4 # Set maximum coordination number, set to 4 since the maximum is 4

    def __init__(self, x: float, y: float, cluster_number: int, mass: int):
        # Instance Variables, every particle that is placed in the lattice will have these variables.
        self.x = x  # Position x of particle
        self.y = y  # Position y of particle
        self.mass = mass  # Mass of particle
        # (only the first particle in a cluster will have mass of cluster)
        self.number = cluster_number # Number of cluster (goes from 0 to total particles - 1)
        self.side_particles = 0 # Number of particles connected to this specific particle

        try:
            Cluster.firstp[self.number_of_clusters] = cluster_number
            self.index = self.number_of_clusters
        except:
            Cluster.firstp.append(cluster_number)
            self.index = self.number_of_clusters

        Cluster.nextp.append(-1)

        try:
            Cluster.lastp[self.number_of_clusters] = cluster_number
        except:
            Cluster.lastp.append(cluster_number)

        Cluster.mass_list[mass] += 1

        Cluster.denominator += 1 / mass

        Cluster.number_of_clusters += 1

    def addSideParticle(self):
        self.side_particles += 1

    #This has to be run before the creation of the first instance of an object
    # of class particle to set all of the possible masses in the dictionary
    @classmethod
    def setMassDictionary(cls, total_particles: int):
        for i in range(1, total_particles + 1):
            cls.mass_list[i] = 0

    # Removes clusters from dictionary and adds the larger cluster to the dictionary
    @classmethod
    def changeMassList(cls, mass1: int, mass2: int):
        cls.mass_list[mass1] -= 1
        cls.mass_list[mass2] -= 1
        cls.mass_list[mass1 + mass2] += 1

    # Updates denominator of constant of normalization when clusters merge
    @classmethod
    def changeDenominator(cls, mass1: int, mass2: int):
        cls.denominator += -(1 / mass1) - (1 / mass2) + (1/(mass1 + mass2))

    # Classmethod utilized to set the value of the constant of normalization A for all clusters
    @classmethod
    def setA(cls):
        cls.A = 1/cls.denominator
        if cls.A <= 0:
            raise ValueError("Constant of integration is negative")

    # Classmethod utilized to reset the values of the class variables to 0, this is done to run the main code many times
    # and reset after every time the code is run.
    @classmethod
    def classReset(cls):
        # Class variables, will control the linked lists and the total number of clusters.
        cls.number_of_clusters = 0  # Total number of clusters present in the Lattice
        cls.firstp = []  # Indicates the first particle of the cluster
        cls.nextp = []  # Linked list that will store the consecutive particles of the clusters
        cls.lastp = []  # Final particle in the cluster
        cls.grid = None # Grid variable will be utilized place and move all particles and will contain all clusters (indexes)
        particle_grid = None  # Grid that contains all particles by number
        cls.mass_list = dict() # Resets dictionary of all masses
        cls.denominator = 0 # Resets denominator of constant of normalization (A)
        cls.A = 0  # Constant of normalization for mass distribution
        cls.MAX_Z = 4 # Set maximum coordination number, set to 4 since the maximum is 4
