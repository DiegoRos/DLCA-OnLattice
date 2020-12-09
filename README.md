# DLCA Off-Lattice Project
  ## Author: Diego Rosenberg

## Description 
This code generates a DLCA cluster through a random collocation of particles and a random movement of these particles (which is dependent on their mass) to generate larger clusters. The movement of these particles is pseodorandom and is dependant on their mass (as would happen in a particle dynamics simulation). Along with this a restriction can be added so that particles can only stick if they have a coordination number of 0 or 1. On these clusters we measure their fractal dimension and the probability of the cluster percolating through the lattice.

## Code
This project is divided into three main files in charge of creting the cluster, one for measuring the fractal dimension of resultant clusters and another file responsible for measuring the probability of the cluster percolating.

- Creation of CCA Cluster:
  - [x] MainDLCA.py
  - [x] FunctionsDLCA.py
  - [x] ClusterClass.py
- Measurment of Fractal Dimension: 
  - [x] FractalDimensionDLCA.py
- Measurment of Fractal Percolation:
  - [x] clusterPercolates function in FunctionsDLCA.py 

## Files
### ClusterClass.py
- Contains the programming description of the Clusters that are formed this .py file will contain the ring size for the particles along with the linked lists that allow clusters to join each other.
- Imported Libraries:
  - NONE
- Classes:
  - Cluster
    - Class Variables:
      - number_of_clusters: Total number of clusters present in the current system.
      - firstp: Indicates the first particle of the cluster **(Note: if viewed in the middle it will repeat values)**
      - nextp: Linked list that will store all of the consecutive values for the clusters.
      - lastp: Last particle for all clusters **(Note: this will repeat values since none are removed after stepping and merging)**
      - cell_list: Linked list that contains the current position for all of the particles in the grid.
      - mass_list: Dictionary that contains the amount of clusters of mass m in the cluster
      - cells: Number of cells horizontally and vertically
      - cells2: Total number of cells in the grid
      - ring_size: Contact size of the particles in the system.
      - denominator: The denominator of the constant of normalization
      - A: Constant of normalization for mass distribuition.
    - Instance Variables:
      - x: Position x of particle
      - y: Position y of particle
      - mass: Mass of cluster **(Note: This value will only be the "real mass" for the leading particle in the cluster)**
      - number: number of cluster (goes from 0 to particles - 1)
      - index: Position in the list firstp of the cluster
      - k: Current position in grid of particle.
    - Init Method:
      - **args**: x, y, cluster_number, mass, L, ring_size
      - **kwargs**: NONE
      - This function is called when creating an instance of the class Cluster. This function will set all of the instance variables to the corresponding particle and will add itself to the corresponding class variables and lists.
    - Instance Methods: (Functions performed over an instance of Cluster class)
      - setK
        - **args**: NONE
        - **kwargs**: NONE
        - Sets k value for a given instance of Cluster.
      - addToCellList
        - Adds instance of Cluster to cell_list utilizing the k value of the instance, if the space is occupied it will add it in the position cells2 + number of the current particle
      - resetCellListElement
        - **args**: NONE
        - **kwargs**: NONE
        - This function will remove the current value from the cell_list and shift the remaining linked values one position to the left.
    - Class Methods:
      - setMassDictionary:
      - **args**: total_particles
      - **kwargs**: NONE
      - This method sets the initial value of the mass_list, to do this it will create a dictionary with keys equal to all of the possible values of the mass of a cluster (ie. 1 to # of total particles).
      - **Note: this function is called only once at the beginning of the code.**
      - changeMassList
        - **args**: mass1, mass2
        - **kwargs**: NONE
        - This function is called when a merge between two clusters is going to occur and therefore removes one cluster with mass $m_1$ and $m_2$ and adds a cluster with a mass of $m_1 + m_2$.
      - changeDenominator
        - **args**: mass1, mass2
        - **kwargs**: NONE
        - This function is also called when a merge between two clusters is about to occur, this removes the corresponding fractions from the denominator to keep the constant of normalization up to date with the changing clusters. In this function from the denominator we perform the following mathematical procedure:
        $denominator = -\frac{1}{m_1} -\frac{1}{m_2} + \frac{1}{m_1 + m_2}$
      - setA
        - **args**: A
        - **kwargs**: NONE
        - Set the constant of normalization for the probability distribuiton. This is equal to: $A = \frac{1}{denominator}$ and therefore will be set utilizing the class variable denominator.
      - classReset:
        - **args**: NONE
        - **kwargs**: NONE
        - Will reset the class values of the cluster (used if multiple iterations are run in the same file)

### FunctionsCCA.py
-  Contains the functions necessary to move the individual clusters and join the clusters together.
- Imported Libraries:
  - numpy, random, ClusterClass
-  Functions:
   -  checkSpot
      -  **args**: current_x, current_y, ring_size, particle_list, k_list
      -  **kwargs**: NONE
      -  Checks the grid to se if any of the particles are within merging distance (ring_size) and if that is the case the particle will not be placed in that position and will look for another. This is similar to diffusing without aggregation until a condition is met.
      -  returns True or False
   - step
     - **args**: L, particle_list, k_list
     - **kwargs**: NONE
     - Selects a random cluster, calculates the probability of that cluster moving and chooses wheter to move or not. It compares the probability of the cluster moving with a random number selected from 0 to 1, if the cluster does not move the function returns the particle_list and False. On the other hand the code will call the move function without any kwargs and the selected direction will be saved from the return value from this function. After this all of the particles in the cluster are run through the checkCluster function and this will return all of the possible clusters that can be connected to the current cluster. From this list of clusters only the one with the smallest distance from one particle to anoter to be connected (this is done in this way to prevent clusters to overlap). Once this minimum distance is found the current cluster will be pushed back in the same direction that it came, utilizing the move function and the direction and distance we want to push back as kwargs, in (multiplied by a value of 1-distance), eliminating the overlap between clusters. 
     - returns particle_list, (True or False)
    - move
    - **args**: L, sclust, particle_list
    - **kwargs**: dir, dist
    - This function is in charge of moving clusters in one of two ways, either it moves the cluster in a random direction (if dir and dist are None) or move back in the direction and magnitude sent in by the kwargs.
   - checkCluster
     - **args**: L, k_list, particle, particle_list
     - **kwargs**: NONE
     - This will take a particle and check its surroundings (with k_list) and check if it is within its ring_size a particle of another cluster. If another particle is within the ring_size the function will merge the clusters by calling joinClusters.
     - returns (particle1 number, particle1 index, particle2 number, particle2 index, distance)
   - returnDist
     - **args**: dist, dir, p0, p1
     - **kwargs**: NONE
     - This function will calculate the magnitude of the distance (on the axis of movement) that the moved cluster has to be returned to remove overlapping.
     - returns rdist
       - **NOTE: This function returns 1 - dist that the particle has to be returned (that is why in move the returned distance is dist - Cluster.ring_size)** 
   - joinClusters
     - **args**: c1, c2, particle_list
     - **kwargs**: NONE
     - This function joins two clusters utilizing the linked lists found in Cluster class and changes the index of the small cluster to match the index of the larger cluster and sets the new mass of the cluster **(Note: only the first element of the cluster will have the mass of the entire cluster)**, along with this the denominator class variable, and the mass_list class variable in Cluster will be changed accordingly to have the ability to calculate the new A. At the end of the function the new A is set with the changed variables
     - returns None
    - posParticlesCluster
      - **args**: index, particle_list
      - **kwargs**: NONE
      - This function can be called to retrieve the position of a given cluster index at any point during the simulation. This is utilized to do manipulations over single clusters during the simulation.
    - klist
      - **args**: L, k
      - **kwargs**: NONE
      - This function will calculate all of the possible combinations for any k (cell value) in the entire grid. For any given k we have the next list:
        - [ Cell bellow to the left **,** Cell bellow **,** Cell bellow to the right **,** Cell to the left **,** Current cell (k) **,** Cell to the right **,** Cell above to the left **,** Cell above **,** Cell above to the right ]
      - returns a single k_list (like the one above)
      - **Note: This function is called only at the begining of the code to generate a long list of lists full of all of the possible combinatios of k in k_list.**
    - centerOfMass
      - **args**: L, particles, *args (list of x and y values or numpy array of x and y values)
      - **kwargs**: NONE
      - This function calculates the center of mass for a cluster existing inside a lattice with periodic boudary conditions utilizing the formula found in article:Calculating Center of Mass in an Unbounded 2D Environment. Authors: Linge Bai and David E. Breen. DOI: 10.1080/2151237X.2008.10129266
      - returns two floating values representing the center of mass in x and the center of mass in y.
    - clusterPercolates
      - **args**: L, particle_list
      - **kwargs**: NONE
      - This function will check every row and every column of the system to revise if a given cluster percolates the system, if all of the columns are filled by at least one particle the cluster will percolate in x, and if all rows are filled by at least one particle the cluster will percolate in the y direction.
    - surroundingSquare
      - **args**: lat_size, x, y
      - **kwargs**: cx (Default value is None), cy (Default value is None)
      - This function thakes the lists of x and y, moves the center of mass of this cluster to the center of the lattice and finds a square that can surround the entire cluster. **Note: This function does not take into account periodic boundary conditions.**
    - moveToCenter
      - **args**: lat_size, x, y
      - **kwargs**: cx (Default value is None), cy (Default value is None)
      - Moves the center of mass of the cluster to the center of the lattice, this can be usefull to plot clusters in a better way.
    - radiusOfGyration
      - **args**: lat_size, particles, x, y
      - **kwargs**: cx (Default value is None), cy (Default value is None)
      - Utilizing the formula $Rg^2 = \frac{1}{N} \sum_{i = i}^{N}(r_i - r_c)^2$, this function calculates the radius of gyration for any given cluster.


### MainCCA.py:
- Contains the main code and is responsible for generating the clusters.
- Imported Libraries:
  - os, shutil, time, matplotlib.pyplot, immageio, ClusterClass, and FunctionsCCA
- Functions:
  - initialize
    - **args**: L, particles, dist
    - **kwargs**: NONE
    - This function creates the necesary constants like particle ring size, the amount of cells horizontally and vertically in the lattice, the total number of cells and sets the initial cell list. Along with that creates k list and places every single particle in the grid and checks if they are adjacent to another particle.
    - returns particle_list and k_list
  - plot
    - **args**: L, particles, particle_list, num_steps
    - **kwargs**: animation (default value is False)
    - Plot is responsable for the animation plotting of each step in the grid (if an animation is required) and of the final plot of the system (which always must have animation = False). All of these will plot circles with radius half the size of the Cluster.ring_size.
    - returns x, y (lists containing all particles in grid)
  - main
    - **args**: lat_size, particles
    - **kwargs**: contact_distance (default value is 1), animation (default value is False), images_per_frame (default value is 10)
    - This function takes care of the recursivity necesary to create the cluster. The first part of the code will create an empty directory called "images" in the path that the MainCCA.py file is placed **(Note: if file explorer is open and animation is True there will be an error)**. Next the system is initialized by calling initialize and next the and call the step function from FunctionsCCA.py until there is only one cluster remaining, and if animation is True it will plot every images_per_frame steps. After that it plots and cretes final animation if animation was True.
    - return None

### FractalDimensionCCA.py
- This code will calculate the fractal dimension for a 2D cluster that is off-lattice and has its center of mass away from the center of the grid. This will be done utilizing three distinct methods and the average from these three methods will be taken. The three methods are:
  1. The slope of the linear regression between logarithm of the box size vs logarithm the average particles per box.
  2. The absolute value of the slope of the linear regression of the logarithm box size vs logarithm number of boxes.
  3. The value of the slope of the linear of the linear regression between the logarithm radius vs logarithm of amount of particles.
- Imoprted Libraries:
  - os, time, matplotlib.pyplot, numpy, LinearReg, FuncitonsCCA
- Functions:
  - boxCount
    - **args**: L, particles, x, y
    - **kwargs**: NONE
    - This function calculates the averagae number of particles per box in different divisions, the size of these boxes, and calculates the amount of boxes that contain particles for a given division size.
    - returns 3 different lists the first contains the boxes counted with different divisons, the second the average amount of particles for each division, and lastly the division sizes.
  - radius
    - **args**: L, x, y, cx, cy
    - **kwargs**: NONE
    - This function calculates the relative radius of all of the particles with respect to their center of mass and calculates the amount (or mass) of particles inside various radii.
    - returns two lists one which contains various radii of importance and another list which contains the amount of particles inside of the radius stated previously.
  - fractalDimension
    - **args**: lat_size, particles
    - **kwargs**: center_x = None, center_y = None
    - Main function that calculates the fractal dimension utilizing the previous functions and the previus techniques. The slope is calculated utilizing a linear regression creted for this project.
    - returns a floating point number which is the average of the three obtained fractal dimensions.
  
### LinearReg.py
- Contains the function to perform a linear regression.
- Imported Libraries:
  - NONE
- Functions:
  - linearReg
    - **args**: xlist, ylist
    - **kwargs**: NONE
    - This function sill find the slope and intercept for a given x and y lists of integers.
    - returns 2 floating point values which are the slope and intercept respectively.