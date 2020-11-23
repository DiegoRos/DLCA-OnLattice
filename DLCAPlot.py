from FunctionsCCA import moveToCenter
import matplotlib.pyplot as plt
import numpy as np


def plot(lat_size, particles, file_name, center = False):
    cluster = np.loadtxt(file_name, delimiter = ',')

    x, y, *_ = np.hsplit(cluster, cluster.shape[1])
    x = x.flatten().astype(int)
    y = y.flatten().astype(int)

    if center:
        x, y = moveToCenter(lat_size, x, y)


    grid = np.full((lat_size, lat_size), -1)

    for xi, yi in zip(x, y):
        grid[int(xi)][int(yi)] = 0

    # Plots final cluster and saves a .csv file and a png for future use
    plt.figure(figsize= (9,9))
    plt.title("CCA Cluster L = {}, Particles = {}".format(lat_size, particles), fontsize = 16)
    plt.xlabel("x", fontsize = 16)
    plt.ylabel("y", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.imshow(grid, cmap = "Greys", origin='lower')

if __name__ == "__main__":
    lat_size = 50
    particles = 200
    file_name = f"Partial Results/Partialcluster{lat_size}, particle{particles}.csv"

    plot(lat_size, particles, file_name)
    plt.show()
