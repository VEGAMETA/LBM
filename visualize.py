import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# 3D-Scatter plot
def scatter_3d(x, y, z, values):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(x, y, z, c=values, cmap='viridis')
    plt.colorbar(sc, label='Values')
    plt.title("3D Visualization")
    plt.savefig("output.png", dpi=300)
    plt.show()
    
def scatter_2d(x, y, values):
    plt.scatter(x, y, c=values, cmap='viridis')
    plt.colorbar(label='Values')
    plt.title("2D Visualization")
    plt.savefig("output.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    data = np.loadtxt('density.dat')
    #x, y, z, values = data.T
    #scatter_3d(x, y, z, values)
    x, y, values = data.T
    scatter_2d(x, y, values)