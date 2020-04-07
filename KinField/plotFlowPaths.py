import numpy as np
import matplotlib.tri as triClass
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_grid_data(path):
    r"""
    Read a grid specification file (fort.14)

    :Note: This function currently only reads in the grid itself and not the 
           boundary specifications.

    :Input:
     - *path* (string) - Path to grid data file.

    :Output:
     - (numpy.ndarray(num_nodes, 2)) - Coordinates of each node.
     - (numpy.ndarray(num_nodes)) - Array containing depths at each node.
     - (numpy.ndarray(num_elements, 3)) - Numpy array containg specifcation of 
       each element's nodal points in counter-clockwise fashion.
    """

    grid_file = open(path, 'r')

    # Read header
    grid_file.readline()
    num_elements, num_nodes = [int(value) for value in grid_file.readline().split()]
    
    # Create necessary array storage
    coords = np.empty((num_nodes,2))
    depth = np.empty((num_nodes))
    triangles = np.empty((num_elements,3),dtype='int32')

    # Read in coordinates of each node
    for n in xrange(num_nodes):
        line = grid_file.readline().split()
        coords[n,0] = float(line[1])
        coords[n,1] = float(line[2])
        depth[n] = float(line[3])

    # Read in triangles
    for n in xrange(num_elements):
        line = grid_file.readline().split()
        # Construct triangles
        triangles[n,0] = int(line[2]) - 1
        triangles[n,1] = int(line[3]) - 1
        triangles[n,2] = int(line[4]) - 1
    # Read in boundary data

    grid_file.close()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(coords[:,0], coords[:,1], -depth)
    plt.show()

    return num_elements, num_nodes, coords, depth, triangles


def plot_flow_paths(gridFilePath, FlowPathFilePath):
	flowvec = np.loadtxt(FlowPathFilePath)
	num_elements, num_nodes, coords, depth, triangles = read_grid_data(gridFilePath)
	triangulation = triClass.Triangulation(coords[:,0], coords[:,1], triangles=triangles)
	bar_x = np.mean(coords[triangles,0], axis=1)
	bar_y = np.mean(coords[triangles,1], axis=1)

	u = flowvec[:,0] - bar_x
	v = flowvec[:,1] - bar_y

	norm = np.sqrt(np.dot(u,u) + np.dot(v,v))
	u = u/norm
	v = v/norm

	u = 0.1*u
	v = 0.1*v

	fig = plt.figure()
	#plt.gca().set_aspect('equal')
	plt.triplot(triangulation)
	plt.quiver(bar_x, bar_y, u, v, units='xy', color='blue', width=1.0, headwidth=3.0, headlength=4.0)
	plt.show()



plot_flow_paths('./ParkingLotInfiltration/ParkingLotMesh.14', './Output/FlowPaths.out')

