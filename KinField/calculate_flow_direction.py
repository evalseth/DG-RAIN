import numpy as np
from matplotlib import cm
import matplotlib.tri as triClass
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
import sys


def calculate_gradient_field(f, triangulation, query_x, query_y):
	myTriInterp = triClass.LinearTriInterpolator(triangulation, f)
	(gx,gy) = myTriInterp.gradient(query_x, query_y)
	return gx, gy

def calculate_local_flow_dir(x, y, z, triangles):
	x = np.array(x)
	y = np.array(y)
	z = np.array(z)
	triangles = np.array(triangles)
	bar_x = np.mean(x[triangles],axis=1)
	bar_y = np.mean(y[triangles],axis=1)
	bar_z = np.mean(z[triangles], axis=1)
	
	triangulation = triClass.Triangulation(x, y, triangles=triangles)
	gx, gy = calculate_gradient_field(z, triangulation, bar_x, bar_y)

	flowDir = []
	u = []
	v = []

	triNum = 0
	for tri in triangles:
		grad_vec = np.array([gx[triNum], gy[triNum]])
		edg_vec = np.empty((2,1))
		dist = np.empty((3,1))
		el_flow_vec =[]
		for i in xrange(3):
			ind1 = i
			ind2 = (i+1)%3
			mid_x = 0.5*(x[tri[ind1]]+x[tri[ind2]])
			mid_y= 0.5*(y[tri[ind1]]+y[tri[ind2]])
			edg_vec[0] = mid_x - bar_x[triNum]
			edg_vec[1] = mid_y - bar_y[triNum]
			norm_edg_vec = norm(edg_vec)
			edg_vec = edg_vec/norm_edg_vec
			cos_angle = np.dot(grad_vec, edg_vec)/norm(grad_vec)
			#angles[i] = np.arccos(np.clip(cos_angle, -1, 1))
			dist[i] = 1 - cos_angle
			el_flow_vec.append(np.copy(edg_vec))
			#print "triNum = ", triNum, " mid_x = ", mid_x, "mid_y = ", mid_y 
			#print "triNum = ", triNum, " i = ", i, "edg_vec = ", edg_vec 
		
		min_ind = np.argmin(dist)
		#print "min_ind = ", min_ind
		flowDir.append(min_ind)

# collect some data for plotting purposes
		u.append(el_flow_vec[min_ind][0])
		v.append(el_flow_vec[min_ind][1])
		
		triNum += 1

## write out the barycenter and the flow direction vectors for plotting purposes
	f = open('./Output/FlowPaths.out', 'w')
	f.write('FlowPaths written by calculate_flow_directions.py\n')
	f.write('bary_x bary_y flowdirvec_x flowdirvec_y\n')
	for i in xrange(triNum):
		f.write('{}\t{}\t{}\t{}\n'.format(bar_x[i], bar_y[i], u[i][0], v[i][0]))
	f.close()

	return flowDir
	
