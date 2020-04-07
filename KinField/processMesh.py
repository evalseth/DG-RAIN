from channel_module import *
from junction_module import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from meshpy import triangle
import matplotlib.tri as tri
import numpy as np
import numpy.linalg as la
import scipy
from scipy.interpolate import interp2d

def draw_junction(junctionPolygon, figName):
	plt.axes()
	xmax = -1000000000000000.0
	ymax = -1000000000000000.0
	xmin = 1000000000000000.0
	ymin = 1000000000000000.0

	points = []
	for vertex in junctionPolygon.vertices:
		points.append([vertex[0], vertex[1]])
		xmax = max(vertex[0], xmax)
		ymax = max(vertex[1], ymax)
		xmin = min(vertex[0], xmin)
		ymin = min(vertex[1], ymin)

	polygon = plt.Polygon(points, fill=None, edgecolor='r')
	plt.gca().add_patch(polygon)
	plt.axis([xmin-1,xmax+1,ymin-1,ymax+1])
	plt.show()
	#plt.savefig(figName)
	plt.close()


def draw_mesh(junctionPolygonList, channelList,coordsDictionary):
	plt.axes()
	xmax = -1000000000000000.0
	ymax = -1000000000000000.0
	xmin = 1000000000000000.0
	ymin = 1000000000000000.0

	for junctionPolygon in junctionPolygonList: 
		points = []
		for vertex in junctionPolygon.vertices:
			points.append([vertex[0], vertex[1]])
			xmax = max(vertex[0], xmax)
			ymax = max(vertex[1], ymax)
			xmin = min(vertex[0], xmin)
			ymin = min(vertex[1], ymin)

		polygon = plt.Polygon(points, fill=None, edgecolor='r')
		plt.gca().add_patch(polygon)

	for channel in channelList:
		points = []
		for node in channel:
			points.append([coordsDictionary[node][0], coordsDictionary[node][1]])
		polygon = plt.Polygon(points, closed=None, fill=None, edgecolor='g')
		plt.gca().add_patch(polygon)

	plt.axis('scaled')
	#plt.savefig('Mesh.png')
	#plt.close()

	plt.show()

def generate_mesh(junctionPolygon):
	def round_trip_connect(start, end):
		return [(i, i+1) for i in range(start, end)] + [(end, start)]

	def connect_two_at_a_time(start, end):
		return[(i, i+1) for i in range(start,end+1,2)]

	vertices = []
	x = []
	y = []
	z = []
	n = []
	for vertex in junctionPolygon.vertices:
		xval = vertex[0]
		yval = vertex[1]
		zval = vertex[2]
		nval = vertex[3].item()
		newvertex = (xval, yval)
		vertices.append(newvertex)
		x.append(xval)
		y.append(yval)
		z.append(zval)
		n.append(nval)
		#print "xval = ", xval, " yval = ", yval, " zval = ", zval

	numVerts = len(vertices)	
	# interpolate the values to calculate the z value at the new
	# node(s) after meshing
	f = interp2d(x, y, z)
	fn = interp2d(x, y, n)

	#print "originalvertex = ", vertices 
	
	old_facets = round_trip_connect(0, junctionPolygon.numVertices-1)
	#channel_facets = connect_two_at_a_time(0,junctionPolygon.numVertices-1)
	
	channel_facets = junctionPolygon.channel_facets

	def needs_refinement(vertices, area):
		bary = np.sum(np.array(vertices), axis=0)/3
		max_area = 0.001 + abs(la.norm(bary, np.inf)-1)*0.01
		return bool(area > max_area)

	info = triangle.MeshInfo()
	info.set_points(vertices)
	info.set_facets(old_facets)

	#mesh = triangle.build(info, refinement_func=needs_refinement)
	mesh = triangle.build(info)
	mesh_points = np.array(mesh.points)
	mesh_tris = np.array(mesh.elements)
	mesh_facets = np.array(mesh.facets)
	newz = list(z)
	newnFriction = list(n)

	for j in xrange(numVerts, len(mesh_points[:,0])):
		newz.append(z[0])
		newnFriction.append(n[0])

	## calculate the interpolated value of z at the new nodes
	#for j in xrange(numVerts, len(mesh_points[:,0])):
	#	# find the original edge this point falls on
	#	my_x = mesh_points[j,0]
	#	my_y = mesh_points[j,1]
	#	found = False
	#	for facet in old_facets:
	#		v1 = facet[0]
	#		v2 = facet[1]
	#		x1 = vertices[v1][0]
	#		x2 = vertices[v2][0]
	#		y1 = vertices[v1][1]
	#		y2 = vertices[v2][1]
	#		z1 = z[v1]
	#		z2 = z[v2]
	#		n1 = n[v1]
	#		n2 = n[v2]
	#		m = (y2-y1)/(x2-x1)
	#		def is_on_line(x, y):
	#			if np.isfinite(m):
	#				return (abs(y - y1 - m*(x-x1)) < 1e-8)
	#			else:
	#				return (abs(x1-x) < 1e-8)
	#		if is_on_line(my_x, my_y):
	#			orig_edg_length = sqrt((x2-x1)**2+(y2-y1)**2)
	#			curr_edg_length = sqrt((my_x-x1)**2 + (my_y-y1)**2)
	#			t = curr_edg_length/orig_edg_length
	#			newz.append(t*z2+ (1-t)*z1)
	#			newnFriction.append(t*n2 + (1-t)*n1)
	#			found = True
	#			break
	#	
	#	if not found:
	#		newz.append(f(mesh_points[j,0], mesh_points[j,1])[0])
	#		newnFriction.append(fn(mesh_points[j,0], mesh_points[j,1])[0])

	# for the new junction edges created, find the channel connected to it
	def are_colinear(a,b,m,n,x,y):
		if ((n-b)*(x-m) == (y-n)*(m-a)):
			return True
		else:
			return False

	def are_equal(a,b):
		if (abs(a-b) < 1e-14):
			return True
		else:
			return False

	my_boundary = []
	#print "new_facets = ", mesh_facets
	for new_facet in mesh_facets:
		newN1 = new_facet[0]
		newN2 = new_facet[1]
		newx1 = mesh_points[newN1,0]
		newx2 = mesh_points[newN2,0]
		newy1 = mesh_points[newN1,1]
		newy2 = mesh_points[newN2,1]

		ch_facet_num = -1
		for facet in channel_facets:
			ch_facet_num += 1
			n1 = facet[0]
			n2 = facet[1]
			x1 = x[n1]
			x2 = x[n2]
			y1 = y[n1]
			y2 = y[n2]
			found = False
			
			# find the line connecting these two points
			m = (y2 - y1)/(x2 - x1)
			def is_on_desired_boundary(x, y):
				return (abs(y - y1 - m*(x -x1)) < 1e-8)
			if (is_on_desired_boundary(newx1, newy1) and is_on_desired_boundary(newx2, newy2)):
					my_boundary.append([newN1,newN2,junctionPolygon.boundary[ch_facet_num][0], junctionPolygon.boundary[ch_facet_num][1]])

	triangulation = tri.Triangulation(mesh_points[:,0], mesh_points[:,1], triangles=mesh_tris)
	fig = plt.figure()
	ax = plt.subplot(1,1,1)
	plt.gca().set_aspect('equal')
	ax = fig.gca(projection='3d')
	plotz = [-val for val in newz]
	ax.plot_trisurf(triangulation, plotz)
	#plt. triplot(mesh_points[:,0], mesh_points[:,1], mesh_tris)
	plt.show()

	#plt.figure()
	#plt.triplot(mesh_points[:,0], mesh_points[:,1], mesh_tris)
	#plt.show()

	return mesh_points, newz , newnFriction, mesh_tris, my_boundary


def main_func(gridFileName, channelNodesFile):
	channelList, channelDictionary = read_channel_nodes(channelNodesFile)
	junctionNodes = find_junction_nodes(channelList)
	channelList = break_channels(junctionNodes, channelList)
	coordsDictionary = read_coordinates(gridFileName)
	calculate_channel_bathymetry(channelList, channelDictionary, coordsDictionary)

	channelList, newToOldNode, junctionSkeletonList = create_junction_polygon_skeleton(channelList, junctionNodes, channelDictionary, coordsDictionary)

	junctionPolygonList = []
	numJunctions = 0
	for junctionSkeleton in junctionSkeletonList:
		junctionPolygon = create_junction_polygon(junctionSkeleton, coordsDictionary, channelDictionary)
		junctionPolygonList.append(junctionPolygon)
		numJunctions += 1

	#draw_mesh(junctionPolygonList, channelList,coordsDictionary)


	channelListWithCoords=[]
	for channel in channelList:
		newChannel =[]
		for node in channel:
			oldNode = newToOldNode.get(node)
			if (oldNode == None):
				oldNode = node
			my_coords= [coordsDictionary[node][0].item(), coordsDictionary[node][1].item(), channelDictionary[node][5].item(), channelDictionary[node][1].item(), channelDictionary[node][2].item(), channelDictionary[node][3].item(), channelDictionary[node][4].item(), channelDictionary[node][0].item(), oldNode]
			newChannel.append(my_coords)
		channelListWithCoords.append(newChannel)
	
	# write Channels.14 with new processed channels
	Cf = open('./Output/Channels.14','w')
	numChannels = len(channelListWithCoords)
	Cf.write('{}\n'.format(numChannels))

	for channel in channelListWithCoords:
		numNodes = len(channel)
		Cf.write('{}\n'.format(numNodes))
		for n in xrange(numNodes):
			my_vals = channel[n]
			Cf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(my_vals[0], my_vals[1], my_vals[2], my_vals[3], my_vals[4], my_vals[5]))
	Cf.close()

	# Create junction meshes and write them out to junction file
	f = open('./Output/JunctionMesh.14', 'w')
	f.write('Junction meshes generated with processMesh.py\n')
	f.write('{}\n'.format(numJunctions))
	myJunctionAttrList = []
	for junctionPolygon in junctionPolygonList:
			mesh_points, z, nFriction, mesh_tris, boundary = generate_mesh(junctionPolygon)
			myJunctionAttrList.append([mesh_points.tolist(),z, nFriction, mesh_tris.tolist(),boundary])
			numVerts = len(mesh_points)
			numTris = len(mesh_tris)
			writeString = "{}\t{}\n".format(numTris, numVerts)
			f.write(writeString)
			for j in xrange(numVerts):
				writeString = "{}\t{}\t{}\t{}\n".format(j+1, mesh_points[j,0], mesh_points[j,1], z[j])
				f.write(writeString)

			for j in xrange(numTris):
				writeString="{}\t{}\t{}\t{}\n".format(j+1, mesh_tris[j,0], mesh_tris[j,1], mesh_tris[j,2])
				f.write(writeString)

			connJuncEdges = len(boundary)
			f.write("Number of Junction Edges Connected to Channels = {}\n".format(connJuncEdges))
			#f.write("JunctionNode JunctionNode ConnectedChannelNumber ConnectionType\n");
			for boundary_edge in boundary:
				writeString="{} {} {} {}\n".format(boundary_edge[0], boundary_edge[1], boundary_edge[2], boundary_edge[3])
				f.write(writeString)

	f.close()

	returnList = [channelListWithCoords, myJunctionAttrList]

	return returnList

#main_func()
