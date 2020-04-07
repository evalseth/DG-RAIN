import numpy as np
from math import sqrt
from math import atan
from math import atan2
from math import cos
from math import sin

class junction_skeleton:
	def __init__(self, centerNode, centerNodeCoords):
		self.centerNode = centerNode
		self.centerNodeCoords = centerNodeCoords
		self.definingNodes = []
		self.connectedChannel = []
		self.connectedChannelBoundaryType = []
		
	def add_connectedChannel(self, connectedChannel):
		self.connectedChannel.append(connectedChannel)

	def add_connectedChannelBoundaryType(self, connectedChannelNodeNumber):
		self.connectedChannelBoundaryType.append(connectedChannelNodeNumber)
	
	def add_definingNode(self, definingNode):
		self.definingNodes.append(definingNode)
				
def create_junction_polygon_skeleton(channelList, junctionNodes, channelDictionary, coordsDictionary):
	r"""This function creates the skeleton of the 2-D junction polygon, which 
	consists of the center node (the junction node)
	and the channel nodes it is connected to."""
	junctionSkeletonList = []
	newToOldNode = {}
	numChannels = len(channelList)
	numCoordinates = len(coordsDictionary)
	for node in junctionNodes:
		centerNodeCoords = coordsDictionary[node]
		junctionSkeleton = junction_skeleton(node, centerNodeCoords)
		for i in xrange(numChannels):
			channel = channelList[i]
			if node in channel:
				# width at the center node of the junction
				width = channelDictionary[node][1]

				# introduce a new node that will be half the width away
				# from the junction node
				newNode = numCoordinates + 1
				numCoordinates = newNode
				
				nodeLocation = np.where(channel==node)[0][0]
				channelLength = len(channel)

				# add the channels and the connected nodes to the junction
				#skeleton. Swap the junction node from the channels
				# for a node that is one width away from the junctionNode
				junctionSkeleton.add_definingNode(newNode)
				junctionSkeleton.add_connectedChannel(i)
				
				newToOldNode[newNode] = channelList[i][nodeLocation]
				channelList[i][nodeLocation]= newNode

				# keep track of the channel is going out from this junction (2)
				# or if it is going into the junction (1)
				nextNode = 0
				if nodeLocation == 0:
					junctionSkeleton.add_connectedChannelBoundaryType(2)
					nextNode = channel[1]
				else:
					junctionSkeleton.add_connectedChannelBoundaryType(1)
					nextNode = channel[channelLength-2]
		
				# Calculate the coordinates of newNode (linearly interpolate x, y, z, b, m1, m2)
				x0 = coordsDictionary[node][0]
				y0 = coordsDictionary[node][1]
				z0 = channelDictionary[node][5]
				b0 = channelDictionary[node][1]
				m10 = channelDictionary[node][2]
				m20 = channelDictionary[node][3]
				nf0 = channelDictionary[node][4]
				d0  = channelDictionary[node][0]
				
				x1 = coordsDictionary[nextNode][0]
				y1 = coordsDictionary[nextNode][1]
				z1 = channelDictionary[nextNode][5]
				b1 = channelDictionary[nextNode][1]
				m11 = channelDictionary[nextNode][2]
				m21 = channelDictionary[nextNode][3]
				nf1 = channelDictionary[nextNode][4]
				d1 = channelDictionary[nextNode][0]

				edg_length = sqrt((x0-x1)**2 + (y0-y1)**2)

				# assume that 2*width is less than edg_length
				t = 6*width/edg_length

				newx = x1*t+(1-t)*x0;
				newx = round(newx, 2) 
				newy = y1*t+(1-t)*y0;
				newy = round(newy, 2)
				newz = z1*t+(1-t)*z0;
				newz = round(newz, 2)
				newb = b1*t+(1-t)*b0;
				newb = round(newb, 2)
				newm1 = m11*t+(1-t)*m10;
				newm1 = round(newm1, 2)
				newm2 = m21*t+(1-t)*m20;
				newm2 = round(newm2,2)
				newnFriction = nf1*t + (1-t)*nf0;
				newnFriction = round(newnFriction, 2)
				newd1 = d1*t + (1-t)*d0;
				newd1 = round(newd1, 2)
			
				# for now take d to be the depth at the center node
				newb = b0
				newz = z0
				newd1 = d0
				
				coordsDictionary[newNode] = np.array([newx, newy, newz-newd1, newnFriction])
				channelDictionary[newNode] = np.array([newd1, newb, newm1, newm2, newnFriction,newz])

		junctionSkeletonList.append(junctionSkeleton)
	
	return channelList, newToOldNode, junctionSkeletonList

def sort_counterclockwise_order(center, vertices, coordsDictionary, connectedChannel):
	vertDict = {}
	connChanDict = {}
	xcenter = coordsDictionary[center][0]
	ycenter = coordsDictionary[center][1]
	angleList = []
	i = 0
	for vertex in vertices:
		x = coordsDictionary[vertex][0]
		y = coordsDictionary[vertex][1]
		angle = atan2(y-ycenter,x-xcenter)
		angleList.append(angle)
		vertDict[angle] = vertex
		connChanDict[angle] = connectedChannel[i]
		i += 1

	angleList = sorted(angleList)
	sortedVertices = []
	sortedChannels = []
	for key in angleList:
		sortedVertices.append(vertDict[key])
		sortedChannels.append(connChanDict[key])
	
	return sortedVertices, sortedChannels

def find_perp_nodes(center, vertices, coordsDictionary, channelDictionary):
	perpNodesX = {}
	perpNodesY = {}
	xcenter = coordsDictionary[center][0]
	ycenter = coordsDictionary[center][1]
	#alpha = atan(0.25)		
	#print "center = ", xcenter, "\t", ycenter
	for vertex in vertices:
		x = coordsDictionary[vertex][0]
		y = coordsDictionary[vertex][1]
		b = channelDictionary[vertex][1]
		dh_sq = (x-xcenter)**2 + (y-ycenter)**2
		h = sqrt(dh_sq+(0.5*b)**2)
		alpha = atan(0.5*b/sqrt(dh_sq))
		#h = 1.0
		theta = atan2(y-ycenter, x-xcenter)
		theta1 = theta - alpha
		theta2 = theta + alpha
		perpNodesX[vertex] = [round(xcenter + h*cos(theta1),2), round(xcenter + h*cos(theta2),2)]
		perpNodesY[vertex] = [round(ycenter + h*sin(theta1),2), round(ycenter + h*sin(theta2),2)]
		#print "node = ", x, "\t", y
		#print "xcenter = ", xcenter, " ycenter = ", ycenter
		#print "perpNodesX = ", perpNodesX[vertex]
		#print "perpNodesY = ", perpNodesY[vertex]
	
	return perpNodesX, perpNodesY

def find_intersecting_node(x1, y1, x2, y2, x11, y11, x21, y21, xc, yc):
	tol = 0.01
	if (abs(x1-xc) > tol and abs(x2-xc) > tol):
		xint = (y11 - y21 - x11*(y1-yc)/(x1-xc)+x21*(y2-yc)/(x2-xc))/((y2-yc)/(x2-xc) - (y1-yc)/(x1-xc))
		yint = y11 - x11*(y1-yc)/(x1-xc) +xint*(y1-yc)/(x1-xc)

	elif (abs(x1 - xc) < tol and abs(x2 - xc) > tol):
		xint = x11
		#print x1, " ", y1, " ", x11, " " , y11, " ", x2, " ", y2, x21, " ", y21, " " , xc, " ", yc
		yint = y21 + (yc-y2)/(xc-x2)*(x11-x21)
	elif (abs(x1 - xc) > tol and abs(x2 - xc) < tol):
		xint = x21
		#print x1, " ", y1, " ", x11, " " , y11, " ", x21, " ", y21
		yint = y11 + (yc - y1)/(xc - x1)*(x21 - x11)
	
	else: # both zero
		xint = None
		yint = None
	
	if xint and yint:
		xint = round(xint, 2)
		yint = round(yint,2)

	return xint, yint

class junction_polygon:
	def __init__(self, junctionSkeleton, coordsDictionary):
		# sort the skeletal vertices in counter-clockwise order
		definingNodes = junctionSkeleton.definingNodes
		#self.connectedChannel = junctionSkeleton.connectedChannel
		self.connectedChannelBoundaryType = junctionSkeleton.connectedChannelBoundaryType
		self.centerNode = junctionSkeleton.centerNode
		self.skeletalNodes, self.connectedChannel = sort_counterclockwise_order(self.centerNode, definingNodes, coordsDictionary, junctionSkeleton.connectedChannel)
		xc = coordsDictionary[self.centerNode][0]
		yc = coordsDictionary[self.centerNode][1]
		#print "xc = ", xc, " yc = ", yc
		#for node in self.skeletalNodes:
		#	print coordsDictionary[node][0], " ", coordsDictionary[node][1]
		self.vertices = []
		self.boundary = []
		self.channel_facets = []
		self.numVertices = 0

	def add_polygon_vertices(self,coordsDictionary, channelDictionary):
		xVerticesDict, yVerticesDict = find_perp_nodes(self.centerNode, self.skeletalNodes, coordsDictionary, channelDictionary)
		numSkeletalVertices = len(self.skeletalNodes)
		vertexNum = 0
		for i in xrange(numSkeletalVertices):
			curr = i
			currSkeletalNode = self.skeletalNodes[curr]
			x1 = xVerticesDict[currSkeletalNode][0]
			y1 = yVerticesDict[currSkeletalNode][0]
			x2 = xVerticesDict[currSkeletalNode][1]
			y2 = yVerticesDict[currSkeletalNode][1]
			z = channelDictionary[currSkeletalNode][5]
			n = channelDictionary[currSkeletalNode][4]
		
			self.vertices.append([x1, y1, z, n])
			self.vertices.append([x2, y2, z, n])
			self.channel_facets.append([vertexNum, vertexNum+1])
			vertexNum += 2
			self.boundary.append([self.connectedChannel[curr], self.connectedChannelBoundaryType[curr]])
			#print "x1 = " ,x1, " y1 = ", y1, " x2 = ", x2, " y2 = ", y2, " z = ", z, " connChan = ", self.connectedChannel[curr]	
		self.numVertices = len(self.vertices)

	def add_polygon_vertices_with_intersection_option(self,coordsDictionary, channelDictionary):
		xVerticesDict, yVerticesDict = find_perp_nodes(self.centerNode, self.skeletalNodes, coordsDictionary, channelDictionary)
		numSkeletalVertices = len(self.skeletalNodes)
		xc = coordsDictionary[self.centerNode][0]
		yc = coordsDictionary[self.centerNode][1]
		zc = channelDictionary[self.centerNode][5]
		vertexNum = 0
		for i in xrange(numSkeletalVertices):
			curr = i
			next_vrt = (curr+1) % numSkeletalVertices
			currSkeletalNode = self.skeletalNodes[curr]
			nextSkeletalNode = self.skeletalNodes[next_vrt]
			x1 = coordsDictionary[currSkeletalNode][0]
			y1 = coordsDictionary[currSkeletalNode][1]
			z1 = channelDictionary[currSkeletalNode][5]
			n1 = channelDictionary[currSkeletalNode][4]
			x2 = coordsDictionary[nextSkeletalNode][0]
			y2 = coordsDictionary[nextSkeletalNode][1]
			z2 = channelDictionary[nextSkeletalNode][5]
			n2 = channelDictionary[currSkeletalNode][4]
			
			x10 = xVerticesDict[currSkeletalNode][0]
			y10 = yVerticesDict[currSkeletalNode][0]
			x11 = xVerticesDict[currSkeletalNode][1]
			y11 = yVerticesDict[currSkeletalNode][1]
			x20 = xVerticesDict[nextSkeletalNode][0]
			y20 = yVerticesDict[nextSkeletalNode][0]
			x21 = xVerticesDict[nextSkeletalNode][1]
			y21 = yVerticesDict[nextSkeletalNode][1]

			self.vertices.append([x11, y11, z1, n1])
			self.channel_facets.append([vertexNum-1, vertexNum])
			self.boundary.append([self.connectedChannel[curr], self.connectedChannelBoundaryType[curr]])
			vertexNum += 1
			
			xint, yint = find_intersecting_node(x1, y1, x2, y2, x11, y11, x20, y20, xc, yc)
			if (xint is not  None and yint is not None):
				#print "intersection nodes"
				#print xc, " ", yc
				#print x1, " ", y1
				#print x2, " ", y2
				#print x11, " ", y11
				#print x21, " ", y21
				#print xint, " ", yint
				self.vertices.append([xint, yint, zc, 0.5*(n1+n2)])
				vertexNum +=1
				#self.boundary.append([-1,-1])
		
			self.vertices.append([x20, y20,z2, n2])
			vertexNum += 1
	
		self.channel_facets[0][0] = vertexNum-1
		self.numVertices = len(self.vertices)

def create_junction_polygon(junctionSkeleton, coordsDictionary, channelDictionary):
	junctionPolygon = junction_polygon(junctionSkeleton, coordsDictionary)
	#junctionPolygon.add_polygon_vertices(coordsDictionary, channelDictionary)

	junctionPolygon.add_polygon_vertices_with_intersection_option(coordsDictionary, channelDictionary)
	return junctionPolygon



