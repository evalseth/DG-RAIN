import numpy as np
from math import sqrt

def read_channel_nodes(filename):
	"This function reads in the channel nodes and returns a list\
	containing the channels represented by their nodes."
	channelList = []
	dataDictionary = {}
	channelNodesFile = open(filename,"r")
	throwawayline = channelNodesFile.readline()
	numChannels = int(channelNodesFile.readline())
	
	for i in range(0, numChannels):
		numNodes = int(channelNodesFile.readline())
		newChannel = []
		for j in range(0, numNodes):
			node, depth, b, m1, m2, nfriction = [float(value) for value in channelNodesFile.readline().split()]
			node = int(node)
			node = node-1
			newChannel.append(node)
			dataArray = np.array([depth, b, m1, m2,nfriction,0.0])
			dataDictionary[node] = dataArray

		channelList.append(np.array(newChannel))
	
	channelNodesFile.close()
	return channelList, dataDictionary

def read_coordinates(filename):
	r""" This function reads the coordinates of the nodes"""

	grid_file = open(filename, 'r')
	grid_file.readline()
	num_elements, num_nodes = [int(value) for value in grid_file.readline().split()]
	
	# Create necessary array storage
	coords =	{} 
	
	# Read in coordinates of each node
	for n in xrange(num_nodes):
	 	line = grid_file.readline().split()
		dataArray = np.array([round(float(line[1]),2), round(float(line[2]),2), round(float(line[3]),2), round(float(line[4]),2)])
	 	coords[int(line[0])-1] = dataArray
	
	grid_file.close()

	return coords

def calculate_channel_bathymetry(channelList, channelDictionary, coordsDictionary):
	"This function calculates the bathymetry from the depth information provided \
	with channels and the topography information of adjacent land"
	numChannels = len(channelList)
	for i in range(0, numChannels):
		channel = channelList[i]
		for node in channel:
			channelDictionary[node][5]=channelDictionary[node][0].item()+coordsDictionary[node][2].item()

def find_junction_nodes(channelList):
	"This function finds the junction nodes"
	numChannels = len(channelList)
	junctionNodes = []

	# compare all channels with each other to find junction nodes
	for i in range(0,numChannels):
		for j in range(i+1,numChannels):
			arr1 = channelList[i]
			arr2 = channelList[j]
			commonNodes = np.intersect1d(arr1, arr2)
			for node in commonNodes:
				# add the node to the junction nodes list if it
				#doesn't already exist
				if node not in junctionNodes:
					junctionNodes.append(node)	
	return junctionNodes

def break_channels(junctionNodes, oldChannelList):
	"This function scans the channel nodes to ensure that a \
	junction doesn't occur in the middle of a channel. If it\
	does, then the function breaks the channels into smaller\
	channels"
	
	channelList = []
	junctionNodes = np.array(junctionNodes)
	# first make sure that a channel is connected to a junction
	# only at the end points
	for channel in oldChannelList:
		# find out how many junction nodes the channel contains
		indices = []
		intersection = np.intersect1d(junctionNodes, channel)
		num_junc_nodes = intersection.size
		if num_junc_nodes :
			for node in intersection:
				indices.append(np.where(channel==node)[0][0])
			
			indices.sort()
			channelLength = len(channel)
			for i in range(0,len(indices)+1):
				if i == 0:
					firstIndex = 0;
				else:
					firstIndex = indices[i-1]
				
				if i == len(indices):
					lastIndex = channelLength
				else:
					lastIndex = indices[i]+1

				subChannel1 = (channel[firstIndex:lastIndex]).tolist()
				#subChannel2 = (channel[indices[i]:channelLength]).tolist()

				if len(subChannel1) > 1:
					channelList.append(subChannel1)

				#if num_junc_nodes == 1:
				#	if len(subChannel2) > 1:
				#		channelList.append(subChannel2)

		else:
			channelList.append(channel)

			
	return channelList




