""" Module for the integration of the TUM structural clustering algorithm
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies: gSpan, java
		
	Please cite the following article if you use the structural clustering procedure or results produced with it in any publication:
	
@InProceedings{ seeland2010,
	title = "Online Structural Graph Clustering Using Frequent Subgraph Mining",
	booktitle = "Proceedings of the ECML/PKDD'10",
	series = "Lecture Notes in Computer Science",
	author = "M. Seeland and T. Girschick and F. Buchwald and S. Kramer",
	editor = "J.L. Balc{\'a}zar and F. Bonchi and A. Gionis and M. Sebag",
	pages = "213--228",
	volume = "6323",
	year = "2010",
	booktitle1 = "Machine Learning and Knowledge Discovery in Databases, European Conference, {ECML} {PKDD} 2010, Barcelona, Spain, September 20-24, 2010, Proceedings, Part {III}"
}
"""

import os
import sys
import time
import shutil
import random

import subprocess
from subprocess import Popen, PIPE
import tempfile
from AZutilities import dataUtilities
import AZOrangeConfig as AZOC

def getStructuralClusters(data, threshold, minClusterSize, minClusterSaveSize = 0, minMolSize = 3, minSaveSDFsize = 0, numThreads=1, timeout=20):
	""" just the clustering
	    returns a list (of clusters) 
		of lists (which contain the smiles string of the cluster members)   
	"""
	clusters = []
	isSuccess = False
	tries = 0

	while (not isSuccess and tries < 10):
		tries += 1

		sdf_tempName = dataUtilities.makeTempSDF(data, smilesAsName=1)
		# create tempdir for usage as 6) outputpath
		temp_dir = tempfile.mkdtemp(prefix="AZorangeTMP_")

		# call clustering routine
		# Example command line call; gspan files are available in the same folder as the jar executable
		#java -jar structuralClustering.jar /home/girschic/proj/AZ/SAR/631.sdf 0.5 3 0 5 /home/girschic/proj/test/ . 2 20
		jarpath = os.path.join(AZOC.STRUCTCLUSTDIR,'structuralClustering.jar')
		opt = '-jar ' + jarpath + ' ' + sdf_tempName + ' ' + str(threshold) + ' ' + str(minMolSize) + ' ' + str(minClusterSaveSize) + ' ' + str(minClusterSize) + ' ' + temp_dir + '/ ' + str(AZOC.STRUCTCLUSTDIR) + ' ' + str(numThreads) + ' ' + str(timeout)
	
		cmd = 'java ' + opt
		p = Popen(cmd, shell=True, close_fds=True, stdout=PIPE)
		stdout = p.communicate()
		
		# parse output 
		outfile = os.path.join(temp_dir,'output_clusters.txt')
		try:
			if os.path.isfile(outfile):
				output = open(outfile, 'r')
				# 1,CCC(C)NC(=O)CSC1=NC2=C(C=CC(=C2)OCC)C=C1C#N	COC1=CC2=C(C=C1)N=C(C(=C2)C#N)SCC(=O)NC3=CC=C(C=C3)S(=O)(=O)N4CCCC4	CCC1=C(N=C2C=C3C(=CC2=C1)OCO3)SCC(=O)NC4=NOC(=C4)C	
				for line in output:
					tmp = line.strip()		
					split = tmp.partition(',')
					smilesList = split[2].split('\t')
					clusters.append(smilesList)
			else:
				print str(outfile) + " does not exist!"
				continue
		
		except IOError as (errno, strerror):
			print "I/O error({0}): {1}".format(errno, strerror) 					
			continue

		shutil.rmtree(temp_dir)
		isSuccess = True

	return clusters




def getReferenceStructures(data,threshold,minClusterSize,minClusterSaveSize = 0, minMolSize = 3, minSaveSDFsize = 0, numThreads=1, timeout=20, representativeMethod="random"):
	""" calls the clustering with the given parameters
		calls the different (now 1) methods to get representatives
	    returns a orange ExampleTable
	"""
	# get clusters
	clusters = getStructuralClusters(data, threshold, minClusterSize, minClusterSaveSize, minMolSize, minSaveSDFsize, numThreads, timeout)
	
	# get representatives
	representatives = []
	if (representativeMethod == 'random'):
		representatives = getRandomRepresentatives(clusters)
	else:
		print "Method for finding cluster representatives unknown!"
		return None
	
	return representatives


def getRandomRepresentatives(clusters):
	""" get random representative for each cluster
            returns list of smiles strings for usage in descriptor calculation
	"""
	reps = []
	for cluster in clusters:
		# generate random index for the cluster at hand
		end = len(cluster) - 1
		idx = random.randint(0,end)
		reps.append(cluster[idx])

	return reps

