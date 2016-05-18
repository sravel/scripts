#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
## @package struc2runClumpak.py
# @author Pierre Gladieux, Sebastien Ravel


##################################################
## Modules
##################################################
## Python modules
import sys, os, glob, argparse
from time import localtime, strftime
from matplotlib import rcParams, colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from MODULES_SEB import *

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='22/04/2016'

##################################################
## Functions

def relativeToAbsolutePath(relative):
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath

def lsExtInDirToList(pathDirectory, extentionFichierKeep):
	"""
	Return a list of 'ext' file's find in directory (exemple ext = "txt" or ["txt","py"])

	:param pathDirectory: a directory Path
	:type pathDirectory: Path
	:param extentionFichierKeep: a list or string with extention
	:type extentionFichierKeep: list or string
	:rtype: list()
	:return: list of 'ext' filename in pathDirectory ( file with extention find in param extentionFichierKeep )

	Example:
		>>> lsDirectory = lsExtInDirToList(path/to/directory/,"txt")
		>>> print(lsDirectory)
		["./out/gemo10_4497_ortho_rename_add.txt", "./out/gemo10_6825_ortho_rename_add.txt", "./out/gemo10_3497_ortho_rename_add.txt"]
	"""

	lsFilesFasta = []
	if pathDirectory[-1] != "/":
		pathDirectory += "/"
	if pathDirectory[-1] != "*":
		pathDirectory += "*"
	directoryFiles=glob.glob(pathDirectory)

	# Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
	for fichier in directoryFiles:
		try:
			if "." in fichier.split("/")[-1]:
				nameFichier = fichier.split("/")[-1].split(".")[0]
				extentionFichier = fichier.split("/")[-1].split(".")[1]
			else:
				nameFichier = fichier.split("/")[-1]
				extentionFichier = ""
		except:
			extentionFichier = "directory"
		if extentionFichier == extentionFichierKeep:
			lsFilesFasta.append(str(fichier))
	return sorted(lsFilesFasta)




##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='plot_barplot_DAPC.py', description='''This Programme Plot DAPC comboplot''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display plot_barplot_DAPC version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	files = parser.add_argument_group('Input info for running')
	files.add_argument('-d', '--directory', metavar="<path/to/directory>", required=True, dest = 'dirPath', help = 'path of result DAPC with csv file')
	files.add_argument('-l', '--label', metavar="<filename>", required=True, dest = 'labelFileParam', help = 'File with LABEL, one per line')
	files.add_argument('-o', '--output', metavar="<filename>", required=True, dest = 'outputFileParam', help = 'Name of output figure file')



	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print "#################################################################"
	print "#          Welcome in plot_barplot_DAPC (Version " + version + ")           #"
	print "#################################################################"
	print 'Start time: ', start_time,'\n'


	# Récupère le fichier de conf passer en argument
	dirPath = relativeToAbsolutePath(args.dirPath)
	labelFileParam = relativeToAbsolutePath(args.labelFileParam)
	outputFileParam = args.outputFileParam

	fileListCSV = lsExtInDirToList(dirPath, "csv")

	#read membership data and store in dictionary with keys=line numbers, which corresponds to the 'count' key in the tags dictionary
	memberships={}
	listK = []
	for file in fileListCSV:
		K=int((file.replace(dirPath+"/",'').replace('K','')).replace('.csv',''))
		if K not in listK:
			listK.append(K)
		IN=open(file,'r')
		IN.readline()		# remove header
		count=0
		for line in IN:
			count+=1
			line=line.strip().replace(',','.')
			tabLine=line.split(';')
			if count not in memberships.keys():
				memberships[count]={}
			memberships[count][K]=[round(float(x),2) for x in tabLine[1:]]
		IN.close()
	listKsorted = sorted(listK)
	minK = listKsorted[0]
	maxK = listKsorted[-1]


	#read strain names
	with open(labelFileParam,'r') as labelFile:
		IDs = [line.rstrip().split()[0] for line in labelFile.readlines()]

	#membership proportions
	Ks={}
	for K in range(minK,maxK+1):
		#print K
		Ks[K]={}
		for k in range(K):
			Ks[K][k]=[]
		for count in memberships:
			for k in range(K):
				Ks[K][k].append(memberships[count][K][k])


	ind=np.arange(len(IDs))
	width = 1.0
	#colors={2:{0:'deepskyblue',1:'forestgreen'},
			#3:{0:'firebrick',1:'forestgreen',2:'deepskyblue'},
			#4:{0:'orchid',1:'firebrick',2:'forestgreen',3:'deepskyblue'},
			#5:{0:'firebrick',1:'forestgreen',2:'deepskyblue',3:'orchid',4:'goldenrod'},
			#6:{0:'deepskyblue',1:'goldenrod',2:'blueviolet',3:'forestgreen',4:'orchid',5:'firebrick'},
			#7:{0:'blueviolet',1:'seagreen',2:'orchid',3:'deepskyblue',4:'goldenrod',5:'forestgreen',6:'firebrick'},
			#8:{0:'blueviolet',1:'forestgreen',2:'seagreen',3:'lightsteelblue',4:'orchid',5:'goldenrod',6:'firebrick',7:'deepskyblue'},
			#9:{0:'seagreen',1:'goldenrod',2:'forestgreen',3:'blueviolet',4:'deepskyblue',5:'lightsteelblue',6:'lightgreen',7:'orchid',8:'firebrick'},
			#10:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick',6:'blueviolet',7:'deepskyblue',9:'lightsteelblue',8:'darkred'}
			#}
	colors={2:{0:'goldenrod',1:'seagreen'},
			3:{0:'goldenrod',1:'seagreen',2:'lightgreen'},
			4:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid'},
			5:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen'},
			6:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick'},
			7:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick',6:'blueviolet'},
			8:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick',6:'blueviolet',7:'deepskyblue'},
			9:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick',6:'blueviolet',7:'deepskyblue',8:'lightsteelblue'},
			10:{0:'goldenrod',1:'seagreen',2:'lightgreen',3:'orchid',4:'forestgreen',5:'firebrick',6:'blueviolet',7:'deepskyblue',8:'lightsteelblue',9:'darkred'}
			}

	#smaller y fonts
	xlabel_size = 8
	ylabel_size = 8
	rcParams['xtick.labelsize'] = xlabel_size
	rcParams['ytick.labelsize'] = ylabel_size

	#lists -> stacked barplot
	#fig = plt.figure(figsize=(18, 15),dpi=1200,facecolor='w', edgecolor='w')
	fig = plt.figure(figsize=(40, 8),dpi=1200,facecolor='w', edgecolor='w')
	gs = gridspec.GridSpec(maxK-1,1)#, width_ratios=[ratios[0]/ratios[6],ratios[1]/ratios[6],ratios[2]/ratios[6],ratios[3]/ratios[6],ratios[4]/ratios[6],ratios[5]/ratios[6],1])
	#y ticks labels
	y = [0,0.5,1]
	#c is the gs index
	c=-1
	for K in range(maxK,1,-1):#loop in reverse order
		c+=1
		ax = plt.subplot(gs[c])
		bottom = np.zeros(len(IDs))
		for k in range(K):
			#ax.bar(ind, Ks[K][k], width,bottom=bottom, color=colors[K][k],edgecolor='white')#,align='center')
			ax.bar(ind, Ks[K][k], width,bottom=bottom, color=colors[K][k],edgecolor=(0, 0, 0, 0))#,align='center')
			bottom += Ks[K][k]
			if K>2:
				ax.xaxis.set_visible(False)
				ax.axis((0,len(Ks[K][k]),0,1))#ensures that I have 74 points, and not 80
				plt.yticks(y, list(y),rotation='vertical')
				ax.set_ylabel('K='+str(K))
				#outside y axis ticks
				plt.tick_params(
					axis='y',          # changes apply to y axis
					which='both',      # both major and minor ticks are affected
					left='on',
					direction='out',
					right='off'#right ticks are off
					)
			else:
				#ax.xaxis.set_visible(False)
				ax.axis((0,len(Ks[K][k]),0,1))#ensures that I have 74 points, and not 80
				ax.set_xticklabels(list(IDs))
				ax.set_ylabel('K='+str(K))
				plt.yticks(y, list(y),rotation='vertical')
				plt.xticks(ind+0.8, IDs, rotation='vertical')
				#outside y axis ticks
				plt.tick_params(
					axis='y',          # changes apply to y axis
					which='both',      # both major and minor ticks are affected
					left='on',
					direction='out',
					right='off'#right ticks are off
					)

	#remove ticks x axis
	plt.tick_params(
		axis='x',          # changes apply to x axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',         # ticks along the top edge are off
		) # labels along the bottom edge are off
	#to make room for the label
	plt.gcf().subplots_adjust(bottom=0.55)

	#plt.show()
	fig.savefig(dirPath+"/"+outputFileParam+".eps",format="eps")#,dpi=1200)
	#os.system("convert "+dirPath+"/"+outputFileParam+".eps "+dirPath+"/"+outputFileParam+".eps"+".pdf")
	fig.savefig(dirPath+"/"+outputFileParam+".pdf",format="pdf",dpi=1200)
