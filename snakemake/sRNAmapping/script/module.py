#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @package module.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
## Python modules
import argparse, pysam, traceback, sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict, OrderedDict
import pandas as pd
from tempfile import NamedTemporaryFile

# environment settings:
pd.set_option('display.max_column',None)
pd.set_option('display.max_rows',None)
pd.set_option('display.max_seq_items',None)
pd.set_option('display.max_colwidth', 500)
pd.set_option('expand_frame_repr', True)
# pd.options.display.width = None




#############################################
## Count library size for n to n size
#############################################
def getLibarySize(sample_path, size_min, size_max, csv_file):
	sampleName = Path(sample_path).name
	listFastq = sorted(Path(sample_path).resolve().rglob('*.fastq'))
	dicoLibrarySize = defaultdict(OrderedDict)
	allCount = 0
	for fastq in listFastq:
		for i in range(int(size_min),int(size_max)+1,1):
			sizeName = f"size{i}"
			if sizeName in fastq.name:
				n = 0
				with open(fastq, "r") as fd:
					lu = fd.read()
					n = lu.count('\n')
					allCount += n
				dicoLibrarySize[sampleName][f"{i:.0f}"] = f"{n/4:.0f}"
	dicoLibrarySize[sampleName][f"{size_min}-{size_max}"] = f"{allCount/4:.0f}"
	dataframeLibrarySize = pd.DataFrame.from_dict(dicoLibrarySize, orient='index')
	with open(csv_file, "w") as libsizeFile:
		print(f"Library size:\n{dataframeLibrarySize}\n")
		dataframeLibrarySize.to_csv(libsizeFile, index=True)

def mergeCSV(csv_files, csv_file, sep=","):
	# dir = Path(csv_files)
	df = (pd.read_csv(f, sep=sep) for f in csv_files)
	df = pd.concat(df)
	df.rename(columns={'Unnamed: 0':'Samples'}, inplace=True)
	with open(csv_file, "w") as libsizeFile:
		print(f"All CSV infos:\n{df}\n")
		df.to_csv(libsizeFile, index=False, sep=sep)

def filterBam(bam_file, size_min, size_max, bam_filter_file, sample):
	if not Path(bam_file+".bai").exists(): pysam.index(bam_file)		# index le bam si besoin
	# Création du fichier des readgroups à garder:
	with NamedTemporaryFile(mode='w',  delete=False) as fp:
	# with open(Path(bam_filter_file).parent.joinpath("readGroup.txt"), mode='w') as fp:
		# print(fp.name)
		for i in range(int(size_min),int(size_max)+1,1):
			print(f"{sample}-size{i}")
			fp.write(f"{sample}-size{i}\n")
	bam = pysam.view("-b", "-h","-R", fp.name, "-o", bam_filter_file, bam_file, catch_stdout=False)	## utilisation de pysam view généré un bam avec les reads groups garder au dessus
	pysam.index(bam_filter_file)

def plotGraph(csv_file_name, sample_name, reference_name, out_file_name, sep = "\t", size_keep_min = 20, size_keep_max = 25, font_size = 18, chunk_size = 80, same_y_scale = False):
	import pandas as pd
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	import numpy as np
	plt.style.use('ggplot')
	mpl.use('Agg')

	df_count_pos_normalized = pd.read_csv(csv_file_name, sep=sep, header=0 ) #
	list_size_keep = range(size_keep_min, size_keep_max+2)

	# initialise the figure. here we share X and Y axis
	fig, axes = plt.subplots(figsize=(45, 35),dpi=150,nrows=len(list_size_keep), ncols=1, sharex=True, sharey=same_y_scale)
	plt.subplots_adjust(top=0.95)
	fig.suptitle(f'Graph of mapping sRNA from {sample_name} to {reference_name}', fontsize=font_size+6)

	# boucle sur les taille 20-25
	for size, i in zip(list_size_keep,range(len(list_size_keep))):
		if i < len(list_size_keep)-1:
			# print(f"plot figure\tsize={size}\ti={i}")
			select_forward = df_count_pos_normalized[f"Forward {size}nt"].to_dict()
			select_reverse = df_count_pos_normalized[f"Reverse {size}nt"].to_dict()
			maxFor = max(select_forward.values())
			minRev = min(select_reverse.values())
			print(f"plot figure\tsize={size}\ti={i}\tmin={minRev}\tmax={maxFor}")
			axes[i].margins(0, 0)																# marge a 0 pour l'intérieur du graph
			ind = np.arange(len(select_reverse))														# axe des x
			width = 10																			# the width of the bars: can also be len(x) sequence
			axes[i].bar(ind, select_forward.values(), width, color='b',label=f"Forward")	# plot du forward
			axes[i].bar(ind, select_reverse.values(), width, color="r",label=f"Reverse")	# plot du reverse
			axes[i].set_ylabel('Expressed in rpm', fontsize=font_size)								# ajoute le label axe y sur chaque plot
			axes[i].set_title(f"{size}-nt viral sRNAs", fontsize=font_size)
			axes[i].set_xlim(1,len(select_reverse))													# Limite entre 0 et n l'axe des x (supprime la marge)
			axes[i].set_ylim(minRev-3,maxFor+3)													# Limite entre 0 et n l'axe des y (supprime la marge)
			axes[i].legend(ncol=2, bbox_to_anchor=(0, 1), loc='lower left', fontsize='medium')	# ajoute la légende loc="best",

		else:
			# print(f"plot figure\tTOTAL\ti={i}")
			select_forward = df_count_pos_normalized[f"Total Forward"].to_dict()
			select_reverse = df_count_pos_normalized[f"Total Reverse"].to_dict()
			maxFor = max(select_forward.values())
			minRev = min(select_reverse.values())
			print(f"plot figure\tTOTAL\ti={i}\tmin={minRev}\tmax={maxFor}")
			axes[i].margins(0, 0)																# marge a 0 pour l'intérieur du graph
			ind = np.arange(len(select_reverse))														# axe des x
			width = 10																			# the width of the bars: can also be len(x) sequence
			axes[i].bar(ind, select_forward.values(), width, color='b',label=f"Forward")	# plot du forward
			axes[i].bar(ind, select_reverse.values(), width, color="r",label=f"Reverse")	# plot du reverse
			axes[i].set_ylabel('Expressed in rpm',fontsize=font_size)													# ajoute le label axe y sur chaque plot
			axes[i].set_xlabel(f'{reference_name} genome size',fontsize=font_size)													# ajoute le label axe x sur chaque plot
			axes[i].set_title(f"{size_keep_min}-{size_keep_max}-nt viral sRNAs",fontsize=font_size)
			axes[i].set_xlim(1,len(select_reverse))													# Limite entre 0 et n l'axe des x (supprime la marge)
			axes[i].set_ylim(minRev-3,maxFor+3)														# Limite entre 0 et n l'axe des y (supprime la marge)
			axes[i].legend(ncol=2, bbox_to_anchor=(0, 1), loc='lower left', fontsize='medium')													# ajoute la légende loc="best",

	plt.xticks(np.arange(1, len(select_reverse),  int(len(select_reverse)/chunk_size)), rotation=70)				# ajoute la grille an 80 partie
	plt.savefig(out_file_name, bbox_inches='tight',pad_inches=1)

def mergeImages(img1, img2, output_name):
	import numpy as np
	from PIL import Image
	from time import sleep
	im1 = Image.open(img1)
	im2 = Image.open(img2)
	imgs_comb = Image.new('RGB', (im1.width + im2.width, min(im1.height, im2.height)))
	imgs_comb.paste(im1, (0, 0))
	imgs_comb.paste(im2, (im1.width, 0))
	imgs_comb.save(output_name)
	sleep(50)
