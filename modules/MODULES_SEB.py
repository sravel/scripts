#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
## @package MODULES_SEB.py
# @author Sebastien Ravel
#__docformat__ = "restructuredtext en"
"""
	The MODULES_SEB module
	======================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 07-06-2016
	:version: 0.2

	Use it to import very handy functions.

	Example:

	>>> from MODULES_SEB import dict2txt
	>>> dico = {"key1":"value1","key2":"value2","key3":"value3"}
	>>> dict2txt(dico)
	key1	value1
	key2	value2
	key3	value3

	Required Module install
	-----------------------

	This module run with Python 3.x and not Python 2.x

	** Include module of python:
		- argparse, os, subprocess, sys, time, glob
	** Module need to install before use:
		- BioPython
		- pyvcf
		- pyfaidx

"""
##################################################
## Modules
##################################################
## Python modules
import argparse, os, subprocess, sys
from time import localtime, strftime, sleep, clock, time
import glob, re
## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet, IUPAC, Gapped
from Bio import AlignIO
from Bio.Nexus import Nexus
from Bio import AlignIO

##ParseGFF modules
from collections import namedtuple
import gzip
import urllib


##################################################
## Variables Globales
version="0.1"
VERSION_DATE='07-06-2016'

##################################################
## Fonctions

def existant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open by default.

	:param x: a file path
	:type x: str()
	:rtype: string
	:return: string

	"""
	if not os.path.exists(x):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))

	return x


def _test():
	"""
	Function to launch doctests use with _test ()
	:note: You need install doctest module.
	"""
	import doctest
	doctest.testmod()


def sort_human(s, _nsre=re.compile('([0-9]+)')):
	""" Sort the list in the way that humans expect, use list.sort(key=sort_human) or sorted(list, key=sort_human)).

	:param s: a python list
	:type s: list()
	:rtype: list()
	:return: liste human sort

	Example:
		>>> listToSorted = ["something1","something32","something17","something2","something29","something24"]
		>>> print(listToSorted.sort(key=sort_human))
		['something1', 'something17', 'something2', 'something25', 'something29', 'something32']
		>>> print(sorted(listToSorted, key=sort_human))
		['something1', 'something17', 'something2', 'something25', 'something29', 'something32']

	"""
	try:
		return [ int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]
	except TypeError:
		if not isinstance(s,int):
			print("WARNNING MODULES_SEB::sort_human : List %s value not understand so don't sort \n" % s)
		return s



def dict2txt(dico):
	"""
	Function that takes a dictionary and returns a tabular string with::

		"key\\tvalue\\n".

	:param dico: a python dictionary
	:type dico: dict()
	:rtype: str()
	:return: string with "key\\tvalue\\n

	Example:
		>>> dico = {"key1":"value1","key2":"value2","key3":"value3"}
		>>> dict2txt(dico)
		key1	value1
		key2	value2
		key3	value3

	:warn: if the value of the value list or dictionary, the display will be plain and without formatting data.
	"""

	txtoutput = ""
	for key in sorted(dico.keys(), key=sort_human):
		value = dico[key]
		txtoutput += "%s\t%s\n" % (str(key),str(value))
	return txtoutput

def dictList2txt(dico):
	"""
	Function that takes a dictionary and returns a tabular string with::

		"key\\tvalue\\n".

	:param dico: a python dictionary
	:type dico: dict()
	:rtype: str()
	:return: string with "key\\tvalue\\n

	Example:
		>>> dico = {"key1":["value1","value1"], "key2":["value2","value2"],"key3":["value3","value3"]}
		>>> dict2txt(dico)
		key1	value1	value1
		key2	value2	value2
		key3	value3	value3

	"""

	txtoutput = ""
	for key in sorted(dico.keys(), key=sort_human):
		value = "\t".join(dico[key])
		txtoutput += "%s\t%s\n" % (str(key),str(value))
	return txtoutput

def dictDict2txt(dico,first="Info"):
	"""
	Function that takes a dictionary and returns a tabular string with::

		"key\\tvalue\\n".

	:param dico: a python dictionary
	:type dico: dict()
	:param first: string for name first column
	:type first: str()
	:rtype: str()
	:return: string with "key\\tvalue\\n

	Example:
		>>> dico = {"Souche1":{"NUM":"171","MIN":"2042","MAX":"3133578","N50 BP":"938544","N50 NUM":"11"},
					"Souche2":{"NUM":"182","MIN":"5004","MAX":"74254","N50 BP":"45245","N50 NUM":"45"}}
		>>> dictDict2txt(dico,"souches")
		souches	NUM	MIN	MAX	N50 BP	N50 NUM
		Souche1	171	2042	3133578	938544	11
		Souche2	182	5004	74254	45245	45

	"""

	txtoutput = ""
	headerc=0
	for key in sorted(dico.keys(), key=sort_human):
		dicoInfosValues = list(dico[key].values())
		dicoInfosValues = [ str(k) for k in dicoInfosValues]

		if headerc == 0:
			header = [ str(k) for k in sorted(dico[key].keys(), key=sort_human) ]
			value = first+"\t" + "\t".join(header)
			txtoutput += "%s\n" % str(value)
			headerc=1

		value = "\t".join([ str(dico[key][key2]) for key2 in header])
		txtoutput += "%s\t%s\n" % (str(key),str(value))
	return txtoutput

def dict2fasta(dico):
	"""
	Function that takes a dictionary with key are ID and value Seq, and returns a fasta string.

	:param dico: a python dictionary
	:type dico: dict()
	:rtype: fasta str()
	:return: string with format fasta

	Example:
		>>> dico = {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
		>>> dict2fasta(dico)
		>Seq1
		ATGCTGCAGTAG
		>Seq2
		ATGCCGATCGATG
		>Seq3
		ATGCTCAGTCAGTAG
	"""

	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">%s\n%s\n" % (str(key),str(value[0]))
	return txtoutput

def fasta2dict(filename):
	"""
	Function that take a file name (fasta), and return a dictionnary of sequence

	:param filename: a fasta file
	:type filename: file in fasta format
	:rtype: record_dict()
	:return: dict() - dictionnary with keys are Id and value SeqRecord() fields
	:requires: this function require ## BIO Python modules: (from Bio import SeqIO,\\n
	from Bio.SeqRecord import SeqRecord \\n
	from Bio.Seq import Seq \\n
	from Bio.Alphabet import SingleLetterAlphabet)

	Example:
		>>> filename = "sequence.fasta"
		>>> fasta2dict(filename)
		{">Seq1":"SeqRecord()"}
	"""

	# chargement du fasta des MGG en mémoire
	handle = open(filename, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()
	return record_dict

def lenSeq2dict(filename):
	"""
	Function that take a file name (fasta), and return a dictionnary with length of sequence

	:param filename: a fasta file
	:type filename: file in fasta format
	:rtype: record_dict()
	:return: dict() - contain length of sequences
	:requires: this function require fasta2dict(filename)

	Example:
		>>> filename = "sequence.fasta"
		>>> lenSeq2dict(filename)
		{">Seq1":20154}
	"""

	dicoLenMGG = {}
	record_dict = fasta2dict(filename)
	for gene in sorted(record_dict.keys(), key=sort_human):
		if record_dict[gene].id not in dicoLenMGG:
			lenseq = len(record_dict[gene].seq)
			dicoLenMGG[gene]=int(lenseq)
	return dicoLenMGG

def nbSeqInFile2dict(pathDirectory):
	"""
	Function  that take a Path Directory and returna dictionnary with number of sequences in fasta file's

	:param pathDirectory: a directory Path
	:type pathDirectory: Path
	:rtype: dict1(), dict2()
	:return: - contient le nombre de sequences dans les fichiers (key = nom de fichier value =  nombre de sequences)\n
			 - contient le nombre de fichier qui ont x sequences (key = nombre de sequence =  nombre de fichier)
	:raise print: print("ERROR: Sequence: "+nameFichier+" allready read") with nameFichier is the current file read.

	Example:
		>>> dico1,dico2 = nbSeqInFile2dict(path/to/directory/)
		>>> print(dict2txt(dico1))
		./out/gemo10_4497_ortho_rename_add.fasta	58
		./out/gemo10_6825_ortho_rename_add.fasta	59
		./out/gemo10_3497_ortho_rename_add.fasta	59
		./out/gemo10_6254_ortho_rename_add.fasta	59
		>>> print(dict2txt(dico2))
		58	1
		59	3
	"""

	dicoNbSeqInFiles = {}
	dicoNbFilesNbSouche = {}
	# récupération des fichiers du repertoire
	if pathDirectory[-1] != "*":
		pathDirectory += "*"
	directoryFiles=glob.glob(pathDirectory)
	# Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
	for fichier in directoryFiles:
		try:
			nameFichier = fichier.split("/")[-1].split(".")[0]
			extentionFichier = fichier.split("/")[-1].split(".")[1]
		except:
			extentionFichier = "directory"
		if extentionFichier in ["fasta", "fa", "fas"]:
			# Ouverture des sequences fasta et chargement dans dictionnaire
			handle = open(fichier, "r")
			record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
			handle.close()
			nbseq=len(record_dict.keys())
			#print(nameFichier+"\t"+str(nbseq))
			if nameFichier not in dicoNbSeqInFiles.keys():
				dicoNbSeqInFiles[nameFichier]=nbseq
			else:
				print("ERROR: Sequence: "+nameFichier+" allready read")

			dicoNbFilesNbSouche[nbseq] = (dicoNbFilesNbSouche.get(nbseq, 1))+1

			#if nbseq not in dicoNbFilesNbSouche.keys():
				#dicoNbFilesNbSouche[nbseq] = 1
			#else:
				#dicoNbFilesNbSouche[nbseq] += 1

	return dicoNbSeqInFiles,dicoNbFilesNbSouche

def readable_dir(prospective_dir):
	"""
	Check if directory exist and if is readable

	:param prospective_dir: a directory Path
	:type prospective_dir: Path

	:raise error: raise argparse.ArgumentTypeError(" :{0} is not a valid path".format(prospective_dir))
	:raise error: raise argparse.ArgumentTypeError(" :{0} is not a readable dir".format(prospective_dir))

	Example:
		>>> parser = argparse.ArgumentParser(prog='make_structure_dir.py', description='''This Programme make arborescence of rep of programme structure''')
		>>> paths = parser.add_argument_group('Input PATH for running')
		>>> paths.add_argument('-p', '--path', metavar="<path/to/>", type = readable_dir, required=True, dest = 'pathParam', help = 'Path to ')
	"""

	if not os.path.isdir(prospective_dir):
		raise argparse.ArgumentTypeError(" :{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK) == False :
		raise argparse.ArgumentTypeError(" :{0} is not a readable dir".format(prospective_dir))

def replace_all(repls, str):
	"""
	Function that take a dictionnary and text variable and return text variable with replace 'Key' from dictionnary with 'Value'.

	:param repls: a python dictionary
	:type repls: dict()
	:param str: a string where remplace some words
	:type str: str()
	:rtype: str()
	:return: - txt with replace 'Key' of dictionnary with 'Value' in the input txt

	Example:
		>>> text =  "i like apples, but pears scare me"
		>>> print(replace_all({"apple": "pear", "pear": "apple"}, text))
		i like pears, but apples scare me
	"""

	return re.sub('|'.join(re.escape(key) for key in repls.keys()),lambda k: repls[k.group(0)], str)


def loadInList(filename):
	"""
	Load file in list() and then remove \\n at end of line

	:param filename: a file
	:type filename: file
	:rtype: list()
	:return: - list of row's file without \\n
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> rows = loadInList(filename)
		>>> rows
		["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
	"""

	list = open(filename,"r").readlines()
	listgood=[line.rstrip() for line in list]
	return listgood

def loadInListCol(filename, col):
	"""
	Load file in list() and then remove \\n at end of line

	:param filename: a file
	:type filename: file
	:rtype: list()
	:return: - list of row's file without \\n
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> rows = loadInListCol(filename, 0)
		>>> rows
		["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
	"""

	list = open(filename,"r").readlines()
	listgood=[line.rstrip().split("\t")[col] for line in list]
	return listgood

def loadInListWithHeader(filename):
	"""
	Load file in list() and then remove \\n at end of line

	:param filename: a file
	:type filename: file
	:rtype: str(), list()
	:return: - header liste\n
			 - list of row's file without
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> header, rows = loadInListWithHeader(filename)
		>>> header
		"head1\thead2\thead3"
		>>> rows
		["i like pears, but apples scare me","i like apples, but pears scare me","End of file"]
	"""

	list = open(filename,"r").readlines()
	header = list[0]
	listgood=[line.rstrip() for line in list[1:]]
	return header, listgood

def loadInDict(filename):
	"""
	Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valueare other column.

	:param filename: a file
	:type filename: file
	:rtype: dict()
	:return: - dict of row's file without \\n with key is first column and value list of other column
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> dico = loadInDict(filename)
		>>> dico
		{
		"col1",["col2","col3"],
		"indiv1",["valeurcol2","valeurcol3"],
		"indiv2",["valeurcol2","valeurcol3"]
		}
	"""

	dicoOut={}
	with open(filename) as filein:
		for line in filein:
			tabLine = line.rstrip().split("\t")
			#print(tabLine[0], tabLine[1])
			if tabLine[0] not in dicoOut.keys():
				dicoOut[tabLine[0]] = []+tabLine[1:]
			#else:
				#dicoOut[tabLine[0]].append(tabLine[1])
	return dicoOut

def loadInDictCol(filename,columnkey, columnvalue):
	"""
	Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valu specify column.

	:param filename: a file
	:type filename: file
	:param columnkey: int of column
	:type columnkey: int
	:param columnvalue: int of column
	:type columnvalue: int
	:rtype: dict()
	:return: - dict of row's file without \\n with key is first column and value column number pass
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> dico = loadInDict(filename,columnkey=1,columnvalue=3 )
		>>> dico
		{
		"col1","col3",
		"indiv1","valeurcol3",
		"indiv2","valeurcol3"
		}
	"""

	dicoOut={}
	with open(filename) as filein:
		for line in filein:
			tabLine = line.rstrip().split("\t")
			#print(tabLine[0], tabLine[1])
			if tabLine[columnkey] not in dicoOut.keys():
				dicoOut[tabLine[columnkey]] = tabLine[columnvalue]
	return dicoOut

def loadInDictLine(filename):
	"""
	Load file in Dict() and then remove \\n at end of line, then add first column in key of dict and valueare other column.

	:param filename: a file
	:type filename: file
	:rtype: dict()
	:return: - dict of row's file without \\n with key is first column and value list of other column
	:warn: Use this function with small file !!! except more RAM are use and crash systeme.

	Example:
		>>> dico = loadInDictLine(filename)
		>>> dico
		{
		"col1",[line1],
		"indiv1",[line2],
		"indiv2",[line3]
		}
	"""

	dicoOut={}
	with open(filename) as filein:
		for line in filein:
			tabLine = line.rstrip().split("\t")
			#print(tabLine[0], tabLine[1])
			if tabLine[0] not in dicoOut.keys():
				dicoOut[tabLine[0]] = line
	return dicoOut


def extractListFromFasta(sequenceFile,FileList ):
	"""
	Function who use 2 files:\\n
			- fasta : file with all Sequences
			- list : file with name of sequence to extract from fasta

	then return dict() with only Sequence in listed file.

	:param sequenceFile: a file with all Sequences
	:type sequenceFile: file
	:param FileList: a file with name of sequence to extract from fasta
	:type FileList: file
	:rtype: dict() and int()
	:return: - dict with only Sequence in listed file and nbtotal seq in file
	:requires: this function require fasta2dict(filename) and loadInList(filename)

	Example:
		>>> dictSequences = extractListFromFasta(sequenceFile, FileList)
		>>> dictSequences
		{"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
	"""

	dicoOutput = {}
	# Ouverture des sequences fasta MGG et chargement dans dictionnaire
	dictSequences = fasta2dict(sequenceFile)
	# ouverture des identifiants a garder
	listKeep = loadInList(FileList)
	keep = 0
	noKeep = 0
	noKeepID=[]
	#for ID, record in dictSequences.items():
	for ID in sorted(dictSequences.keys(), key=sort_human):
		record = dictSequences[ID]
		if ID in listKeep:
			keep +=1
			dicoOutput[record] = record
		else:
			noKeepID.append(ID)
			noKeep += 1
	#print("seq keep:"+str(keep))
	#print("seq nokeep:"+str(noKeep))
	#print("ID nokeep:"+str(noKeepID))
	total = noKeep+keep
	return dicoOutput, total

def extractInverseListFromFasta(sequenceFile,FileList ):
	"""
	Function who use 2 files:\\n
			- fasta : file with all Sequences
			- list : file with name of sequence to extract from fasta

	then return dict() with only Sequence not in listed file.

	:param sequenceFile: a file with all Sequences
	:type sequenceFile: file
	:param FileList: a file with name of sequence to not extract from fasta
	:type FileList: file
	:rtype: dict() and int()
	:return: - dict with only Sequence in listed file and nbtotal seq in file
	:requires: this function require fasta2dict(filename) and loadInList(filename)

	Example:
		>>> dictSequences = extractListFromFasta(sequenceFile, FileList)
		>>> dictSequences
		{"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
	"""

	dicoOutput = {}
	# Ouverture des sequences fasta MGG et chargement dans dictionnaire
	dictSequences = fasta2dict(sequenceFile)
	# ouverture des identifiants a garder
	listnotKeep = loadInList(FileList)
	keep = 0
	noKeep = 0
	noKeepID=[]
	#for ID, record in dictSequences.items():
	for ID in sorted(dictSequences.keys(), key=sort_human):
		record = dictSequences[ID]
		if ID not in listnotKeep:
			keep +=1
			dicoOutput[record] = record
		else:
			noKeepID.append(ID)
			noKeep += 1
	total = noKeep+keep
	#print("seq keep:"+str(keep))
	#print("seq nokeep:"+str(noKeep))
	#print("ID nokeep:"+str(noKeepID))
	return dicoOutput, total

def lsDirToList(pathDirectory):
	"""
	Return a list of file and directory find in directory

	:param pathDirectory: a directory Path
	:type pathDirectory: Path
	:rtype: list()
	:return: list of filename in pathDirectory

	Example:
		>>> lsDirectory = lsDirToList(path/to/directory/)
		>>> print(lsDirectory)
		["./out/gemo10_4497_ortho_rename_add.fasta", "./out/gemo10_6825_ortho_rename_add.fasta", "./out/gemo10_3497_ortho_rename_add.fasta", "./out/rename/"]
	"""

	if pathDirectory[-1] != "/":
		pathDirectory += "/"
	if pathDirectory[-1] != "*":
		pathDirectory += "*"
	lsFiles=glob.glob(pathDirectory)
	return lsFiles

def lsFastaInDirToList(pathDirectory):
	"""
	Return a list of fasta file's find in directory ("fasta", "fa", "fas")

	:param pathDirectory: a directory Path
	:type pathDirectory: Path
	:rtype: list()
	:return: list of fasta filename in pathDirectory ( file with extention "fa", "fasta", "fas" )

	Example:
		>>> lsDirectory = lsFastaInDirToList(path/to/directory/)
		>>> print(lsDirectory)
		["./out/gemo10_4497_ortho_rename_add.fasta", "./out/gemo10_6825_ortho_rename_add.fasta", "./out/gemo10_3497_ortho_rename_add.fasta"]
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
			nameFichier = fichier.split("/")[-1].split(".")[0]
			extentionFichier = fichier.split("/")[-1].split(".")[-1]
		except:
			extentionFichier = "directory"
		if extentionFichier in ["fasta", "fa", "fas"]:
			lsFilesFasta.append(fichier)

	return sorted(lsFilesFasta)

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
			lsFilesFasta.append(fichier)

	return sorted(lsFilesFasta)


def concatFastasFiles(pathDirectory):
	"""
	Return a fasta dictionnary of concatenation fasta file's find in directory ("fasta", "fa", "fas")

	:param pathDirectory: a directory Path
	:type pathDirectory: Path
	:rtype: dict()
	:return: dict of concatenation fasta filename in pathDirectory ( file with extention "fa", "fasta", "fas" )

	Example:
		>>> dico_concate = concatFastasFiles(path/to/directory/)
		>>> print(dico_concate)
			{"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
	"""
	outputDicoSeqs = {}
	listFiles = []

	if pathDirectory[-1] != "/":
		pathDirectory += "/"
	if pathDirectory[-1] != "*":
		pathDirectory += "*"

	directoryFiles=glob.glob(pathDirectory)
	# Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
	for fichier in directoryFiles:
		try:
			nameFichier = fichier.split("/")[-1].split(".")[0]
			extentionFichier = fichier.split("/")[-1].split(".")[1]
		except:
			extentionFichier = "directory"
		if extentionFichier in ["fasta", "fa", "fas"]:
			listFiles.append(fichier)

	try:
		for fasta in listFiles:
			#Ouverture des sequences fasta et chargement dans dictionnaire
			file = "%s" %(fasta)
			#print(file)
			record_dict = fasta2dict(file)
			for ind in record_dict.keys():
				if ind not in outputDicoSeqs.keys():
					name=record_dict[ind].id
					#print(name)
					seq = record_dict[ind].seq
					outputDicoSeqs[name]=seq
				else:
					name=record_dict[ind].id
					#print(name)
					seq = record_dict[ind].seq
					outputDicoSeqs[name]+=seq
	except:
		print(fasta)
		print(name)
		print(outputDicoSeqs[name])
	return outputDicoSeqs



def convertFasta2Nexus(pathDirectoryIn, pathDirectoryOut):
	"""
	Return the number of fasta file's convert  find in directory ("fasta", "fa", "fas") where are converted

	:param pathDirectoryIn: a directory Path
	:type pathDirectoryIn: Path
	:param pathDirectoryOut: a directory Path
	:type pathDirectoryOut: Path
	:rtype: int()
	:return: pathDirectory with file's converted and nb file converted( file with extention "nex" )

	Example:
		>>> nbFileConvert = convertFasta2Nexus(path/to/directory/, path/to/directory/)
		>>> print(nbFileConvert)
			"4172"
	"""
	countFilesConvert = 0
	fastaFiles = lsFastaInDirToList(pathDirectoryIn)
	if pathDirectoryIn[-1] != "/":
		pathDirectoryIn += "/"
	if pathDirectoryOut[-1] != "/":
		pathDirectoryOut += "/"

	# general variables:
	minimal_record = "#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=%s; end;" % "dna"

	for fileFasta in fastaFiles:
		try:
			countFilesConvert += 1
			filefastaName = fileFasta.split("/")[-1].split(".")[0]
			alignment = AlignIO.read(fileFasta, format='fasta')
			lenAlignement = int(alignment.get_alignment_length())
			n = Nexus.Nexus(minimal_record)
			n.alphabet = alignment._alphabet
			for record in alignment:
				n.add_sequence(record.id.replace("-","_"), str(record.seq))
			n.write_nexus_data(pathDirectoryOut+filefastaName+".nex", interleave=False)
		except ValueError:
			print(fileFasta)
			print(lenAlignement)


	return countFilesConvert

		#input_handle = open(fileFasta, "rU")
		#filefastaName = fileFasta.split("/")[-1].split(".")[0]
		#output_handle = open(pathDirectoryOut+filefastaName+".nex", "w")

		#alignments = AlignIO.read(input_handle, "fasta", alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
		#AlignIO.write(alignments, output_handle, "nexus", interleave="FALSE")

		#outputfilename.close()
		#input_handle.close()


def printcolor(txt, color, noprint=1):
	"""	Return the printed color txt format

	:param txt: a string
	:type txt: string
	:param color: a color value
	:type color: string
	:type noprint: int 0=noprint 1=print (default)
	:rtype: string()
	:return: string with acci color for printed
	:warn: List of avail color: reset, hicolor, underline, inverse, fblack, fred, fgreen, fyellow, fblue, fmagenta, fcyan, fwhite, bblack, bred, bgreen, byellow, bblue, bmagenta, bcyan, bwhite

	Example:
		>>> printcolor("il fait beau aujourd'hui","bgreen")
			"\\033[36mil fait beau aujourd'hui"
		>>> txtcolor = printcolor("il fait beau aujourd'hui","bgreen", 0)

	"""
	dicoColor = {
	"reset":"\033[0m","hicolor":"\033[1m","underline":"\033[4m","inverse":"\033[7m",
	"fblack":"\033[30m","fred":"\033[31m","fgreen":"\033[32m","fyellow":"\033[1;33m",
	"fblue":"\033[34m","fmagenta":"\033[35m","fcyan":"\033[36m","fwhite":"\033[37m",
	"bblack":"\033[40m","bred":"\033[41m","bgreen":"\033[42m","byellow":"\033[43m",
	"bblue":"\033[44m","bmagenta":"\033[45m","bcyan":"\033[46m","bwhite":"\033[47m",
	}
	if color in dicoColor.keys():
		txtout = dicoColor[color]+txt
		if noprint == 0:
			return txtout
		else:
			print(txtout)
	else:
		txtout = "Error, color value non exist, please check other color\n\n"+txt

def relativeToAbsolutePath(relative):
	"""	Return the absolutPath

	:param relative: a string path
	:type relative: string
	:rtype: string()
	:return: absolutePath
	:warn: need subprocess::check_output

	Example:
	>>>print(relative)
	../test
	>>> pathDirectory = relativeToAbsolutePath(relative)
	>>>print(pathDirectory)
	/home/sebastien/test

	"""
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = subprocess.check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:						# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath



#################################################
# CLASS
#################################################

class parseGFF():

	def __init__(self, filename):
		#Initialized GeneInfo named tuple. Note: namedtuple is immutable
		self.filename = filename
		self.gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
		self.GFFRecord = namedtuple("GFFRecord", self.gffInfoFields)

	def parseGFFAttributes(self, attributeString):
		"""Parse the GFF3 attribute column and return a dict"""
		if attributeString == ".": return {}
		ret = {}
		for attribute in attributeString.split(";"):
			key, value = attribute.split("=")
			ret[urllib.unquote(key)] = urllib.unquote(value)
		return ret

	def parseGFF3(self):
		"""
		A minimalistic GFF3 format parser.
		Yields objects that contain info about a single GFF3 feature.

		Supports transparent gzip decompression.
		"""
		#Parse with transparent decompression
		openFunc = gzip.open if self.filename.endswith(".gz") else open
		with openFunc(self.filename) as infile:
			for line in infile:
				if line.startswith("#"): continue
				parts = line.strip().split("\t")
				#If this fails, the file format is not standard-compatible
				assert len(parts) == len(self.gffInfoFields)
				#Normalize data
				normalizedInfo = {
					"seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
					"source": None if parts[1] == "." else urllib.unquote(parts[1]),
					"type": None if parts[2] == "." else urllib.unquote(parts[2]),
					"start": None if parts[3] == "." else int(parts[3]),
					"end": None if parts[4] == "." else int(parts[4]),
					"score": None if parts[5] == "." else float(parts[5]),
					"strand": None if parts[6] == "." else urllib.unquote(parts[6]),
					"phase": None if parts[7] == "." else urllib.unquote(parts[7]),
					"attributes": self.parseGFFAttributes(parts[8])
				}
				#Alternatively, you can emit the dictionary here, if you need mutability:
				#	yield normalizedInfo
				yield self.GFFRecord(**normalizedInfo)


class printCol():
	"""
	printCol() CLASS
	===============

	Classe qui ajoute des méthodes à print pour afficher de la couleur

	Example:

	>>> printCol.red("j'affiche en rouge")
	j'affiche en rouge

	"""

	RED = '\033[91m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	LIGHT_PURPLE = '\033[94m'
	PURPLE = '\033[95m'
	END = '\033[0m'

	@classmethod
	def red(cls, s):
		print(cls.RED + s + cls.END)

	@classmethod
	def green(cls, s):
		print(cls.GREEN + s + cls.END)

	@classmethod
	def yellow(cls, s):
		print(cls.YELLOW + s + cls.END)

	@classmethod
	def lightPurple(cls, s):
		print(cls.LIGHT_PURPLE + s + cls.END)

	@classmethod
	def purple(cls, s):
		print(cls.PURPLE + s + cls.END)

class AutoVivification(dict):
	"""
	AutoVivification(dict) CLASS
	============================

	Implementation of perl's autovivification feature.

	Example:

	>>> a = AutoVivification()
	>>> a[1][2][3] = 4
	>>> a[1][3][3] = 5
	>>> a[1][2]['test'] = 6
	>>> print a
	>>> {1: {2: {'test': 6, 3: 4}, 3: {3: 5}}}

	"""
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value

#*********************************************** Classe directory *******************
class directory(str):
	"""
	directory() CLASS
	=================

	Class which derives from string.
	Checks that the string is and path to valid directory and not empty

	Returns object which able to return basename and list file (with an extention, else consider as folder)

	Example:

	>>> inDirectory=directory("/home/sravel/Documents")
	>>> print(inDirectory.pathDirectory())
	>>> /home/sravel/Documents/

	>>> print(inDirectory.listFiles())
	>>> ["File1.txt","File2.pl","file.toto"]
	"""
	def __init__(self, pathDirectory = None):
		"""
			Initialise variable
		"""
		self.listPath = []			# all in the path
		self.listDir = []			# only directory in path
		self.listFiles = []			# only files in path

		self.current_dir = os.path.dirname(os.path.abspath(__file__))

		#Change relative path to absolute path
		self.pathDirectory = relativeToAbsolutePath(pathDirectory)

		# appel les fonctions
		self.testDirExist()
		self.lsInDir()
		self.splitFilesDir()

	def __repr__(self):
		"""Fonction qui permet de formater le text de sortie lors du print du dictionnaire"""
		txtOut = """
pathDirectory=%s\n
listPath=%s\n
listDir=%s\n
listFiles=%s\n
""" % (self.pathDirectory, str(self.listPath), str(self.listDir), str(self.listFiles))
		return txtOut

	def testDirExist(self):
		"""Test l'existance du répertoire"""
		if os.path.isdir(self.pathDirectory) != True :
			print("ERROR MODULES_SEB::Class-directory : path '%s' is not valide path" % self.pathDirectory )
			exit()

	def lsInDir(self):
		"""List all in directory"""
		self.testDirExist()
		if self.pathDirectory[-1] != "/":
			self.pathDirectory += "/"
		if self.pathDirectory[-1] != "*":
			pathDirectoryList = self.pathDirectory+"*"
		self.listPath=glob.glob(pathDirectoryList)

	def lsExtInDirToList(self, ext):
		"""List specific extention file in directory"""
		lsFilesFasta =[]
		for fichier in self.listPath:
			try:
				if "." in fichier.split("/")[-1]:
					nameFichier = fichier.split("/")[-1].split(".")[0]
					extentionFichier = fichier.split("/")[-1].split(".")[1]
				else:
					nameFichier = fichier.split("/")[-1]
					extentionFichier = ""
			except:
				extentionFichier = "directory"
			if extentionFichier == ext:
				lsFilesFasta.append(fichier)
		return sorted(lsFilesFasta, key=sort_human)

	def splitFilesDir(self):
		"""list files and list directory"""
		self.lsInDir()
		self.listDir = []		# only directory in path
		self.listFiles = []		# only files in path
		# Ouverture des fichiers du repertoire
		for fichier in self.listPath:
			try:
				if os.path.isdir(fichier) == True :		# C'est un dossier
					self.listDir.append(str(fichier))
				elif os.path.exists(fichier) == True :	# C'est un fichier
					self.listFiles.append(str(fichier))
			except:
				print("ERROR MODULES_SEB::Class-directory : path '%s' is not valide path contain other type (not files or directory)" % self.pathDirectory)
				exit()



#*********************************************** Classe StrList *******************
class StrList(list):
	"""
	StrList() CLASS
	===============

	Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien dans la liste "expectedValue"

	Example:

	>>> x=StrList(expectedValue=['true'])
	>>> x.append("true")
	>>> x
	['true']
	>>> x.append("false")
	Traceback (most recent call last):
		...
	ValueError: Input value is not in the following list: ['true']
	"""
	def __init__(self, expectedValue = None, defaultvalue = None):
		"""
		StrList(expectedValue=[])

		:rtype: list()
		:return: list check if not value of past list
		"""
		self.expectedValue = expectedValue
		self.defaultValue = defaultvalue
		super(StrList, self).__init__()

	def append(self, value):
		"""Redefine append() function to do testItem before append in list"""
		self.testItem(value)
		super(StrList, self).append(value)

	def insert(self, at, value):
		"""Redefine insert() function to do testItem before insert in list"""
		self.testItem(value)
		return list.insert(self, at, value)

	def testItem(self, value):
		"""Méthode qui test la valeur a ajouter est dans la list des expected"""
		if (value != None) and (value != "") and (self.expectedValue != None):
			if value in self.expectedValue:
				pass
			else:
				raise ValueError("Input value is not in the following list: %s" % self.expectedValue)

# *********************************************** Classe NumberList *******************
class NumberList(list):
	"""Classe qui dérive de list et qui vérifie que le nombre ajouté dans l'objet (list) est bien du type défini dans "dataType"
	Returns list

	Example:

	>>> x=NumberList(datatype = int, minValue = 2, maxValue = 100)
	>>> x.append("32")
	>>> x
	['32']
	>>> x.append("0.5")
	Traceback (most recent call last):
		...
	ValueError: "0.5" is not type "int"
	>>> x.append("1")
	Traceback (most recent call last):
		...
	ValueError: "1" is less than "2" ( 2 < value < 100 )
	>>> x.append("1000")
	Traceback (most recent call last):
		...
	ValueError: "1000" is greater than "100" ( 2 < value < 100 )
	"""
	def __init__(self, datatype = None, minValue = None, maxValue = None, defaultvalue = None):
		"""
		:rtype: list()
		:return: list check if not value of past list
		"""
		self.dataType = datatype
		self.min = minValue
		self.max = maxValue
		self.defaultValue = defaultvalue
		super(NumberList, self).__init__()

	def append(self, value):
		"""L.append(object) -- append object to end after do testItem"""
		self.testItem(value)
		super(NumberList, self).append(value)

	def insert(self, at, value):
		self.testItem(value)
		return list.insert(self, at, value)

	def testItem(self, value):
		"""Méthode qui test la valeur a ajouter"""

		if (value != None) and (value != ""):
			try:
				self.dataType(value)
			except :
				raise ValueError("\"%s\" is not type \"%s\"" % (value , self.dataType.__name__))
			try:
				if self.min < self.dataType(value):
					pass
				else:
					raise ValueError("\"%s\" is less than \"%s\" ( %s < value < %s )" % (value , self.min, self.min, self.max))
				if self.dataType(value) < self.max:
					pass
				else:
					raise ValueError("\"%s\" is greater than \"%s\" ( %s < value < %s )" % (value , self.max, self.min, self.max))
			except ValueError as infos:
				raise ValueError("%s" % infos)

# ************************************************************* Classe StrSpaceList *******************
class StrSpaceList(list):
	""""Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien du type défini dans "dataType"
	pour une chaine de plusieurs valeurs. test également si la valeur et comprise entre 2 bornes; et le nobre de valeur attendu.
	Returns list

	Example:

	>>> x=StrSpaceList(datatype = float, minValue = 0, maxValue = 1, nbInListMin = 4, nbInListMax = 4, sumOfValues = 1)
	>>> x.append("0.1 0.2 0.3 0.4")
	>>> x
	['0.1 0.2 0.3 0.4']
	>>> x.append("0.5")
	Traceback (most recent call last):
		...
	ValueError: You must enter values more than 1 ​​separated by a space ( 4 < nbValue < 4 )
	>>> x.append("2 0.2 0.3 0.4")
	Traceback (most recent call last):
		...
	ValueError: "2" is greater than "1" ( 0 < value < 1 )
	>>> x.append("-1 0.2 0.3 0.4")
	Traceback (most recent call last):
		...
	ValueError: "-1" is less than "0" ( 0 < value < 1 )
	"""
	def __init__(self, datatype = None, minValue = None, maxValue = None, nbInListMin = None, nbInListMax = None, sumOfValues = None, defaultvalue = None):
		"""
		:rtype: list()
		:return: list check if not value of past list
		"""
		self.dataType = datatype
		self.min = minValue
		self.max = maxValue
		self.nbInListMin = nbInListMin
		self.nbInListMax = nbInListMax
		self.sumOfValues = sumOfValues
		self.defaultValue = defaultvalue

	def append(self, value):
		"""L.append(object) -- append object to end after do testItem"""
		self.testItem(value)
		super(StrSpaceList, self).append(value)

	def insert(self, at, value):
		self.testItem(value)
		return list.insert(self, at, value)

	def testItem(self,value):
		"""Méthode qui test la valeur a ajouter"""
		if (value != None) and (value != ""):
			l = value.split(' ')

			for num in l:
				try:
					self.dataType(num)
				except :
					raise ValueError("\"%s\" is not type \"%s\"" % (num , self.dataType.__name__))
				try:
					if self.min < self.dataType(num):
						pass
					else:
						raise ValueError("\"%s\" is less than \"%s\" ( %s < value < %s )" % (num , self.min, self.min, self.max))
					if self.dataType(num) < self.max:
						pass
					else:
						raise ValueError("\"%s\" is greater than \"%s\" ( %s < value < %s )" % (num , self.max, self.min, self.max))
				except ValueError as infos:
					raise ValueError("%s" % infos)
			try:
				if len(l) < self.nbInListMin:
					raise ValueError("You must enter values more than %s ​​separated by a space ( %s < nbValue < %s )" %(len(l), self.nbInListMin, self.nbInListMax))
				if len(l) > self.nbInListMax:
					raise ValueError("You must enter value less than %s ​​separated by a space ( %s < nbValue < %s )" %(len(l), self.nbInListMinn, self.nbInListMax))
			except ValueError as infos:
				raise ValueError("%s" % infos)

			if self.sumOfValues != None:
				s = 0
				for num in l:
					s+= self.dataType(num)
				if s != self.sumOfValues:
					raise ValueError("The values must sum up to %s." % self.sumOfValues)

# ************************************************************* Classe StrSpace *******************
class LetterList(list):
	"""Class which derives from list. Checks that the string added to the object (list) contains only the letters from the list "Letters"
	Returns list

	Example:

	>>> x=LetterList(letters=['A', 'C', 'T', 'G'])
	>>> x.append("cagtcgatgcatgctagctagtcagtcat")
	>>> x
	['cagtcgatgcatgctagctagtcagtcat']
	>>> x.append("CGATCGATCGATCGT")
	>>> x
	['cagtcgatgcatgctagctagtcagtcat', 'CGATCGATCGATCGT']
	>>> x.append("CGAGFBDCGATCGATCGT")
	Traceback (most recent call last):
		...
	ValueError: "FBD" is not in "['A', 'C', 'T', 'G']"
	"""
	def __init__(self, letters = None, defaultvalue = None):
		"""
		:rtype: list()
		:return: list check if not value of past list
		"""
		self.Letter = letters
		self.defaultValue = defaultvalue

	def append(self, value):
		"""L.append(object) -- append object to end after do testItem"""
		self.testItem(value)
		super(LetterList, self).append(value)

	def insert(self, at, value):
		self.testItem(value)
		return list.insert(self, at, value)

	def testItem(self,value):
		"""Méthode qui test la valeur a ajouter"""
		if (value != None) and (value != ""):
			value = value.upper()
			ln = ""
			for l in value:
				if not l in self.Letter:
					ln+=l
			if ln != "":
				raise ValueError("\"%s\" is not in \"%s\"" % (ln , self.Letter))
