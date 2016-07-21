#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package getdistrgenotypeTAB.py
# @author Sebastien Ravel

"""
	The getdistrgenotypeTAB script
	==============================
	:author: Sebastien Ravel
	:contact: sebastien.ravel@cirad.fr
	:date: 08/07/2016
	:version: 0.1

	Script description
	------------------

	This Programme take proteineOrtho output and take file with Orthologue 1/1

	Example
	-------

	>>> getdistrgenotypeTAB.py -i SNP_table.tab

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display getdistrgenotypeTAB.py version number and exit

	Input mandatory infos for running:
		- \-i <filename>, --input <filename>
						tab file with SNP

"""

##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')
from MODULES_SEB import relativeToAbsolutePath,existant_file

## Python modules
import argparse
from time import localtime, strftime

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='18/05/2016'
debug="False"
#debug="True"


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())
# 	start=time.clock()
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='getdistrgenotypeTAB.py', description='''This Programme count allele ref/snp for position from Tab file''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display getdistrgenotypeTAB version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-i', '--input', metavar="<filename>",type=existant_file, required=True, dest = 'paramfile', help = 'tab file with SNP')

	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("#################################################################")
	print("#            Welcome in getdistrgenotypeTAB (Version " + version + ")              #")
	print("#################################################################")
	print('Start time: ', start_time,'\n')

	# Récupère le fichier de conf passer en argument
	basename = args.paramfile.split(".")[0]
	tabFile = relativeToAbsolutePath(args.paramfile)
	workingDir = "/".join(tabFile.split("/")[:-1])
	print("\t - Working directory is: %s" % workingDir)


	nblignetotal = float(os.popen("wc -l "+tabFile).read().rstrip().split(" ")[0])


	print("Parse tab file with %i lines" % nblignetotal)

	# ajoute à la variable current_dir le chemin ou est executer le script
	current_dir = os.path.dirname(os.path.abspath(__file__))

	# Utilisation du VCF
	# lecture et ecriture du header dans le fichier de output

	ctr=0
	dicoCount={}
	snponly={}
	dicoCountU={}
	dicoCountF={}
	dicoCountUF={}
	dicoCountUFR={}
	dicoCountR={}
	dicoCountSNP={}
	dicoCountUSOUCHE={}
	dicoCountRSOUCHE={}
	dicoCountSNPSOUCHE={}
	dicoCountFSOUCHE={}
	nbligneremove=0

	with open(tabFile, "r") as tabFileOpen:

		header = tabFileOpen.readline()
		lTAB = header.rstrip().split("\t")
		tabnamesouche=lTAB[3:]
		NBsouche = len(tabnamesouche)
		for line in tabFileOpen:
			ctr += 1
			if ((ctr % 10000 == 0) and (ctr != 0)) or (float(ctr) == nblignetotal):
				percent = (float(ctr)/float(nblignetotal))*100
				sys.stdout.write("\rProcessed up to %0.2f %%..." % percent)
				sys.stdout.flush()
			try:
	#CHROM	POS	reference	BD0024	BR0026	CD0073	CD0203	CH0052	CH0063	CH0092	CH0333	CH0549	CL0026	CL3-6-7	HN001	IN0072	IN0082	IN0094	JP0010	MC0016	MD0929	ML0025	PH0103	PH0118	PR0009	SP0005	US0031	US0098	TH0012-Gandalf	BR32-gemo	CD156-gemo	PH14-gemo	FR13-gemo	TH16-gemo	GY11-gemo	CH532-Terauchi	CH1019-Terauchi	CH689-Terauchi	CH999-Terauchi	LA5-Terauchi	CH860-Terauchi	CH328-Terauchi	LA21-Terauchi	CH680-Terauchi	CH595-Terauchi	NP52-Terauchi	TH17-Terauchi	CH701-Terauchi	TH14-Terauchi	NP37-Terauchi	NP41-Terauchi	NP61-Terauchi	CH1016-Terauchi
				ref, unmap, filtered, snp, indice = 0 ,0 ,0 ,0 ,0
				tabindice=[]
				lTAB = line.rstrip().split("\t")
				tabsouche = lTAB[3:]
				for genotype in tabsouche:
					souche=tabnamesouche[indice]
					if genotype == "U":
						unmap+=1
						tabindice.append(0)
						dicoCountUSOUCHE[souche] = (dicoCountUSOUCHE.get(souche, 1))+1

					elif genotype == "F":
						filtered+=1
						dicoCountFSOUCHE[souche] = (dicoCountFSOUCHE.get(souche, 1))+1
						tabindice.append(0)

					elif genotype == "R":
						ref+=1
						tabindice.append(0)
						dicoCountRSOUCHE[souche] = (dicoCountRSOUCHE.get(souche, 1))+1

					elif genotype == "A" or genotype == "C" or genotype == "T" or genotype == "G":
						snp+=1
						tabindice.append(1)
					else:
						print("ERROR : file contain '%s' so can be use by the script !!" % genotype)
						exit()
					indice+=1
				#print(lTAB[1])
				if snp == 1:
					idsoucheonesnp = tabindice.index(1)
					soucheonesnp=tabnamesouche[idsoucheonesnp]

					snponly[soucheonesnp] = (snponly.get(soucheonesnp, 1))+1
					dicoCountSNPSOUCHE[soucheonesnp] = (dicoCountSNPSOUCHE.get(soucheonesnp, 1))+1

				else:
					i=0
					for idsouchesnp in tabindice:
						if idsouchesnp == 1:
							souchesnp=tabnamesouche[i]
							dicoCountSNPSOUCHE[souchesnp] = (dicoCountSNPSOUCHE.get(souchesnp, 1))+1
						i+=1

				if snp == 0:
					nbligneremove+=1
					print(lTAB[0],lTAB[1])
				nballele=snp+ref+unmap+filtered
				if nballele != int(NBsouche):
					print("Error : more or less SNP than souches")
					print(nballele)
					exit()
				UF=unmap+filtered
				UFR=UF+ref
				# Ajouter dico SNP
				if snp in dicoCountSNP.keys():
					dicoCountSNP[snp] = (dicoCountSNP.get(snp, 1))+1

				# Ajouter dico Unmap
				if unmap in dicoCountU.keys():
					dicoCountU[unmap] = (dicoCountU.get(unmap, 1))+1

				# Ajouter dico Filtered
				if filtered in dicoCountF.keys():
					dicoCountF[filtered] = (dicoCountF.get(filtered, 1))+1

				# Ajouter dico Unmap-Filter
				if UF in dicoCountUF.keys():
					dicoCountUF[UF] = (dicoCountUF.get(UF, 1))+1

				# Ajouter dico Unmap-Filter-Ref
				if UFR in dicoCountUFR.keys():
					dicoCountUFR[UFR] = (dicoCountUFR.get(UFR, 1))+1

				# Ajouter dico Ref
				if ref in dicoCountR.keys():
					dicoCountR[ref] = (dicoCountR.get(ref, 1))+1

			except ValueError as infos:
				print(line)

	print("\n\nBuild tableau 1")
	outfile = open(workingDir+"/"+basename+"_Stats.txt", "w")
	outfile.write("#######################\n## Stats Tableau 1\n#######################\n\n")

	dicoListPositionInfo={}
	headListPositionInfo=[]
	for i in range(0,int(NBsouche)+1):
		dicoListPositionInfo[i]=[]

	for key, dico in {"NbPosCountU":dicoCountU,"NbPosCountF":dicoCountF,"NbPosCountUF":dicoCountUF,"NbPosCountUFR":dicoCountUFR,"NbPosCountR":dicoCountR,"NbPosCountSNP":dicoCountSNP}.items():
		total=0
		#print("\n"+key)
		headListPositionInfo.append(key)
		for key in range(0,int(NBsouche)+1):
			if key in dico.keys():
				dicoListPositionInfo[key].append(str(dico[key]))
			else:
				dicoListPositionInfo[key].append("NA")
	#print(dicoListPositionInfo)

	head="NbSouche\t"
	head+="\t".join(headListPositionInfo)
	outfile.write(head+"\n")
	for key in sorted(dicoListPositionInfo.keys()):
		value="\t".join(dicoListPositionInfo[key])
		#print(str(key)+"\t"+value)
		outfile.write(str(key)+"\t"+value+"\n")

	#exit()
	##########################################################
	### Tableau 2
	##########################################################
	print("Build tableau 2")
	dicoListPositionInfoSouche={}
	headListPositionInfoSouche=[]
	for souche in tabnamesouche:
		dicoListPositionInfoSouche[souche]=[]


	outfile.write("\n\n#######################\n## Stats Tableau 2\n#######################\n\n")
	for key, dico in {"NbPosCountSNPSOUCHE":dicoCountSNPSOUCHE,"NbPosCountRSOUCHE":dicoCountRSOUCHE,"NbPosCountUSOUCHE":dicoCountUSOUCHE,"NbPosCountSNPonlySOUCHE":snponly,"NbPosCountFSOUCHE":dicoCountFSOUCHE}.items():
		total=0
		#print("\n"+key)
		headListPositionInfoSouche.append(key)

		for souche in tabnamesouche:
			if souche in dico.keys():
				dicoListPositionInfoSouche[souche].append(str(dico[souche]))
			else:
				dicoListPositionInfoSouche[souche].append("NA")



		#for key in sorted(dico.keys()):
			#print(key,dico[key])
			#if key in dicoListPositionInfoSouche.keys():
				#dicoListPositionInfoSouche[key].append(str(dico[key]))
			#else:
				#dicoListPositionInfoSouche[key] = []
				#dicoListPositionInfoSouche[key].append(str(dico[key]))

	#print(dicoListPositionInfoSouche)

	head="Souche\t"
	head+="\t".join(headListPositionInfoSouche)
	outfile.write(head+"\n")
	for key in sorted(dicoListPositionInfoSouche.keys()):
		#print(dicoListPositionInfoSouche[key])
		value="\t".join(dicoListPositionInfoSouche[key])
		#print(str(key)+"\t"+value)
		outfile.write(str(key)+"\t"+value+"\n")

	#print("\ndicoCount-SNP-only-SOUCHE")
	#for key in sorted(snponly.keys()):
		#print(key,snponly[key])

	outfile.write("\n\nnbligne total = "+str(nblignetotal))
	print("total nbligne sans SNP = "+str(nbligneremove))

	outfile.close()

	print("\n\nExecution summary:")
	print("  - Files outputting: ")
	print("  \t- "+tabFile.split(".")[0]+"_Stats.txt")

	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
