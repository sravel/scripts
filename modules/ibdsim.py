#!/usr/bin/python
# -*- coding: utf-8 -*-
## @package python_interface.ibdsim
# @author Sebastien Ravel

# variable à changer
VERSION='2.0'
VERSION_DATE='28/02/2012'
ACTIVE_GRAPH = True             # valeur True pour activer, False désactiver
namePDF = 'IBDSimV2.0_UserGuide.pdf'



from utils.cbgpUtils import Documentator, log
import copypython27 as copy
import dataPath
import output

import sys
import os
from subprocess import Popen

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import  *
from PyQt4.QtGui import *
from PyQt4 import uic

if ACTIVE_GRAPH == True:
    from PyQt4.Qwt5 import QwtPlot, QwtPlotCurve, Qwt


# les .app démarrent python avec '/' comme cwd
if "darwin" in sys.platform : # and ".app/" in sys.argv[0]:
    # pour aller a l'interieur du .app
    mycwd = sys.argv[0].split(".app/")[0] + ".app/Contents/Resources/"
    os.chdir(mycwd)




formIbdsim,baseIbdsim = uic.loadUiType("./uis/IBDSim.ui")
# *********************************************** Classe DictWithRollback *******************
class DictWithRollback(dict):
    """Class dérivant de "dict"
    Fonctionne sur le principe de la base de donnée: les données sont ajouter dans le dictionnaire mais ne sont valider
    que si aucune erreur n'est lever. En cas d'erreur, effectue un _rollback grâce à une copie antérieur à l'ajout"""
    def __init__(self,  keyvalues):
        self._keyValues = keyvalues
        super(DictWithRollback, self).__init__()
        self.__mode = "standard"

    def enterInSessionMode(self):
        """Fonction qui permet de faire une copie du Self afin de sauvegarder les données en cas d'erreur"""
        if self.__mode == "standard":
            self.__mode = "saving"
            self.backupDict = {}
            for key, value in self.items() :
                self.backupDict[key] = copy.deepcopy(value)
        else:
            raise ValueError, "It's not possible to make enterInSessionMode in saving mode"

    def _commit(self):
        """Fonction de suppression de la copie de DictWithRollback si auccune erreur"""
        if self.__mode == "saving":
            del self.backupDict
            self.__mode = "standard"
        elif self.__mode == "standard":
            raise ValueError, "It's not possible to make a commit in standard mode"

    def _rollback(self):
        """Fonction qui permet de faire un retour arrière grâce à la sauvegarde en cas d'erreurs"""
        if self.__mode == "saving":
            for key , value in self.backupDict.items():
                self[key]= value
            self.__mode = "standard"
        elif self.__mode == "standard":
            raise ValueError, "It's not possible to make a _rollback in standard mode"

    def __setitem__(self,x,y, obj = None):
#        """%s""" % dict.__doc__
        # dict.__item__.__doc__

        return super(DictWithRollback, self).__setitem__(x, y )



# *********************************************** Classe StrList *******************
class StrList(list):
    """Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien dans la liste "expectedValue"
    Returns list

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
        self.expectedValue = expectedValue
        self.defaultValue = defaultvalue
        super(StrList, self).__init__()

    def append(self, value):
        """L.append(object) -- append object to end after do testItem"""
        self.testItem(value)
        super(StrList, self).append(value)

    def insert(self, at, value):
        self.testItem(value)
        return list.insert(self, at, value)

    def testItem(self, value):
        """Méthode qui test la valeur a ajouter"""
        if (value != None) and (value != "") and (self.expectedValue != None):
            if value in self.expectedValue:
                pass
            else:
                raise ValueError, "Input value is not in the following list: %s" % self.expectedValue

# *********************************************** Classe NumberList *******************
class NumberList(list):
    """"Classe qui dérive de list et qui vérifie que le nombre ajouté dans l'objet (list) est bien du type défini dans "dataType"
    Returns list

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
                raise ValueError, "\"%s\" is not type \"%s\"" % (value , self.dataType.__name__)
            try:
                if self.min < self.dataType(value):
                    pass
                else:
                    raise ValueError, "\"%s\" is less than \"%s\" ( %s < value < %s )" % (value , self.min, self.min, self.max)
                if self.dataType(value) < self.max:
                    pass
                else:
                    raise ValueError, "\"%s\" is greater than \"%s\" ( %s < value < %s )" % (value , self.max, self.min, self.max)
            except ValueError as infos:
                raise ValueError, "%s" % infos

# ************************************************************* Classe StrSpaceList *******************
class StrSpaceList(list):
    """"Classe qui dérive de list et qui vérifie que la variable ajouter dans l'objet (list) est bien du type défini dans "dataType"
    pour une chaine de plusieurs valeurs. test également si la valeur et comprise entre 2 bornes; et le nobre de valeur attendu.
    Returns list

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
                    raise ValueError, "\"%s\" is not type \"%s\"" % (num , self.dataType.__name__)
                try:
                    if self.min < self.dataType(num):
                        pass
                    else:
                        raise ValueError, "\"%s\" is less than \"%s\" ( %s < value < %s )" % (num , self.min, self.min, self.max)
                    if self.dataType(num) < self.max:
                        pass
                    else:
                        raise ValueError, "\"%s\" is greater than \"%s\" ( %s < value < %s )" % (num , self.max, self.min, self.max)
                except ValueError as infos:
                    raise ValueError, "%s" % infos
            try:
                if len(l) < self.nbInListMin:
                    raise ValueError, "You must enter values more than %s ​​separated by a space ( %s < nbValue < %s )" %(len(l), self.nbInListMin, self.nbInListMax)
                if len(l) > self.nbInListMax:
                    raise ValueError, "You must enter value less than %s ​​separated by a space ( %s < nbValue < %s )" %(len(l), self.nbInListMinn, self.nbInListMax)
            except ValueError as infos:
                raise ValueError, "%s" % infos

            if self.sumOfValues != None:
                s = 0
                for num in l:
                    s+= self.dataType(num)
                if s != self.sumOfValues:
                    raise ValueError, "The values must sum up to %s." % self.sumOfValues

# ************************************************************* Classe StrSpace *******************
class LetterList(list):
    """Class which derives from list. Checks that the string added to the object (list) contains only the letters from the list "Letters"
    Returns list

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
    def __init__(self, letters = None, defaultvalue = None):#
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
                raise ValueError, "\"%s\" is not in \"%s\"" % (ln , self.Letter)

# ************************************************************* Classe DictIBDSimConfig Stockage de l'information récuperer dans l'interface/controle des données *******************
class DictIBDSimConfig(DictWithRollback):
    """Classe de configuration des paramètres d'IBDSim
    Contient les listes de dictionnaires des paramètres d'IBDSim, Chaque liste correspond à un onglet dans l'interface.
    Attention les nom de la liste sons aussi présent dans la class IBDSim pour la gestion de l'affichage"""
    def __init__(self):
        """Initialise le dictionnaire DictWithRollback avec les paramètres des listes"""
        super(DictIBDSimConfig, self).__init__(self)

        self.allModels = ["KAM","IAM", "SMM", "TPM", "GSM", "SNP","ISM", "JC69", "K80", "F81", "HKY85", "TN93"]
        
        self.listParametersModel = [
                # Markers Parameters:
                      {'Locus_Number': NumberList(datatype= int, minValue = 0, maxValue = 999999999)}, 
                      {'Mutation_Rate': NumberList(datatype= float, minValue = 0, maxValue = 1)},
                      {'Mutation_Model': StrList( expectedValue = self.allModels )},
                      {'Variable_Mutation_Rate': StrList(expectedValue = ['False', 'True'])},
                      {'Polymorphic_Loci_Only': StrList(expectedValue = ['False', 'True'])},
                      {'Minor_Allele_Frequency':NumberList(datatype= float, minValue = -0.0001, maxValue = 0.5001)},
                      {'Min_Allele_Number': NumberList(datatype= int, minValue = 0, maxValue = 999999999)},
                      {'Allelic_Lower_Bound': NumberList(datatype= int, minValue = 0, maxValue = 999999999)},
                      {'Allelic_Upper_Bound': NumberList(datatype= int, minValue = 0, maxValue = 999999999)},
                      {'Allelic_State_MRCA': NumberList(datatype= int, minValue = -0.5, maxValue = 999999999)},
                      {'Repeated_Motif_Size': NumberList(datatype= int, minValue = 0, maxValue = 999999999)},
                      {'SMM_Probability_In_TPM': NumberList(datatype= float, minValue = -0.00001, maxValue = 1.00001)},
                      {'Geometric_Variance_In_TPM': NumberList(datatype= float, minValue = 0, maxValue = 999999999)},
                      {'Geometric_Variance_In_GSM': NumberList(datatype= float, minValue = -0.00001, maxValue = 999999999)},
                      {'Ploidy': StrList(expectedValue = ['Diploid','Haploid'])},
                # Sequence Specific settings
                      {'MRCA_Sequence': LetterList(letters = ['A','C','T','G'])},
                      {'Sequence_Size': NumberList(datatype= int, minValue = 0, maxValue = 999999999)},
                      {'Transition_Transversion_ratio': NumberList(datatype= float, minValue = -0.00001, maxValue = 1.00000001)},
                      {'Transition1_Transition2_ratio': NumberList(datatype= float, minValue = -0.00001, maxValue = 1.00000001)},
                      {'Equilibrium_Frequencies': StrSpaceList( datatype = float, minValue = 0, maxValue = 1, nbInListMin = 4, nbInListMax= 4, sumOfValues = 1)},
                    ]
        self.listParametersDemographic = [
                    # Demographic Change
                      {'New_Demographic_Phase_At': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                # Demographic Options
                    # Lattice
                      {'Lattice_Boundaries': StrList(expectedValue=['Absorbing', 'Circular','Reflecting'])},
                      {'Total_Range_Dispersal': StrList(expectedValue= ['True', 'False'])},
                      {'Lattice_SizeX': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 10)},
                      {'Lattice_SizeY': NumberList(datatype=int, minValue = 0, maxValue = 999999999 , defaultvalue = 1)},
                      {'Ind_Per_Pop': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Void_Nodes': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Specific_Density_Design': StrList(expectedValue= ['False', 'True'])},#'true', 'false'
                      {'Zone': StrList(expectedValue=['False', 'True'])},
                      {'Void_Nodes_Zone': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Ind_Per_Pop_Zone': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Min_Zone_CoordinateX': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Max_Zone_CoordinateX': NumberList(datatype=int, minValue = 0, maxValue = 999999999 , defaultvalue = 1)},
                      {'Min_Zone_CoordinateY': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Max_Zone_CoordinateY': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                    # Sample
                      {'Specific_Sample': StrList(expectedValue=['True', 'False'])},
                      {'Sample_SizeX': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Sample_SizeY': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Min_Sample_CoordinateX': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Min_Sample_CoordinateY': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Sample_Coordinates_X': StrSpaceList( datatype = int, minValue = 0, maxValue = 999999999, nbInListMin = 0, nbInListMax= 999999999, sumOfValues = None)},
                      {'Sample_Coordinates_Y': StrSpaceList( datatype = int, minValue = 0, maxValue = 999999999, nbInListMin = 0, nbInListMax= 999999999, sumOfValues = None)},
                      {'Ind_Per_Pop_Sampled': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                      {'Void_Sample_Node': NumberList(datatype=int, minValue = 0, maxValue = 999999999, defaultvalue = 1)},
                    # Dispersale
                      {'Dispersal_Distribution': StrList(expectedValue=['Stepping Stone', 'Geometric','Pareto','Sichel' , '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])},
                      {'Immigration_Control': StrList(expectedValue=['Simple_1D_Product', '1D_Product_Without_m0'])},
                      {'Total_Emigration_Rate': NumberList(datatype=float, minValue = -0.0001, maxValue = 1.00000001, defaultvalue = 0.5)},
                      {'Dist_max': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Pareto_Shape': NumberList(datatype=float, minValue = 0, maxValue = 999999999, defaultvalue = 0.5)},
                      {'Geometric_Shape': NumberList(datatype=float, minValue = -0.00001, maxValue = 1.00001, defaultvalue = 0.5)},
                      {'Sichel_Gamma': NumberList(datatype=float, minValue = -100, maxValue = 999999999, defaultvalue = -2.15)},
                      {'Sichel_Xi': NumberList(datatype=float, minValue = 0, maxValue = 999999999, defaultvalue = 100)},
                      {'Sichel_Omega': NumberList(datatype=float, minValue = -100, maxValue = 999999999, defaultvalue = -1)},
                    # Continuous Change in Density
                      {'Random_Translation': StrList(expectedValue=['True', 'False'])},
                      {'Continuous_Deme_Size_Variation': StrList(expectedValue = ['None', 'Linear', 'Exponential', 'Logistic'])},
                      {'Logistic_Growth_Rate': NumberList(datatype=float, minValue = -0.00001, maxValue = 999999999)},
                      {'Continuous_Lattice_Size_Variation': StrList(expectedValue = ['None', 'Linear', 'Exponential', 'Logistic'])},
                      {'Lattice_Logistic_Growth_Rate': NumberList(datatype=float, minValue = -0.00001, maxValue = 999999999)},
                      {'Barrier': StrList(expectedValue=['True', 'False'])},
                      {'X1_Barrier': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'X2_Barrier': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Y1_Barrier': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Y2_Barrier': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Barrier_Crossing_Rate': NumberList(datatype=float, minValue = -0.0001, maxValue = 1.0001)},
                      
                      ]
        self.listParametersComputation = [
                    # Various Computation options
                      {'DiagnosticTables': StrList(expectedValue=['Iterative_Identity_Probability', 'Hexp', 'Fis', 'Seq_stats', 'Prob_Id_Matrix', 'Effective_Dispersal', 'Iterative_Statistics', 'Allelic_Variance'])},
                      ]
        self.listParametersSimulAndOutput = [
                    # Simulation
                      {'Data_File_Name': StrList(expectedValue = None)},
                      {'Gene_Pop_File_Extension': StrList(expectedValue= None)},
                      {'Run_Number': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Random_Seeds': NumberList(datatype=int, minValue = 0, maxValue = 999999999)},
                      {'Pause': StrList(expectedValue= ['Default','Final', 'OnError'])},
                      # Output format
                      {'Genepop': StrList(expectedValue= ['False', 'True'])},
                      {'Migrate': StrList(expectedValue= ['False', 'True'])},
                      {'Migrate_Letter': StrList(expectedValue= ['False', 'True'])},
                      {'Nexus_File_Format': StrList(expectedValue= ['Haplotypes_only', 'Haplotypes_and_Individuals'])},
                                            ]
        self.__updateFromList(listparam = self.listParametersModel)
        self.__updateFromList(listparam = self.listParametersDemographic)
        self.__updateFromList(listparam = self.listParametersComputation)
        self.__updateFromList(listparam = self.listParametersSimulAndOutput)

    def __updateFromList(self, listparam = None):
        """mise a jour du dictionnaire a patir de la liste"""
        for dico in listparam:
            self.update(dico)

    def commitIbdsimWithNone(self, useList = None):
        """Fonction qui ajoute les valeurs None à la liste fournie dans "useList" pour les paramètres non remplis.
         appel du commit de DictWithRollback pour la suppression de la copie"""
        self._commit()
        lenMax = 0
        for dico in eval("self.%s" % useList):
            for key in dico.keys():
                setNumber = len(self[key])
                if setNumber > lenMax:
                    lenMax = setNumber
        for dico in eval("self.%s" % useList):
            for key in dico.keys():
                if len(self[key]) == lenMax:
                    i = lenMax
                    while i > 0:
                        if self[key][i-1] == '':
                            self[key][i-1] = None
                        i-=1
                else:
                    self[key].append(None)


    def deleteValuesInDict(self, deletedNumber = None, useList = None):
        """Fonction qui permet de supprimer une colonne dans la liste des objets du dictionnaire 
        Attend en paramètre la colonne et la liste des paramètres à supprimer"""
        for dico in eval("self.%s" % useList):
            for key in dico.keys():
                if len(self[key]) > 0:
                    del( self[key][deletedNumber])
    
    def resetDemographicParam(self):
        """Réinitialise les paramètres Démographiques"""
        for dico in self.listParametersDemographic:
            for key in dico.keys():
                nb = len(self[key])
                while 0 < nb:
                    del( self[key][nb-1])
                    nb -= 1

    def __repr__(self):
        """Fonction qui permet de formater le text de sortie lors du print du dictionnaire"""
        
        text = "%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%\n"
        
        for dico in self.listParametersSimulAndOutput:
            for key in dico.keys() :
                if len(self[key]) > 0:
                    if self[key][0] != None:
                        text += "%s = %s\n" % (key,  self[key][0])
        
        
        text += "\n\n%%%%%%%%%% MARKERS PARAMETERS %%%%%%%%%%%%\n"

        for dico in self.listParametersModel:
            for key in dico.keys() :
#                text += "%s = %s\n" % (key,  self[key])
                nbNone = self[key].count(None)
                nbValues = len(self[key])
                if nbNone != nbValues:
                    text += "%s = " % (key)
                    if key in ['Equilibrium_Frequencies','Transition1_Transition2_ratio','Transition_Transversion_ratio','Sequence_Size','MRCA_Sequence','Ploidy']:
                        while self[key][nbValues-1] == None:
                            nbValues = nbValues -1
                        text += self[key][nbValues-1]+"\n"
                    else:
                        for v in self[key]:
                            if (v == None) or (v == ""):
                                text += "Na, "
                            else:
                                text += v+", "
                        text = text[0:-2]+"\n"

        text += "\n%%%%%%%%%% VARIOUS COMPUTATION OPTIONS %%%%%%%%%%%%\n"
        
        for dico in self.listParametersComputation:
                for key in dico.keys() :
                    if len(self[key]) != 0:
                        text += "%s = " % key
                        for n in range(len(self[key])):
                            text += "%s, " % (self[key][n])
                        text = text[0:-1]
                    
        text += "\n\n%%%%%%%%%% DEMOGRAPHIQUE OPTIONS %%%%%%%%%%%%\n"
                        
        for dico in self.listParametersDemographic:
                for key in dico.keys() :
                    if key != "Specific_Sample":
    #                    text += "%s = %s\n" % (key,  self[key])
                        taille = len(self[key])-1
                        if taille != -1:
                            if (self[key][0] != None) and (self[key] != []):
                                text += "%s = %s\n" % (key,  self[key][0])

        if taille > 0:
            if key != "Specific_Sample":
                n = 1
                while n <= taille:
                    text += "\n\n%% %s DEMOGRAPHIQUE CHANGE\n" % str(n)
                    for dico in self.listParametersDemographic:
                        for key in dico.keys() :
                            if (self[key][n] != None) and (self[key] != []):
                                if self[key][n] != self[key][n-1]:
                                    text += "%s = %s\n" % (key,  self[key][n])
                    n = n + 1

        return text


class ComboBoxAbstract(QtGui.QComboBox):
    """Classe qui dérive de Combobox. Prend une liste lors de l'initialisation.
    La liste contient les paramètres à ajouter à la combobox"""
    def __init__(self, contents = None):
        super(ComboBoxAbstract, self).__init__()
        self.insertItems(1,contents)

class QTableWidgetItemAbstract(QtGui.QTableWidgetItem):
    """Classe qui dérive de QtGui.QTableWidgetItem."""
    def __init__(self, name = None, content = None):
        self.name = name
        
        super(QTableWidgetItemAbstract, self).__init__(content)

# ************************************* CLASSE IBDSim Gestion de l'affichage et de la récupération de donnée ****************************************************
## @class IBDSim
# @brief Classe principale, fenêtre principale
class IBDSim( formIbdsim, baseIbdsim ):
    """ Classe principale qui est aussi la fenêtre principale de la GUI
    L'initialisation permet de creer un objet "DictIBDSimConfig"
    Elle contient également les dictionnaires:
    - dicoForTypeInputKeys: Qui permet de remplir les cases du tableau des paramétres démographiques avec le type d'objet attendu (Combo ou CheckBox)
    - actualizeCellDisabled: contient les éléments a giser en fonction de l'élément en cour.
    et la liste:
    - listKeyDemographicParam créer à partir des clés des dictionnaires du DictIBDSimConfig
    """
    def __init__(self,app,parent=None):
        super(IBDSim,self).__init__(parent)
        self.app = app

        # initialisation des variables
        self.ibdsimConfig = DictIBDSimConfig()
        self.__nbModel = 0
        #self.__nbDemo = 0

        self.dicoForTypeInputKeys = {
        'Lattice_Boundaries':'Combo', 'Dispersal_Distribution': 'Combo', 'Immigration_Control': 'Combo',
         'Continuous_Deme_Size_Variation': 'Combo', 'Continuous_Lattice_Size_Variation': 'Combo',
        'Total_Range_Dispersal': 'CheckBox', 'Specific_Density_Design':'CheckBox', 'Zone': 'CheckBox',
         'Random_Translation': 'CheckBox', 'Specific_Sample': 'CheckBox', 'Barrier': 'CheckBox'
        }
        
        self.listKeyDemographicParam = []
        for dico in self.ibdsimConfig.listParametersDemographic:
            for key in dico.keys():
                self.listKeyDemographicParam += eval("['%s']" % key)
        
        self.actualizeCellDisabled = {
        'Stepping Stone':['Dist_max','Geometric_Shape', 'Pareto_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        'Geometric':['Dist_max','Pareto_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        'Pareto':['Dist_max', 'Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        'Sichel':['Dist_max','Geometric_Shape', 'Pareto_Shape'],
        '0':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '1':['Dist_max','Geometric_Shape', 'Pareto_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '2':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '3':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '4':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '5':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '6':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '7':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '8':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'],
        '9':['Geometric_Shape', 'Sichel_Gamma', 'Sichel_Xi', 'Sichel_Omega'], 
        'Zone':['Void_Nodes_Zone', 'Ind_Per_Pop_Zone', 'Min_Zone_CoordinateX', 'Max_Zone_CoordinateX', 'Min_Zone_CoordinateY', 'Max_Zone_CoordinateY'],
        'None':['Logistic_Growth_Rate'],
        'Linear':['Logistic_Growth_Rate'],
        'Exponential':['Logistic_Growth_Rate'],
        'Logistic':[],
        'Specific_SampleUnCheck':['Sample_Coordinates_X', 'Sample_Coordinates_Y'],
        'Specific_SampleCheck':['Sample_SizeX','Sample_SizeY','Min_Sample_CoordinateX','Min_Sample_CoordinateY'],
        'Barrier': ['X1_Barrier','X2_Barrier','Y1_Barrier','Y2_Barrier', 'Barrier_Crossing_Rate'],
        }
        self.actualizeCellDisabled2 = {
        'None':['Lattice_Logistic_Growth_Rate'],
        'Linear':['Lattice_Logistic_Growth_Rate'],
        'Exponential':['Lattice_Logistic_Growth_Rate'],
        'Logistic':[],
                                       }

        self.valueForCell = {
        'Stepping Stone':{'Total_Emigration_Rate': 0.5},
        'Geometric':{'Total_Emigration_Rate': 0.5, 'Geometric_Shape': 0.5},
        'Pareto':{'Total_Emigration_Rate': 0.5, 'Pareto_Shape': 0.5},
        'Sichel':{'Sichel_Gamma': -2.15, 'Sichel_Omega': -1, 'Sichel_Xi': 100, 'Total_Emigration_Rate': 0.5},
        '0':{'Dist_max': 15, 'Total_Emigration_Rate': 0.599985 ,'Pareto_Shape':2.51829},
        '1':{'Total_Emigration_Rate': 2/3},
        '2':{'Dist_max': 49, 'Total_Emigration_Rate': 0.599985,'Pareto_Shape':3.79809435},
        '3':{'Dist_max': 48, 'Total_Emigration_Rate': 0.6,'Pareto_Shape':1.246085},
        '4':{'Dist_max': 48, 'Total_Emigration_Rate': 0.824095,'Pareto_Shape':4.1078739681},
        '5':{'Dist_max': 48, 'Total_Emigration_Rate': 0.913378,'Pareto_Shape':4.43153111547},
        '6':{'Dist_max': 48, 'Total_Emigration_Rate': 0.719326,'Pareto_Shape':2.0334337244},
        '7':{'Dist_max': 49, 'Total_Emigration_Rate': 0.702504,'Pareto_Shape':2.313010658},
        '8':{'Dist_max': 48, 'Total_Emigration_Rate': 0.678842,'Pareto_Shape':4.1598694692},
        '9':{'Dist_max': 48, 'Total_Emigration_Rate': 0.700013,'Pareto_Shape':2.74376568},
            }
        
        
        self.createWidgets()
        self.createdDocForWhatIsThis()
        self.createDemograficTable()
        self.createComputationListe()
        self.createdConnectionSimulation()


    def createWidgets(self):
        """Mise en place du masquage des frames et des connections de l'interface"""
        self.ui = self
        self.ui.setupUi(self)
        
        # Destruction de l'onglet graph si desactiver
        
        if ACTIVE_GRAPH == False:
            self.ui.tabWidget.removeTab(5)

        pic = QPixmap("./uis/icons/IBDSIM.png")
#        pic.scaled(10,10)
        self.ui.imgLabel.setPixmap(pic)
        txt = str(self.ui.welcomeLabel.text())
        txt = txt.replace('vvv',VERSION).replace('ddd',VERSION_DATE)
        self.ui.welcomeLabel.setText(txt)
        

        self.setWindowTitle('IBDSim %s (%s) Graphic Interface' % (VERSION,VERSION_DATE))
        
        self.setWindowIcon(QIcon("./uis/icons/IBDSIM.png"))
                # about window
        self.aboutWindow = uic.loadUi("uis/about.ui")
        self.aboutWindow.parent = self
        self.aboutWindow.setWindowTitle('About IBDSim')
        ui = self.aboutWindow
        ui.logoLabel.setPixmap(QPixmap("./uis/icons/IBDSIM.png"))
        txt = str(self.aboutWindow.infoLabel.text())
        txt = txt.replace('vvv',VERSION).replace('ddd',VERSION_DATE)
        self.aboutWindow.infoLabel.setText(txt)
        QObject.connect(ui.okButton,SIGNAL("clicked()"),self.aboutWindow.close)

        # Initialisation des frames et bouttons
        self.ui.newModelFrame.hide()
        self.ui.paramFrame.hide()
        self.ui.displayErrorEdit.hide()
        
        for model in self.ibdsimConfig.allModels:
            eval("self.ui.%sframe.hide()" % model)
        self.ui.saveModelPushButton.setEnabled (False)
        self.ui.deleteModelPushButton.setEnabled (False)
        self.ui.editModelCombo.setEnabled (False)
        ### Various Computation Options
        self.ui.deleteSelectParamCompPushButton.setEnabled(False)
        self.ui.nexusFormatFrame.hide()
        ###Æ Graph
        self.ui.viewGraphPushButton.setEnabled(False)
        self.ui.hideGraphPushButton.setEnabled(False)
        self.ui.graphFrame.hide()

        # ajout des panneau attention pour la sauvegarde
        pic = QPixmap("./uis/icons/attention.gif")
#        pic.scaled(10,10)
        self.ui.attentionModelLabel.setPixmap(pic)
        self.ui.attentionDemographicLabel.setPixmap(pic)
        self.ui.attentionModelFrame.hide()
        self.ui.attentionDemographicFrame.hide()

        # Edition des connect:
        QObject.connect(self.ui.tabWidget,SIGNAL("currentChanged (int)"),self.changeTab)
                # Model TAB
        QObject.connect(self.ui.newModelPushButton,SIGNAL("clicked()"),self.newModel)
        QObject.connect(self.ui.saveModelPushButton,SIGNAL("clicked()"),self.saveModel)
        QObject.connect(self.ui.deleteModelPushButton,SIGNAL("clicked()"),self.deleteModel)
        QObject.connect(self.ui.Mutation_ModelCombo,SIGNAL("currentIndexChanged(int)"),self.actualizeModel)
                #Computational Options
        QObject.connect(self.ui.selectCompPushButton,SIGNAL("clicked()"),self.addCompItem)
        QObject.connect(self.ui.deleteSelectParamCompPushButton,SIGNAL("clicked()"),self.deleteCompItem)
            # Simulation and Output
        QObject.connect(self.ui.saveConfigFilePushButton,SIGNAL("clicked()"),self.saveFile)
        QObject.connect(self.ui.runIBDSimPushButton,SIGNAL("clicked()"),self.runIBDSim)
        
    # Menu toolbar
        file_menu = self.ui.menubar.addMenu("&File")
        self.file_menu = file_menu
        self.setStyleSheet("QMenu {\
                color: black;\
                background-color: white;\
                margin: 2px; /* some spacing around the menu */\
                }\
                \
                QMenu::item {\
                padding: 2px 25px 2px 20px;\
                border: 1px solid transparent; /* reserve space for selection border */\
                }\
                \
                QMenu::item:selected {\
                border-color: darkblue;\
                /*background: rgba(100, 100, 100, 150);*/\
                background-color: #FFD800;\
                }\
                \
                QMenu::item:disabled {\
                /*border-color: darkblue;*/\
                /*background: rgba(100, 100, 100, 150);*/\
                background-color: #ECECEC;\
                color: gray;\
                }\
                \
                QMenu::icon:checked { /* appearance of a 'checked' icon */\
                background: gray;\
                border: 1px inset gray;\
                position: absolute;\
                top: 1px;\
                right: 1px;\
                bottom: 1px;\
                left: 1px;\
                }\
                \
                QMenu::separator {\
                height: 2px;\
                background: lightblue;\
                margin-left: 10px;\
                margin-right: 5px;\
                }\
                \
                QMenu::indicator {\
                width: 13px;\
                height: 13px;\
                }")

        self.ui.file_menu.addAction(QIcon("./uis/icons/save.png"),"&Save",self.saveFile,QKeySequence(Qt.CTRL + Qt.Key_S))
        self.ui.file_menu.addAction(QIcon("./uis/icons/gtk-save-as.png"),"&Save As",self.saveAs,QKeySequence(Qt.CTRL + Qt.Key_E))
        self.ui.file_menu.addAction(QIcon("./uis/icons/IBDSIM.png"),"&Run IBDSim Program",self.runIBDSim,QKeySequence(Qt.CTRL + Qt.Key_R))
        self.ui.file_menu.addAction(QIcon("./uis/icons/window-close.png"),"&Exit",self.close,QKeySequence(Qt.CTRL + Qt.Key_Q))
        help_menu = self.ui.menubar.addMenu("&Help")
        self.help_menu = help_menu
        help_menu.addAction(QIcon("./uis/icons/dialog-question.png"),"&About IBDSim",self.aboutWindow.show,QKeySequence(Qt.CTRL + Qt.Key_A))
        help_menu.addAction(QIcon("./uis/icons/book_open.png"),"&Open User Guide pdf",self.openDoc,QKeySequence(Qt.CTRL + Qt.Key_H))

        # TOOLBAR
        saveButton = QPushButton(QIcon("./uis/icons/save.png"),"Save",self)
        self.saveButton = saveButton
        saveButton.setObjectName("saveButton")
        saveButton.setToolTip("Save current file in current path with the name IbdSettings.txt")
        saveButton.setMaximumSize(QSize(66, 22))
        #saveButton.setMinimumSize(QSize(16, 18))
        saveButton.setFlat(True)
        QObject.connect(saveButton,SIGNAL("clicked()"),self.saveFile)
        self.ui.toolBar.addWidget(saveButton)
        #self.ui.toolBar.addSeparator()

        saveAsButton = QPushButton(QIcon("./uis/icons/gtk-save-as.png"),"Save As",self)
        self.saveAsButton = saveAsButton
        saveAsButton.setObjectName("saveAs")
        saveAsButton.setToolTip("Save current file in own path with own name (advisable: IbdSettings.txt)")
        saveAsButton.setMaximumSize(QSize(85, 22))
        #saveButton.setMinimumSize(QSize(16, 18))
        saveAsButton.setFlat(True)
        QObject.connect(saveAsButton,SIGNAL("clicked()"),self.saveAs)
        self.ui.toolBar.addWidget(saveAsButton)

        runIBDsimButton = QPushButton(QIcon("./uis/icons/IBDSIM.png"),"Run IBDSim",self)
        self.runIBDsimButton = runIBDsimButton
        runIBDsimButton.setObjectName("runIBDsimButton")
        runIBDsimButton.setToolTip("Run IBDSim programme with the current file settings")
        runIBDsimButton.setMaximumSize(QSize(110, 22))
        #saveButton.setMinimumSize(QSize(16, 18))
        runIBDsimButton.setFlat(True)
        QObject.connect(runIBDsimButton,SIGNAL("clicked()"),self.runIBDSim)
        self.ui.toolBar.addWidget(runIBDsimButton)
        #self.ui.toolBar.addSeparator()

# bouton pour cacher la prévisualisation du fichier de config
        self.setVisiblePreviewFile = False
        hideConfFileButton = QPushButton(QIcon("./uis/icons/show_frame.png"),"Hide/Show Preview",self)
        self.hideConfFileButton = hideConfFileButton
        hideConfFileButton.setObjectName("hideConfFileButton")
        hideConfFileButton.setToolTip("Hide/Show Preview configuration file\nWarning to see the errors, must be show")
        #hideConfFileButton.setMaximumSize(QSize(85, 22))
        #hideConfFileButton.setMinimumSize(QSize(16, 18))
        hideConfFileButton.setFlat(True)
        QObject.connect(hideConfFileButton,SIGNAL("clicked()"),self.hidePreviewFileConf)
        self.ui.toolBar.addWidget(hideConfFileButton)

        exitButton = QPushButton(QIcon("./uis/icons/window-close.png"),"Exit",self)
        self.exitButton = exitButton
        exitButton.setObjectName("exitButton")
        exitButton.setToolTip("Exit")
        exitButton.setMaximumSize(QSize(85, 22))
        #saveButton.setMinimumSize(QSize(16, 18))
        exitButton.setFlat(True)
        QObject.connect(exitButton,SIGNAL("clicked()"),self.close)
        self.ui.toolBar.addWidget(exitButton)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.ui.toolBar.addWidget(spacer)

        wtButton = QPushButton(QIcon("./uis/icons/whats.png"),"What's this ?",self)
        wtButton.setDisabled(False)
        wtButton.setToolTip("Click on this button and then on another object to get the documentation")
        self.wtButton = wtButton
        wtButton.setMaximumSize(QSize(136, 22))
        #wtButton.setMinimumSize(QSize(16, 18))
        wtButton.setFlat(True)
        QObject.connect(wtButton,SIGNAL("clicked()"),self.enterWhatsThisMode)
        self.ui.toolBar.addWidget(wtButton)
        #wtButton.hide()

        for but in [saveButton,saveAsButton,wtButton,runIBDsimButton, exitButton, hideConfFileButton]:
            but.setStyleSheet("QPushButton:hover { background-color: #FFD800;  border-style: outset; border-width: 1px; border-color: black;border-style: outset; border-radius: 5px; } QPushButton:pressed { background-color: #EE1C17; border-style: inset;} ")

            # Table Demo Param
        QObject.connect(self.ui.newDemoColPushButton,SIGNAL("clicked()"),self.insertColumn)
        QObject.connect(self.ui.saveDemoColPushButton,SIGNAL("clicked()"),self.saveColumn)
        QObject.connect(self.ui.deleteDemoColPushButton,SIGNAL("clicked()"),self.deleteColumn)
        QObject.connect(self.ui.modifyDemoColPushButton,SIGNAL("clicked()"),self.modifyDemoColumn)
            # Graph TAB
        if ACTIVE_GRAPH == True:
            QObject.connect(self.ui.viewGraphPushButton,SIGNAL("clicked()"),self.activeGraph)
            QObject.connect(self.ui.hideGraphPushButton,SIGNAL("clicked()"),self.hideGraph)

    def openDoc(self):
        """ permet d'ouvrir le fichier PDF d'aide"""
# Pour l'utilisation sous Linux:
        if "linux" in sys.platform:
            os.system('evince ./docs/%s' % namePDF)
# Pour l'utilisation sous windows:
        if "win" in sys.platform and "darwin" not in sys.platform:
            os.startfile('.\docs\%s'% namePDF) 
# Pour l'utilisation sous Mac OS:
        if "darwin" in sys.platform:
            os.system('open ./docs/%s ' % namePDF)

    def createdDocForWhatIsThis(self):
        """ """
        try:
            self.documentator = Documentator("docs/documentation/index.html")
            self.updateDoc()
        except Exception as e:
            log(1,"Documentation error : %s"%e)
            self.documentator = None

    def enterWhatsThisMode(self):
        """ Change le style du curseur de souris et attend un clic
        sur un objet pour afficher son "whatsThis"
        """
        QWhatsThis.enterWhatsThisMode()

    def updateDoc(self,obj=None):
        """ Met à jour le "what's this" pour tous les
        objets fils de obj documentés dans Documentator
        """
        l = []
        if obj == None:
            obj = self
#            log(3,"Updating documentation of "+"%s"%obj.objectName())
            for typ in [QLineEdit,QTableWidgetItem,QPushButton,QPlainTextEdit,QListWidget,QLineEdit,QRadioButton,QComboBox,QProgressBar,QCheckBox]:#,QLabel
                l.extend(obj.findChildren(typ))
        else:
            l.append(obj)
        for e in l:
            if isinstance(e,QTableWidgetItemAbstract) :
                objnamestr = e.name
                objname_debug = 1
            else:
                objnamestr = "%s\n\n"%e.objectName()
                objname_debug = 1
            try:
                if isinstance(e,QTableWidgetItemAbstract):
                    doc_dico = self.documentator.getDocHashByTags("%s" %e.name)
                else:
                    doc_dico = self.documentator.getDocHashByTags("%s"%e.objectName())

                # on remplace les SRC pour que les images soient bien chargées
                for tag in doc_dico.keys():
                    doc_dico[tag] = doc_dico[tag].replace('SRC="','SRC="%s/'%dataPath.DOCPATH)
                docstr = ""
                # on n'encadre pas le default tag
                if doc_dico.has_key("default_tag"):
                    docstr += "<table border='1'><tr><td>" +doc_dico["default_tag"]+"</td></tr></table>\n"
                # on encadre les textes pour chaque tag
                for tag in doc_dico.keys():
                    if tag != "default_tag":
                        docstr += "<table border='1'><tr><th><font color='orange'>%s</font>"%tag
                        docstr += "</th></tr><tr><td>" + doc_dico[tag]+"</td></tr></table>\n"
                if objname_debug:
                    e.setWhatsThis(output.whatsthis_header + objnamestr + "<br/>" + docstr + output.whatsthis_footer)
                else:
                    e.setWhatsThis(output.whatsthis_header + docstr + output.whatsthis_footer)
                
#                if isinstance(e,QTableWidgetItemAbstract):
#                    log(4,"Adding documentation of "+"%s"%e.name)
#                else:
#                    log(4,"Adding documentation of "+"%s"%e.objectName())
            except Exception: # as ex:
#                print ex
                if objname_debug:
                    e.setWhatsThis(output.whatsthis_header + objnamestr + output.whatsthis_footer)
                pass
#    pass

    def createDemograficTable(self) :
        """Méthode pour la création de la table des paramétres démographiques:
        Créer un QTableWidget chaque colonne correspondant à une phase démographique
        Les noms des lignes sont créer à partir de la liste listKeyDemographicParam """
        # Création du QTableWidget
        self.ui.demographicTableWidget = QtGui.QTableWidget(self.ui.tableDemoFrame)
        self.ui.demographicTableWidget.setObjectName("demographicTableWidget")
        self.ui.demographicTableWidget.setColumnCount(0)
        nbrow = len(self.ibdsimConfig.listParametersDemographic)
        self.ui.demographicTableWidget.setRowCount(nbrow)
        row = 0

        for key in self.listKeyDemographicParam:
            # Ajout du label des lignes

            exec('self.ui.%sItem = QTableWidgetItemAbstract(content = "%s", name = "%sItem")' % (key,key,key))
            eval('self.updateDoc(obj = self.ui.%sItem)' % key)
            
            self.ui.demographicTableWidget.setVerticalHeaderItem(row,eval("self.ui.%sItem" % key))
            row = row + 1


        # Applique le QTableWidget à la frame
        self.ui.gridLayoutDemo.addWidget(self.ui.demographicTableWidget, 0, 0, 1, 1)
        self.ui.verticalLayoutDemo.addWidget(self.ui.tableDemoFrame)
        self.ui.demographicTableWidget.show()
        
        # appel fonction pour ajouter la première colonne
#        self.ui.newDemoColPushButton.emit(SIGNAL("clicked()"))
        self.insertColumn()

        QObject.connect(self.ui.demographicTableWidget,SIGNAL('currentCellChanged(int, int, int, int)'),self.saveInQuitCell)
#        QObject.connect(self.ui.demographicTableWidget,SIGNAL('itemChanged()'),self.saveInQuitCell)
#        QObject.connect(self.ui.demographicTableWidget,SIGNAL('cellChanged()'),self.saveInQuitCell)
        # Pour les bouttons de l'interface
        self.ui.newDemoColPushButton.setEnabled (False)
        self.ui.saveDemoColPushButton.setEnabled (True)
        self.ui.deleteDemoColPushButton.setEnabled (False)
        self.ui.modifyDemoColPushButton.setEnabled (False)

    def deleteColumn(self):
        """Méthode qui permet de supprimer la colonne du tableau et les valeurs dans DictIBDSimConfig"""
        # récupère le nombre de colonne et supprime la derniere
        QObject.disconnect(self.ui.demographicTableWidget,SIGNAL('currentCellChanged(int, int, int, int)'),self.saveInQuitCell)
        colAtRemove = self.ui.demographicTableWidget.columnCount()
        self.ui.demographicTableWidget.removeColumn(colAtRemove-1)
        self.ibdsimConfig.deleteValuesInDict(deletedNumber = colAtRemove-1, useList= "listParametersDemographic")

        self.checkCompatibilityOfOption()
        # pour griser les boutons si 1 colonne
        if colAtRemove-1 == 1:
            self.ui.deleteDemoColPushButton.setEnabled (False)
        self.displayParamText()
        QObject.connect(self.ui.demographicTableWidget,SIGNAL('currentCellChanged(int, int, int, int)'),self.saveInQuitCell)

    def modifyDemoColumn(self):
        """Permet la modification des valeurs entrer dans le tableau"""
        self.ui.demographicTableWidget.setEnabled(True)
        self.ui.newDemoColPushButton.setEnabled (False)
        self.ui.saveDemoColPushButton.setEnabled (True)
        self.ui.deleteDemoColPushButton.setEnabled (False)
        self.ui.modifyDemoColPushButton.setEnabled (False)

    def insertColumn(self):
        """Méthode qui permet l'ajout d'une phase démographique: ajoute une colonne dans le tableau
        Utilise le dictionnaire dicoForTypeInputKeys pour connaitre le type d'objet attendu"""
        nbcol = self.ui.demographicTableWidget.columnCount()
        # regarde le nombre de colonne présente dans le tableau puis ajoute une colonne aprés et les compossant
        self.ui.demographicTableWidget.insertColumn(nbcol)
        row = 0
        for dico in self.ibdsimConfig.listParametersDemographic:
            for key, value in dico.items():
                # Ajout des Combobox
                if key in self.dicoForTypeInputKeys.keys():
                    if self.dicoForTypeInputKeys[key] == "CheckBox":
                        exec("self.ui.%sCheckBox%s = QCheckBox()" % (key, nbcol))   # Ajoute le nom a l'ui
                        self.ui.demographicTableWidget.setCellWidget(row,nbcol ,eval("self.ui.%sCheckBox%s" % (key, nbcol) ) )  #ajoute au tableau
                        # creer le signal pour acutaliser
                        exec(  "QObject.connect(self.ui.%sCheckBox%s,SIGNAL('toggled(bool)'),self.actualizeDemoTable)" %(key, nbcol)  )

                    if self.dicoForTypeInputKeys[key] == "Combo": # idem au dessus mais pour les combobox
                        exec("self.ui.%sCombo%s = ComboBoxAbstract(contents=value.expectedValue)" % (key, nbcol))
                        self.ui.demographicTableWidget.setCellWidget(row,nbcol ,eval("self.ui.%sCombo%s" % (key, nbcol) ) )
                        exec(  "QObject.connect(self.ui.%sCombo%s,SIGNAL('currentIndexChanged(int)'),self.actualizeDemoTable)" %(key, nbcol)  )

                else:   # c'est des string attendu
                    try:    # récupere la valeur précédente dnas le dictionnaire ibdsim pour ajouter dans la nouvelle colonne
                        self.ibdsimConfig[key][nbcol-1]
                        if self.ibdsimConfig[key][nbcol-1] != None:
                            self.ui.demographicTableWidget.setItem(row,nbcol ,QtGui.QTableWidgetItem(self.ibdsimConfig[key][nbcol-1]))
                        else:
                            self.ui.demographicTableWidget.setItem(row,nbcol ,QtGui.QTableWidgetItem(""))
                    except:
                        self.ui.demographicTableWidget.setItem(row,nbcol ,QtGui.QTableWidgetItem("%s" % value.defaultValue))
                row = row + 1
        # modifie la taille pour que tout rentre
        self.ui.demographicTableWidget.setColumnWidth(nbcol, 160)
        self.actualizeDemoTable()
        
                #gestion de l'interface
        self.ui.attentionDemographicFrame.show()
        self.ui.newDemoColPushButton.setEnabled (False)
        self.ui.saveDemoColPushButton.setEnabled (True)
        self.ui.deleteDemoColPushButton.setEnabled (False)
        self.ui.modifyDemoColPushButton.setEnabled (False)
        self.ui.demographicTableWidget.setEnabled(True)

    def displayError(self, error):
        """ affiche les erreurs dans la zone de text"""
        self.ui.displayErrorEdit.show()
        self.ui.displayErrorEdit.setPlainText(error)

    def setErrorCell(self, itempass = None):
        """Methode qui permet mettre la case en rouge en cas d'erreur"""
        brush = QtGui.QBrush(QtGui.QColor(255, 50, 50))
        brush.setStyle(QtCore.Qt.SolidPattern)
        itempass.setBackground(brush)

    def setErrorLineEdit(self, lineEdit = None):
        """Methode qui permet mettre les LineEdit en rouge en cas d'erreur"""
        lineEdit.setStyleSheet("background-color: rgb(255, 50, 50);")


    def setDisabledCell(self, itempass = None):
        """Methode qui permet de griser la case et de la rendre non editable, efface la valeur présente"""
        brush = QtGui.QBrush(QtGui.QColor(223, 223, 223))
        brush.setStyle(QtCore.Qt.SolidPattern)
        itempass.setBackground(brush)
        itempass.setText('')
        itempass.setFlags(QtCore.Qt.ItemIsEnabled)

    def setEnabledCell(self, itempass = None):
        """Methode qui permet de degriser la case et la rendre editable"""
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        itempass.setBackground(brush)
        itempass.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable)

    def actualizeDemoTable(self):
        """pour l'affichage du tableau: Parcour tout le tableau pour dégriser les cases,
        puis regarde l'état de chaque Combo et CheckBox pour regriser en appelant leurs fonctions d'actualisations spécifiques
        Cette méthode et appeler à chaque modification dans le tableau"""
        col = 0
        while col < self.ui.demographicTableWidget.columnCount():
            for key in self.listKeyDemographicParam:
                if key not in self.dicoForTypeInputKeys.keys():
                    row = self.listKeyDemographicParam.index(key)
                    self.setEnabledCell(itempass=self.ui.demographicTableWidget.item(row,col))      # Dégrise toute les cases

            # Pour imiter le signal et actualiser les case a griser
            for key in self.dicoForTypeInputKeys.keys():
                if self.dicoForTypeInputKeys[key] == "Combo":
                    eval("self.actualize%s(colonne = col)" % key)
                if self.dicoForTypeInputKeys[key] == "CheckBox":
                    eval( "self.actualize%s(colonne = col)" % key)


            if eval("self.ui.ZoneCheckBox%s.isChecked()" % col) or eval("self.ui.Specific_Density_DesignCheckBox%s.isChecked()" % col):
                eval("self.ui.Continuous_Deme_Size_VariationCombo%s.setEnabled (False)" % col)
            else:
                eval("self.ui.Continuous_Deme_Size_VariationCombo%s.setEnabled (True)" % col)



            # pour la gestion de parametre non editable a la 1ere colonne ou colonne suivante:
            if col == 0:
                for key in self.listKeyDemographicParam:
                    if key in ['New_Demographic_Phase_At', 'Random_Translation']:
                        row = self.listKeyDemographicParam.index(key)
                        self.ui.demographicTableWidget.removeCellWidget(row,col)     # on supprime les widget
                        newitem=QtGui.QTableWidgetItem('')
                        self.ui.demographicTableWidget.setItem(row,col ,newitem)
                        self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,col))
            else:
                for key in self.listKeyDemographicParam:
                    if key in ['Total_Range_Dispersal','Lattice_Boundaries', 'Specific_Sample', 'Sample_SizeX', 'Sample_SizeY', 'Min_Sample_CoordinateX',\
                               'Min_Sample_CoordinateY', 'Sample_Coordinates_X', 'Sample_Coordinates_Y', 'Ind_Per_Pop_Sampled', 'Void_Sample_Node']:
                        row = self.listKeyDemographicParam.index(key)
                        self.ui.demographicTableWidget.removeCellWidget(row,col)     # on supprime les widget
                        newitem=QtGui.QTableWidgetItem('')
                        self.ui.demographicTableWidget.setItem(row,col ,newitem)
                        self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,col))
            col += 1


    def actualizeLattice_Boundaries(self, colonne = None):
        """   """
        pass

    def actualizeTotal_Range_Dispersal(self, colonne = None):
        """   """
        pass

    def actualizeSpecific_Density_Design(self, colonne = None):
        """   """
        if eval("self.ui.Specific_Density_DesignCheckBox%s.isChecked()" % colonne):
            self.displayError(error = "Be careful with \"Specific_Density_Design\" you must have the file \"DensityMatrix.txt\" at path: %s" % os.getcwd())
        else:
            self.displayParamText()
            
    def actualizeRandom_Translation(self, colonne = None):
        """ """
        pass

    def actualizeImmigration_Control(self, colonne = None):
        """   """
        pass

    def actualizeSpecific_Sample(self, colonne = None):
        """   """
        colonne = 0
        lRowisCheck = []
        lRowisUnCheck = []
        
        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled["Specific_SampleCheck"]:
                lRowisCheck.append(self.listKeyDemographicParam.index(key))
            if key in self.actualizeCellDisabled["Specific_SampleUnCheck"]:
                lRowisUnCheck.append(self.listKeyDemographicParam.index(key))
                
        if eval("self.ui.Specific_SampleCheckBox%s.isChecked()" % colonne):
            for row in lRowisCheck:
                self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))
        else:
            for row in lRowisUnCheck:
                self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))

    def actualizeDispersal_Distribution(self, colonne = None):
        """Permet d'actualiser les cases a griser en fonction du contenu de la combobox Dispersal Distribution"""
        lRowDisabled = []
        currentValue = str(eval("self.ui.Dispersal_DistributionCombo%s.currentText()" % colonne))

        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled[currentValue]:
                lRowDisabled.append(self.listKeyDemographicParam.index(key))

            if key in self.valueForCell[currentValue].keys():
                row = self.listKeyDemographicParam.index(key)
                self.ui.demographicTableWidget.removeCellWidget(row,colonne)     # on supprime les widget
                newitem=QtGui.QTableWidgetItem(str(self.valueForCell[currentValue][key]))
                self.ui.demographicTableWidget.setItem(row,colonne ,newitem)


        for row in lRowDisabled:
            self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))

    def actualizeZone(self, colonne = None):
        """ Permet d'actualiser les cases a griser en fonction de la checkbox zone"""
        lRowDisabled = []
        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled["Zone"]:
                lRowDisabled.append(self.listKeyDemographicParam.index(key))
        if eval("self.ui.ZoneCheckBox%s.isChecked()" % colonne):
            for row in lRowDisabled:
                self.setEnabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))
        else:
            for row in lRowDisabled:
                self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))

    def actualizeBarrier(self, colonne = None):
        """ Permet d'actualiser les cases a griser en fonction de la checkbox zone"""
        lRowDisabled = []
        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled["Barrier"]:
                lRowDisabled.append(self.listKeyDemographicParam.index(key))
        if eval("self.ui.BarrierCheckBox%s.isChecked()" % colonne):
            for row in lRowDisabled:
                self.setEnabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))
        else:
            for row in lRowDisabled:
                self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))


    def actualizeContinuous_Deme_Size_Variation(self, colonne = None):
        """ permet d'actualiser les cases a griser en fonction du contenu de la combobox ContinuousDemeSizeVariantion """
        lRowDisabled = []
        currentValue = str(eval("self.ui.Continuous_Deme_Size_VariationCombo%s.currentText()" % colonne))
        
        if currentValue in ['Linear', 'Exponential']:
            self.displayError(error = "Be careful with \"%s\" in \"Continuous_Deme_Size_Variation\" you need to do a new demographic change" % currentValue)

        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled[currentValue]:
                lRowDisabled.append(self.listKeyDemographicParam.index(key))

        for row in lRowDisabled:
            self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))

    def actualizeContinuous_Lattice_Size_Variation(self, colonne = None):
        """ permet d'actualiser les cases a griser en fonction du contenu de la combobox ContinuousLatticeSizeVariantion """
        lRowDisabled = []
        currentValue = str(eval("self.ui.Continuous_Lattice_Size_VariationCombo%s.currentText()" % colonne))
        
        if currentValue in ['Linear', 'Exponential']:
            self.displayError(error = "Be careful with \"%s\" in \"Continuous_Lattice_Size_Variation\" you need to do a new demographic change" % currentValue)

        for key in self.listKeyDemographicParam:
            if key in self.actualizeCellDisabled2[currentValue]:
                lRowDisabled.append(self.listKeyDemographicParam.index(key))

        for row in lRowDisabled:
            self.setDisabledCell(itempass=self.ui.demographicTableWidget.item(row,colonne))

    def saveInQuitCell(self):
        """Fonction qui permet de vérifier toute les données lorsque l'utilisateur quite une cellule"""
        try:
            self.ibdsimConfig.enterInSessionMode()
            self.ui.displayErrorEdit.hide()
            col = 0
            # parcour les colonnes puis ligne par ligne
            while col < self.ui.demographicTableWidget.columnCount():
                row = 0
                for key in self.listKeyDemographicParam:
                    if isinstance(self.ui.demographicTableWidget.cellWidget(row, col), ComboBoxAbstract):
                        self.ibdsimConfig[key].append( str(eval("self.ui.%sCombo%s.currentText()" % (key, col)) ) )
                    elif isinstance(self.ui.demographicTableWidget.cellWidget(row, col), QCheckBox):
                        if eval("self.ui.%sCheckBox%s.checkState()" % (key,col)) == 0:
                            self.ibdsimConfig[key].append(value = "False" )
                        else:
                            self.ibdsimConfig[key].append(value = "True" )
                    else:
                        self.ibdsimConfig[key].append( str(self.ui.demographicTableWidget.item(row, col).text()) )
                    row += 1
                col += 1
            self.ibdsimConfig._rollback()
            self.ui.paramTextEdit.setPlainText("")
            self.actualizeDemoTable()
            self.displayParamText()

        except ValueError as infos:
            self.actualizeDemoTable()
            self.setErrorCell(itempass=self.ui.demographicTableWidget.item(row,col))
            self.ibdsimConfig._rollback()
            self.displayError(error = "Value Error for \"%s\" : %s " % (key , infos))

    def saveColumn(self):
        """Méthode qui permet de lire les valeurs entrer dans les colonnes par l'utilisateur et de les stocker dans le dictionnaire si correcte
        vérifie toute les colonnes, en écrasant les valeurs déja stocker et vérifiant celle en cour dans l'interface"""
        self.actualizeDemoTable()
        self.ui.displayErrorEdit.hide()
        try:
            self.ibdsimConfig.enterInSessionMode()
            # Vide les paramètres deja entrer
            self.ibdsimConfig.resetDemographicParam()
            col = 0
            # parcour les colonnes puis ligne par ligne
            while col < self.ui.demographicTableWidget.columnCount():
                row = 0
                for key in self.listKeyDemographicParam:
                    if isinstance(self.ui.demographicTableWidget.cellWidget(row, col), ComboBoxAbstract):
                        self.ibdsimConfig[key].append( str(eval("self.ui.%sCombo%s.currentText()" % (key, col)) ) )
                    elif isinstance(self.ui.demographicTableWidget.cellWidget(row, col), QCheckBox):
                        if eval("self.ui.%sCheckBox%s.checkState()" % (key,col)) == 0:
                            self.ibdsimConfig[key].append(value = "False" )
                        else:
                            self.ibdsimConfig[key].append(value = "True" )
                    else:
                        self.ibdsimConfig[key].append( str(self.ui.demographicTableWidget.item(row, col).text()) )
                    row += 1
                col += 1    

            self.checkCompatibilityOfOption()
            self.ibdsimConfig.commitIbdsimWithNone(useList = "listParametersDemographic")
            self.displayParamText()

            # gestion des boutons interface
            self.ui.attentionDemographicFrame.hide()
            self.ui.newDemoColPushButton.setEnabled (True)
            self.ui.saveDemoColPushButton.setEnabled (False)
            self.ui.modifyDemoColPushButton.setEnabled (True)
            if self.ui.demographicTableWidget.columnCount() == 1:
                self.ui.deleteDemoColPushButton.setEnabled (False)
            else:
                self.ui.deleteDemoColPushButton.setEnabled (True)
            self.ui.demographicTableWidget.setEnabled(False)    

        except ValueError as infos:
            print infos
            self.actualizeDemoTable()
            self.setErrorCell(itempass=self.ui.demographicTableWidget.item(row,col))
            self.ibdsimConfig._rollback()
            self.displayError(error = "Value Error for \"%s\" : %s " % (key , infos))

        except AttributeError as inst:
            self.actualizeDemoTable()
            infos, keys, col = inst.args
            for key in keys:
                self.setErrorCell(itempass=self.ui.demographicTableWidget.item(int(self.listKeyDemographicParam.index(key)),col))
            self.ibdsimConfig._rollback()
            self.displayError(error = "%s" %  infos)

    def actualizeModel(self):
        """actualise l'affichage de l'ui en fonction des parametres courants de l'interface"""
        #*************** Markers Tab**************
        modelc = self.ui.Mutation_ModelCombo.currentText()
        # Regarde le model courant et applique la frame corespondante
        for model in self.ibdsimConfig.allModels:
            eval("self.ui.%sframe.hide()" % model)
        eval("self.ui.%sframe.show()" % modelc)

    def newModel(self):
        """Active la frame pour créer un nouveau model et entrer les valeurs de paramètres"""
        self.actualizeModel()
        self.ui.newModelFrame.show()
        self.__nbModel += 1 

        self.ui.attentionModelFrame.show()
        
        self.ui.editModelCombo.insertItem(self.__nbModel,str(self.__nbModel))
        self.ui.editModelCombo.setFocus( self.__nbModel )                  # a faire fonctionner...
        self.ui.newModelPushButton.setEnabled (False)
        self.ui.saveModelPushButton.setEnabled (True)
        self.ui.deleteModelPushButton.setEnabled (False)
        self.ui.editModelCombo.setEnabled (False)

    def deleteModel(self):
        """suppression du model selectionner"""
        nb = int( self.ui.editModelCombo.currentText() )-1
        self.__nbModel -= 1
        self.ibdsimConfig.deleteValuesInDict(deletedNumber = nb, useList= "listParametersModel")
        self.checkCompatibilityOfOptionModel()
        self.ui.editModelCombo.clear()
        for i in range(self.__nbModel):
            self.ui.editModelCombo.insertItem(i,str(i+1))

        if self.ui.editModelCombo.currentText() == "":
            self.ui.editModelCombo.setEnabled (False)
            self.ui.deleteModelPushButton.setEnabled (False)
        self.displayParamText()

    def saveModel(self):
        """sauvegarde des donnees du model"""
        l = []
        for typ in [QLineEdit,QComboBox,QCheckBox]:# QListWidget,QPlainTextEdit 
            l.extend(self.ui.MarkersTab.findChildren(typ))
        l = [str(e.objectName()) for e in l ]
        modelc = self.ui.Mutation_ModelCombo.currentText()

        self.ui.displayErrorEdit.hide()
        try:
            self.ibdsimConfig.enterInSessionMode()
            for param in self.ibdsimConfig.keys():
                if param+"Edit" in l:
                    lineEdit = "self.ui.%sEdit" % param
                    eval( """%s.setStyleSheet("background-color: rgb(255, 255, 255);") """ % lineEdit)
                    self.ibdsimConfig[param].append(value = str(eval("self.ui.%sEdit.text()" % param )) )

                if param+"%sEdit" % modelc in l:
                    lineEdit = "self.ui.%s%sEdit" % (param, modelc)
                    eval( """%s.setStyleSheet("background-color: rgb(255, 255, 255);")""" % lineEdit)
                    self.ibdsimConfig[param].append(value = str(eval("self.ui.%sEdit.text()" % (param+modelc) )) )

                if param+"Combo" in l:
                    self.ibdsimConfig[param].append(value = str(eval("self.ui.%sCombo.currentText()" % param)) )

                if param+"CheckBox" in l:
                    if eval("self.ui.%sCheckBox.checkState()" % param) == 0:
                        self.ibdsimConfig[param].append(value = "False" )
                    else:
                        self.ibdsimConfig[param].append(value = "True" )

            self.checkCompatibilityOfOptionModel()
            self.ibdsimConfig.commitIbdsimWithNone(useList="listParametersModel")
            self.displayParamText()

        # Gestion de l'affichage
            self.ui.newModelFrame.hide()                        # Cache la frame de selection model
            eval("self.ui.%sframe.hide()" % modelc)        # Cache la frame du model
            self.ui.attentionModelFrame.hide()
            self.ui.newModelPushButton.setEnabled (True)
            self.ui.saveModelPushButton.setEnabled (False)
            self.ui.deleteModelPushButton.setEnabled (True)
            self.ui.editModelCombo.setEnabled (True)
        except AttributeError as inst:
            infos, keys = inst.args
            for key in keys:
                self.setErrorLineEdit(lineEdit = eval("self.ui.%sEdit" % key))
            self.ibdsimConfig._rollback()
            self.displayError(error = "%s" %  infos)
        except ValueError as infos:
            self.ibdsimConfig._rollback()
            try:
                self.setErrorLineEdit(lineEdit = eval("self.ui.%sEdit" % param))
            except:
                self.setErrorLineEdit(lineEdit = eval("self.ui.%s%sEdit" % (param, modelc) ))
            self.displayError(error = "Value Error for \"%s\" : %s " % (param , infos))



    def displayParamText(self):
        """Affichage des paramètres dans la zone de text de l'interface"""
        self.ui.paramTextEdit.setPlainText(self.ibdsimConfig.__repr__())

    def createComputationListe(self):
        """ Méthode qui lit la liste des "Various Computation Options" et ajoute les valeurs dans la liste """
        for item in self.ibdsimConfig['DiagnosticTables'].expectedValue:
            self.ui.availableCompListWidget.addItem(item)
            
    def addCompItem(self):
        """Methode qui permet d'ajouter un élément de la liste des "Various Computation Options" """
        row = self.ui.availableCompListWidget.currentRow() 
        item = self.ui.availableCompListWidget.takeItem(row).text()
        self.ui.selectListCompWidget.addItem(item)

        nbAvailable = self.ui.availableCompListWidget.count()
        if nbAvailable == 0:
            self.ui.selectCompPushButton.setEnabled(False)
        else:
            self.ui.selectCompPushButton.setEnabled(True)
        nbSelect = self.ui.selectListCompWidget.count()
        if nbSelect == 0:
            self.ui.deleteSelectParamCompPushButton.setEnabled(False)
        else:
            self.ui.deleteSelectParamCompPushButton.setEnabled(True)

        self.ibdsimConfig["DiagnosticTables"].append(str(item))
        self.displayParamText()

    def deleteCompItem(self):
        """Methode qui permet de supprimer un élément de la liste des "Various Computation Options" """
        row = self.ui.selectListCompWidget.currentRow() 
        item = self.ui.selectListCompWidget.takeItem(row).text()
        self.ui.availableCompListWidget.addItem(item)

        nbAvailable = self.ui.availableCompListWidget.count()
        if nbAvailable == 0:
            self.ui.selectCompPushButton.setEnabled(False)
        else:
            self.ui.selectCompPushButton.setEnabled(True)
        nbSelect = self.ui.selectListCompWidget.count()
        if nbSelect == 0:
            self.ui.deleteSelectParamCompPushButton.setEnabled(False)
        else:
            self.ui.deleteSelectParamCompPushButton.setEnabled(True)
        
        self.ibdsimConfig["DiagnosticTables"].remove(item)
        self.displayParamText()

    def createdConnectionSimulation(self):
        """ Ajoute les connections pour les objets présents dans l'interface à l'onget simulation and output"""
        l = []
        for typ in [QLineEdit,QComboBox,QCheckBox]:# QListWidget,QPlainTextEdit 
            l.extend(self.ui.Simulationtab.findChildren(typ))
        self.objInSimulationTab = [str(e.objectName()) for e in l ]
        
        for item in self.objInSimulationTab:
            item = str(item)
            if item.find("Edit") != -1:
                eval("""QObject.connect(self.ui.%s,SIGNAL("editingFinished()"),self.saveSimulation)""" % (item))
            if item.find("CheckBox") != -1:
                exec("""QObject.connect(self.ui.%s,SIGNAL("toggled(bool)"),self.saveSimulation)""" % (item))
            if item.find("Combo") != -1:
                exec("""QObject.connect(self.ui.%s,SIGNAL("currentIndexChanged(int)"),self.saveSimulation)""" % (item))

    def saveSimulation(self):
        """ Sauvegarde les paramètres de la tabulation Simulation"""
        self.ibdsimConfig.deleteValuesInDict(deletedNumber= 0, useList = "listParametersSimulAndOutput")
        self.actualizeNexusCheckBox()
        self.ui.displayErrorEdit.hide()
        try:
            self.ibdsimConfig.enterInSessionMode()
            for param in self.ibdsimConfig.keys():
                if param+"Edit" in self.objInSimulationTab:
                    lineEdit = eval("self.ui.%sEdit" % param)
                    lineEdit.setStyleSheet("background-color: rgb(255, 255,255);")
                    self.ibdsimConfig[param].append(value = str(eval("self.ui.%sEdit.text()" % param )) )

                elif param+"Combo" in self.objInSimulationTab:
                    if (param == "Nexus_File_Format"):
                        if (self.ui.NexusCheckBox.isChecked()):
                            self.ibdsimConfig[param].append(value = str(eval("self.ui.%sCombo.currentText()" % param)) )
                    if (param != "Nexus_File_Format"):
                        self.ibdsimConfig[param].append(value = str(eval("self.ui.%sCombo.currentText()" % param)) )

                elif param+"CheckBox" in self.objInSimulationTab:
                    if eval("self.ui.%sCheckBox.checkState()" % param) == 0:
                        self.ibdsimConfig[param].append(value = "False" )
                    else:
                        self.ibdsimConfig[param].append(value = "True" )

            self.ibdsimConfig.commitIbdsimWithNone(useList = "listParametersSimulAndOutput")
            self.displayParamText()
        except ValueError as infos:
            self.ibdsimConfig._rollback()
            self.setErrorLineEdit(lineEdit = eval("self.ui.%sEdit" % param))
            self.displayError(error = "Value Error for \"%s\" : %s " % (param , infos))

    def actualizeNexusCheckBox(self):
        """ Affiche la frame nexus"""
        if self.ui.NexusCheckBox.isChecked():
            self.ui.nexusFormatFrame.show()
        else:
            self.ui.nexusFormatFrame.hide()

    def saveFile(self):
        """Méthode qui permet de sauvegarder la zone de text des paramètres dans un fichier"""
# Pour l'utilisation sous Mac OS:
        if "darwin" in sys.platform :
            filename = "../../../IbdSettings.txt"
    # Pour l'utilisation sous Linux:   
        if "linux" in sys.platform:
            filename = "./IbdSettings.txt"
    # Pour l'utilisation sous windows:
        if "win" in sys.platform and "darwin" not in sys.platform:
            filename = ".\IbdSettings.txt"
    
        f = open(filename, 'w') 
        filedata = self.ui.paramTextEdit.toPlainText() 
        f.write(filedata)
        f.close()
        if "darwin" in sys.platform :
            pathSave = sys.argv[0].split("GUI_IBDSim")[0][:0]
            self.ui.statusbar.showMessage (QString("The file is save in %s" % pathSave +"/ with name IbdSettings.txt"), 7200)
        else:
            self.ui.statusbar.showMessage (QString("The file is save in %s" % os.getcwd()+"/ with name IbdSettings.txt"), 7200)

    def saveAs(self):
        """Méthode qui permet de sauvegarder la zone de text des paramètres dans un fichier"""
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save the File", os.getcwd(), "Text file *.txt","Text file *.txt")
        f = open(filename, 'w') 
        filedata = self.ui.paramTextEdit.toPlainText() 
        f.write(filedata) 
        f.close()
        self.ui.statusbar.showMessage (QString("The file is save in %s" % filename),7200)


    def runIBDSim(self):
        """  lance IBDSim avec les paramètres de l'interfaces"""
        self.ui.statusbar.showMessage (QString("IBDSim is running, please view terminal"), 0)
# Pour l'utilisation sous Linux:
        if "linux" in sys.platform :
            filename = "./IbdSettings.txt"
            f = open(filename, 'w') 
            filedata = self.ui.paramTextEdit.toPlainText() 
            f.write(filedata)
            f.close()
            Popen(args= ["./IBD_Sim"],shell=True)
            self.ui.statusbar.showMessage (QString("IBDSim has finished running "), 7200)
            

# Pour l'utilisation sous windows:
        if "win" in sys.platform and "darwin" not in sys.platform :
            filename = "./IbdSettings.txt"
            f = open(filename, 'w') 
            filedata = self.ui.paramTextEdit.toPlainText() 
            f.write(filedata)
            f.close()
            Popen(args= ["IBDSim.exe"],shell=True)
            self.ui.statusbar.showMessage (QString("IBDSim has finished running "), 7200)

# Pour l'utilisation sous Mac OS:
        if "darwin" in sys.platform :
            filename = "../../../IbdSettings.txt"
            f = open(filename, 'w') 
            filedata = self.ui.paramTextEdit.toPlainText() 
            f.write(filedata)
            f.close()
            Popen(args= ["../../../IBD_Sim"],shell=True)
            self.ui.statusbar.showMessage (QString("IBDSim has finished running "), 7200)
        

#***********************************************************************************************************************
    def checkCompatibilityOfOptionModel(self):
        """ Fonction qui vérifie que les valeurs entrer sont bien correcte """
#        pour toutes les catégories/set de marqueurs vérifier que :
#        les equilibrium fréquencies (Equilibrium_Frequencies) doivent etre positive ou nulle
#        la somme doit faire 1.0        -----> Réglé dans l'objet

#        la longueur des séquences simulées (Sequence_Size) doit etre la meme que celle de la sequence du MRCA (MRCA_Sequence)
        for i in range(len(self.ibdsimConfig["MRCA_Sequence"])):
            if self.ibdsimConfig["MRCA_Sequence"][i] != None:
                model = self.ibdsimConfig["Mutation_Model"][i]
                if len(self.ibdsimConfig["Sequence_Size"][i]) == 0:
                    raise AttributeError("Please input a sequence for Sequence_Size",["Sequence_Size%s" % model])
                if len(self.ibdsimConfig["MRCA_Sequence"][i]) != 0:
                    if int( len(self.ibdsimConfig["MRCA_Sequence"][i])) != int( self.ibdsimConfig["Sequence_Size"][i]):
                        raise AttributeError("The specified \"Sequence_Size\" and \"MRCA_Sequence\" do not correspond.",["MRCA_Sequence%s" % model,"Sequence_Size%s" % model])


#        
#        pour des marqueurs alléliques (KAM, SMM, GSM, TPM),
#        Allelic_State_MRCA doit etre compris entre Allelic_Lower_Bound et Allelic_Upper_Bound """"""si diff 0
        for i in range(len(self.ibdsimConfig["Mutation_Model"])):
            model = self.ibdsimConfig["Mutation_Model"][i]
#        POur des marqueurs SNP, minAlleleNumber doit etre 1 ou 2            
            if model in ["SNP"]:
                if self.ibdsimConfig["Min_Allele_Number"][i] not in ['1', '2']:
                    raise AttributeError("With \"SNP\" Mutation Model, Min_Allele_Number must be equal at 1 or 2 only",["Min_Allele_Number"])
                
            
            
            
            if model in ['KAM', 'SMM', 'GSM', 'TPM']:
                if len(self.ibdsimConfig["Allelic_State_MRCA"][i]) == 0:
                    raise AttributeError("Please input a value for Allelic_State_MRCA",["Allelic_State_MRCA%s" % model])
                if len(self.ibdsimConfig["Allelic_Upper_Bound"][i]) == 0:
                    raise AttributeError("Please input a value for Allelic_Upper_Bound",["Allelic_Upper_Bound%s" % model])
                if len(self.ibdsimConfig["Allelic_Lower_Bound"][i]) == 0:
                    raise AttributeError("Please input a value for Allelic_Lower_Bound",["Allelic_Lower_Bound%s" % model])
                if int(self.ibdsimConfig["Allelic_State_MRCA"][i]) != 0:
                    if int(self.ibdsimConfig["Allelic_State_MRCA"][i]) > int(self.ibdsimConfig["Allelic_Upper_Bound"][i]):
                        raise AttributeError("Allelic_State_MRCA\" not between \"Allelic_Upper_Bound\" and \"Allelic_Lower_Bound\"",["Allelic_State_MRCA%s" % model,"Allelic_Upper_Bound%s" % model,"Allelic_Lower_Bound%s" % model])
                    if int(self.ibdsimConfig["Allelic_Lower_Bound"][i]) > int(self.ibdsimConfig["Allelic_State_MRCA"][i]):
                        raise AttributeError("Allelic_State_MRCA\" not between \"Allelic_Upper_Bound\" and \"Allelic_Lower_Bound\"",["Allelic_State_MRCA%s" % model,"Allelic_Upper_Bound%s" % model,"Allelic_Lower_Bound%s" % model])
#        Gestion du Polymorphic_Loci_Only
            if int(self.ibdsimConfig["Min_Allele_Number"][i]) >= 2:
                if self.ibdsimConfig["Polymorphic_Loci_Only"][i] == "False":
                    raise AttributeError("With Min_Allele_Number >= 2 you must check Polymorphic_Loci_Only",["Min_Allele_Number"])
            elif int(self.ibdsimConfig["Min_Allele_Number"][i]) == 1:
                if self.ibdsimConfig["Polymorphic_Loci_Only"][i] == "True":
                    raise AttributeError("With Min_Allele_Number = 1 you can't check Polymorphic_Loci_Only",["Min_Allele_Number"])

#        Minor_Allele_Frequency doit etre entre 0 et 0.5   ------> Dans IBDSIMCONFIG
        # Si >2 model pas de format Migrate
#        if len(self.ibdsimConfig["Mutation_Model"]) > 1:
#            self.ui.MigrateCheckBox.setChecked(False)
#            self.ui.Migrate_LetterCheckBox.setChecked(False)
#            self.ibdsimConfig["Migrate"][0] = "False"
#            self.ibdsimConfig["Migrate_Letter"][0] = "False"
#            self.ui.MigrateCheckBox.setEnabled(False)
#            self.ui.Migrate_LetterCheckBox.setEnabled(False)
#        else:
#            self.ui.MigrateCheckBox.setEnabled(True)
#            self.ui.Migrate_LetterCheckBox.setEnabled(True)
            
            
################## DEMOGRAPHIC Option ##########################################
    def checkCompatibilityOfOption(self):
        """ Fonction qui vérifie que les valeurs entrer sont bien correcte """
        self.ui.demographicTableWidget.setEnabled(True)
        self.ui.newDemoColPushButton.setEnabled (False)
        self.ui.saveDemoColPushButton.setEnabled (True)
        self.ui.modifyDemoColPushButton.setEnabled (False)

#        LatticeSizeX doit etre >= LatticeSizeY (contrainte débile je suis d'accord, mais c'est comme ca!)
        for i in range(len(self.ibdsimConfig["Lattice_SizeX"])):
            if i!=0:
                if len(self.ibdsimConfig["New_Demographic_Phase_At"][i]) == 0:
                    raise AttributeError("Please input a sequence for New_Demographic_Phase_At", ["New_Demographic_Phase_At"], i)
            if len(self.ibdsimConfig["Lattice_SizeX"][i]) == 0:
                raise AttributeError("Please input a sequence for Lattice_SizeX", ["Lattice_SizeX"], i)
            if len(self.ibdsimConfig["Lattice_SizeY"][i]) == 0:
                raise AttributeError("Please input a sequence for Lattice_SizeY", ["Lattice_SizeY"], i)
            if int(self.ibdsimConfig["Lattice_SizeX"][i]) < int(self.ibdsimConfig["Lattice_SizeY"][i]):
                raise AttributeError("Found Lattice_SizeX < Lattice_SizeY", ["Lattice_SizeX","Lattice_SizeY"], i)


#        une "NewDemographicPhase" doit commencer (="NewDemographicPhaseAt=X") apres la précédente (Xn<Xn-1>..>X2>X1>0)
        nbDemografiqueChange = len( self.ibdsimConfig["New_Demographic_Phase_At"])-1
        comp = nbDemografiqueChange
        if nbDemografiqueChange >= 2:
            while comp >= 2:
                if self.ibdsimConfig["New_Demographic_Phase_At"][comp] <= self.ibdsimConfig["New_Demographic_Phase_At"][comp-1]:
                    raise AttributeError("New demographic phase must begin at later time than previous one", ["New_Demographic_Phase_At"], comp)
                comp -= 1


#        Le nombre de coordonnées de l'échantillon en X (taille du vecteur "Sample_Coordinates_X")
#         doit etre égale au nombre de coordonnées de l'échantillon en Y
        if self.ibdsimConfig["Sample_Coordinates_X"][0] != None:
            if len(self.ibdsimConfig["Sample_Coordinates_X"][0].split(' ')) != len(self.ibdsimConfig["Sample_Coordinates_Y"][0].split(' ')):
                raise AttributeError("The number of coordinates of 'Sample_Coordinates_X' is not the same as the number of coordinates of 'Sample_Coordinates_Y'",["Sample_Coordinates_X","Sample_Coordinates_Y"],0)

#        quand l'échantillon est "non specifique" = avec SampleSizeX et Y... on doit vérifier toutes les conditions suivantes :
        if (self.ibdsimConfig["Specific_Sample"][0] == "False") :
            if len(self.ibdsimConfig["Sample_SizeX"][0]) == 0:
                raise AttributeError("Please input a sequence for Sample_SizeX", ["Sample_SizeX"], 0)
            if len(self.ibdsimConfig["Sample_SizeY"][0]) == 0:
                raise AttributeError("Please input a sequence for Sample_SizeY", ["Sample_SizeY"], 0)
            if len(self.ibdsimConfig["Void_Nodes"][0]) == 0:
                raise AttributeError("Please input a sequence for Void_Nodes", ["Void_Nodes"], 0)
            if len(self.ibdsimConfig["Void_Sample_Node"][0]) == 0:
                raise AttributeError("Please input a sequence for Void_Sample_Node", ["Void_Sample_Node"], 0)
            if len(self.ibdsimConfig["Min_Sample_CoordinateX"][0]) == 0:
                raise AttributeError("Please input a sequence for Min_Sample_CoordinateX", ["Min_Sample_CoordinateX"], 0)
            if len(self.ibdsimConfig["Min_Sample_CoordinateY"][0]) == 0:
                raise AttributeError("Please input a sequence for Min_Sample_CoordinateY", ["Min_Sample_CoordinateY"], 0)

#            LatticeSizeX et Y >= (SampleSizeX et Y-1)*void_sampled_nodes+1
            if int(self.ibdsimConfig["Lattice_SizeX"][0]) < ( (int(self.ibdsimConfig["Sample_SizeX"][0])-1) * (int(self.ibdsimConfig["Void_Sample_Node"][0] )) +1  ):
                res = (int(self.ibdsimConfig["Sample_SizeX"][0])-1*(int(self.ibdsimConfig["Void_Sample_Node"][0])+1))
                raise AttributeError("Habitat dimension LatticeSizeX = %s < sample dimension = %s\n as implied by SampleSizeX = %s and Void_Sample_Node = %s"\
                % (self.ibdsimConfig["Lattice_SizeX"][0],res, self.ibdsimConfig['Sample_SizeX'][0],self.ibdsimConfig["Void_Sample_Node"][0] ), ["Lattice_SizeX", "Sample_SizeX", "Void_Sample_Node"],0)
                
            if int(self.ibdsimConfig["Lattice_SizeY"][0]) < ( (int(self.ibdsimConfig["Sample_SizeY"][0])-1) * (int(self.ibdsimConfig["Void_Sample_Node"][0] )) +1  ):
                res = (int(self.ibdsimConfig["Sample_SizeY"][0])-1*(int(self.ibdsimConfig["Void_Sample_Node"][0])+1))
                raise AttributeError("Habitat dimension LatticeSizeY = %s < sample dimension = %s\n as implied by SampleSizeY = %s and Void_Sample_Node = %s"\
                % (self.ibdsimConfig["Lattice_SizeY"][0],res, self.ibdsimConfig['Sample_SizeY'][0],self.ibdsimConfig["Void_Sample_Node"][0] ), ["Lattice_SizeY", "Sample_SizeY", "Void_Sample_Node"],0)

#             Min_Sample_CoordX et Y doivent etre un mutliple de voidNodes
            if (int(self.ibdsimConfig["Min_Sample_CoordinateX"][0]) % int(self.ibdsimConfig["Void_Nodes"][0])) != 0:
                raise AttributeError("The spacing between samples (Min_Sample_CoordinateX) is not a multiple of the spacing between occupied nodes (Void_Nodes)", ["Min_Sample_CoordinateX", "Void_Nodes"],0)
            if (int(self.ibdsimConfig["Min_Sample_CoordinateY"][0]) % int(self.ibdsimConfig["Void_Nodes"][0])) != 0:
                raise AttributeError("The spacing between samples (Min_Sample_CoordinateY) is not a multiple of the spacing between occupied nodes (Void_Nodes)", ["Min_Sample_CoordinateY", "Void_Nodes"],0)

#             LatticeSizeX et Y doivent etre >= à ( (Min_Sample_CoordX et Y) + (DimSampleX et Y-1)*VoidSampleNodes)
#                if int(self.ibdsimConfig["Lattice_SizeX"][0]) < (int(self.ibdsimConfig["Min_Sample_CoordinateX"][0]) + ((int(self.ibdsimConfig["Sample_SizeX"][0])-1)) * int(self.ibdsimConfig["Void_Sample_Node"][0])):
#                    raise AttributeError(""Habitat dimension LatticeSizeX< maximum sample coordinate ="<<xmin_sample+(dim_sample1-1)*vide_sampleX;
#            cerr<<"\n    as implied by MinSampleCoordinateX="<<xmin_sample<<" ,SampleSizeX="<<dim_sample1<<" and void_sample_node=", ["Lattice_SizeX", "Min_Sample_CoordinateX", "Sample_SizeX", "Void_Sample_Node"],0)
            if int(self.ibdsimConfig["Lattice_SizeY"][0]) < (int(self.ibdsimConfig["Min_Sample_CoordinateY"][0]) + ((int(self.ibdsimConfig["Sample_SizeY"][0])-1)) * int(self.ibdsimConfig["Void_Sample_Node"][0])):
                raise AttributeError("error", ["Lattice_SizeY", "Min_Sample_CoordinateY", "Sample_SizeY", "Void_Sample_Node"],0)

#            Min_Sample_CoordX et Y doivent etre >0        -------> Dans DictIBDSimConfig

#            VoidSampleNodes doit etre un multiple de VoidNodes
            if (int(self.ibdsimConfig["Void_Sample_Node"][0]) % int(self.ibdsimConfig["Void_Nodes"][0])) != 0:
                raise AttributeError("The spacing between Void_Sample_Node is not a multiple of the spacing between occupied Void_Nodes",["Void_Sample_Node", "Void_Nodes"], 0)


#        si l'échantillon est spécifique = avec Sample_CoordX et Y, on doit vérifier les conditions suivantes :
#            et doivent etre <= LatticeSize X et Y
        if (self.ibdsimConfig["Specific_Sample"][0] == "True") :
#            toutes les coordonnées en X et Y doivent etre un multiple de voideNodes a toutes les phases demo
            for coorx in self.ibdsimConfig["Sample_Coordinates_X"][0].split(' '):
                if (int(coorx) % int(self.ibdsimConfig["Void_Nodes"][i])) != 0:
                    raise AttributeError("The spacing between Sample_Coordinates_X at position %s is not a multiple of the spacing between occupied Void_Nodes" % self.ibdsimConfig["Sample_Coordinates_X"][i].split(' ').index(coorx, ),["Sample_Coordinates_X", "Void_Nodes"], 0)
                if int(coorx) > int(self.ibdsimConfig["Lattice_SizeX"][0]):
                    raise AttributeError("Sample_Coordinates_X at position %s is greater than Lattice_SizeX" % self.ibdsimConfig["Sample_Coordinates_X"][0].split(' ').index(coorx, ),["Sample_Coordinates_X", "Lattice_SizeX"], 0)
            
            for coory in self.ibdsimConfig["Sample_Coordinates_Y"][0].split(' '):
                if (int(coory) % int(self.ibdsimConfig["Void_Nodes"][0])) != 0:
                    raise AttributeError("The spacing between Sample_Coordinates_Y at position %s is not a multiple of the spacing between occupied Void_Nodes" % self.ibdsimConfig["Sample_Coordinates_Y"][0].split(' ').index(coory, ),["Sample_Coordinates_Y", "Void_Nodes"], 0)
                if int(coory) > int(self.ibdsimConfig["Lattice_SizeY"][0]):
                    raise AttributeError("Sample_Coordinates_Y at position %s is greater than Lattice_SizeY" % self.ibdsimConfig["Sample_Coordinates_Y"][0].split(' ').index(coory, ),["Sample_Coordinates_Y", "Lattice_SizeY"], 0)

#        pour chaque phase demo il faut vérifier que : 
#            Lattice Szie X et Y soient des multiples de voidNodes   
        for i in range( len( self.ibdsimConfig["Void_Nodes"])):
            if (int(self.ibdsimConfig["Lattice_SizeX"][i]) % int(self.ibdsimConfig["Void_Nodes"][i])) != 0:
                raise AttributeError("The spacing between Lattice_SizeX is not a multiple of the spacing between occupied Void_Nodes",["Lattice_SizeX", "Void_Nodes"], i)
            if (int(self.ibdsimConfig["Lattice_SizeY"][i]) % int(self.ibdsimConfig["Void_Nodes"][i])) != 0:
                raise AttributeError("The spacing between Lattice_SizeY is not a multiple of the spacing between occupied Void_Nodes",["Lattice_SizeY", "Void_Nodes"], i)


#        Si l'échantillon est spécifique, alors on ne peut pas calculer les différentes stats de DiagnosticTables, il doit etre vide
#            if (self.ibdsimConfig["Specific_Density_Design"][i] == "True"):
#                self.ibdsimConfig["DiagnosticTables"] = ""
            
#        si l'option "1DProductwithoutM0" est true alors il faut que la dispersion soit géometrique ou Sichel, pour les auytres ca ne marche pas (pour l'instant)

#        Si il a des changements continus (ContinuousSizeChange), il ne peut pas y avoir de  de "zone" ni de "SpecificDensity Design"

        # si le Dispersal Distribution est dans la liste active le graph
        if self.ibdsimConfig["Dispersal_Distribution"][0] in ['Stepping Stone', 'Pareto', '0', '2', '3', '6', '7', '9', 'Geometric']:
            self.ui.viewGraphPushButton.setEnabled(True)

        # si barrier true alors vérifier:
        

        for i in range(len(self.ibdsimConfig["Barrier"])):
            if self.ibdsimConfig["Barrier"][i] == "True":
                if len(self.ibdsimConfig["X1_Barrier"][i]) == 0:
                    raise AttributeError("Please input a value for X1_Barrier", ["X1_Barrier"], i)
                if len(self.ibdsimConfig["X2_Barrier"][i]) == 0:
                    raise AttributeError("Please input a value for X2_Barrier", ["X2_Barrier"], i)
                if len(self.ibdsimConfig["Y1_Barrier"][i]) == 0:
                    raise AttributeError("Please input a value for Y1_Barrier", ["Y1_Barrier"], i)
                if len(self.ibdsimConfig["Y2_Barrier"][i]) == 0:
                    raise AttributeError("Please input a value for Y2_Barrier", ["Y2_Barrier"], i)
                if len(self.ibdsimConfig["Barrier_Crossing_Rate"][i]) == 0:
                    raise AttributeError("Please input a value for Barrier_Crossing_Rate", ["Barrier_Crossing_Rate"], i)
                
                if int(self.ibdsimConfig["X1_Barrier"][i]) > int(self.ibdsimConfig["Lattice_SizeX"][i]):
                    raise AttributeError("Found \"X1_Barrier\" > \"Lattice_SizeX\"",["X1_Barrier", "Lattice_SizeX"], i)
                if int(self.ibdsimConfig["X2_Barrier"][i]) > int(self.ibdsimConfig["Lattice_SizeX"][i]):
                    raise AttributeError("Found \"X2_Barrier\" > \"Lattice_SizeX\"",["X2_Barrier", "Lattice_SizeX"], i)
                if int(self.ibdsimConfig["Y1_Barrier"][i]) > int(self.ibdsimConfig["Lattice_SizeY"][i]):
                    raise AttributeError("Found \"Y1_Barrier\" > \"Lattice_SizeY\"",["Y1_Barrier", "Lattice_SizeY"], i)
                if int(self.ibdsimConfig["Y2_Barrier"][i]) > int(self.ibdsimConfig["Lattice_SizeY"][i]):
                    raise AttributeError("Found \"Y2_Barrier\" > \"Lattice_SizeX\"",["Y2_Barrier", "Lattice_SizeY"], i)
                
                # Les barrier doivent etre vertical ou horizontal soit X1 = X2 ou Y1 = Y2
                if int(self.ibdsimConfig["X1_Barrier"][i]) != int(self.ibdsimConfig["X2_Barrier"][i]) or int(self.ibdsimConfig["Y1_Barrier"][i]) != int(self.ibdsimConfig["Y2_Barrier"][i]):
                    raise AttributeError("The barrier is not horizontal or vertical line because X1 != X2 or Y1 != Y2, so IBDSim can not consider sush barriers",['X1_Barrier','X2_Barrier', 'Y1_Barrier', 'Y2_Barrier'],i)
                if int(self.ibdsimConfig["X1_Barrier"][i]) > int(self.ibdsimConfig["X2_Barrier"][i]) :
                    raise AttributeError("Exchanging 'X1_Barrier' and 'X2_Barrier' values, because X1 > X2",['X1_Barrier', 'X2_Barrier'],i)
                if int(self.ibdsimConfig["Y1_Barrier"][i]) > int(self.ibdsimConfig["Y2_Barrier"][i]) :
                    raise AttributeError("Exchanging 'Y1_Barrier' and 'Y2_Barrier' values, because Y1 > Y2",['Y1_Barrier', 'Y2_Barrier'],i)
                # test pour les barrier en point ou sur toute la longeur:
                
                
                
                if ( int(self.ibdsimConfig["X1_Barrier"][i]) == int(self.ibdsimConfig["X2_Barrier"][i]) ) and (float(self.ibdsimConfig["Barrier_Crossing_Rate"][i]) > 0.99999):
                    if (int(self.ibdsimConfig["Y1_Barrier"][i]) == int(self.ibdsimConfig["Y2_Barrier"][i]) ) and (float(self.ibdsimConfig["Barrier_Crossing_Rate"][i]) > 0.99999):
                        raise AttributeError("Barrier seetings will not have any effect because X1 = X2 and Y1 =Y2 or Barrier_Crossing_Rate > 0.99999",['X1_Barrier','X2_Barrier','Y1_Barrier','Y2_Barrier','Barrier_Crossing_Rate'],i)

                if (float(self.ibdsimConfig["Barrier_Crossing_Rate"][i] < 0.0000000001)):
                    if ( int(self.ibdsimConfig["X1_Barrier"][i]) == 1) and (int(self.ibdsimConfig["X2_Barrier"][i]) == int(self.ibdsimConfig["Lattice_SizeX"][i])):
                        if ( int(self.ibdsimConfig["Y1_Barrier"][i]) == 1) and (int(self.ibdsimConfig["Y2_Barrier"][i]) == int(self.ibdsimConfig["Lattice_SizeY"][i])):
                            raise AttributeError("Barrier_Grossing_Rate is < 0.0000000001 with a complete vertical or horizontal barrier ! \
                            with such configuration, either IBDSim will never stop a run because there is no possible MRCA if sample is taken across the barrier\
                            or the barrier will only cut the lattice into a reacheable part and a non-reacheable one\
                            We suspect it is not the desired configuration",['X1_Barrier','X2_Barrier','Y1_Barrier','Y2_Barrier','Barrier_Crossing_Rate'],i)


        self.ui.demographicTableWidget.setEnabled(False)
        self.ui.newDemoColPushButton.setEnabled (True)
        self.ui.saveDemoColPushButton.setEnabled (False)
        self.ui.modifyDemoColPushButton.setEnabled (True)
#        pour les changements linéaires et exponentiel, il doit y avoir une autre phase demo apres.
#        pas vrai pour logistique
#        t = len(self.ibdsimConfig["Continuous_Deme_Size_Variation"])-1
#        currentValue = str(eval("self.ui.Continuous_Deme_Size_VariationCombo%s.currentText()" % t))
#        if currentValue in ['Linear', 'Exponential']:
#            self.ui.newDemoColPushButton.emit(SIGNAL("clicked()"))
    if ACTIVE_GRAPH == True:
        def activeGraph(self):
            titleGraph = self.ibdsimConfig["Dispersal_Distribution"][0]
            xmax = self.ibdsimConfig["Lattice_SizeX"][0]
            M = self.ibdsimConfig["Total_Emigration_Rate"][0]
            n = None
            if titleGraph in ['Pareto', '0', '2', '3', '6', '7', '9']:
                n = self.ibdsimConfig["Pareto_Shape"][0]
            if titleGraph == "Geometric":
                n = self.ibdsimConfig["Geometric_Shape"][0]
            try:
                self.ui.graphicDispersal.deleteLater()
            except:
                pass
            self.ui.graphFrame.show()
            self.ui.graphicDispersal = Graph(latticeSizeX= xmax, totalEmigrationRate = M, Shape = n, nameGraph = titleGraph)
            self.ui.graphicDispersal.setObjectName("self.ui.graphicDispersal") 
            self.ui.graphGridLayout.addWidget(self.ui.graphicDispersal, 1, 1, 1, 1)
            self.ui.graphicDispersal.show()
            self.ui.hideGraphPushButton.setEnabled(True)
        def hideGraph(self):
            self.ui.graphicDispersal.deleteLater()
            self.ui.graphFrame.hide()
            self.ui.hideGraphPushButton.setEnabled(False)

    def changeTab(self):
        """ """
        currentTab = self.ui.tabWidget.currentWidget().objectName()
        if currentTab in ["welcome", "graphic_tab"]:
            self.ui.paramFrame.hide()
            self.setVisiblePreviewFile = False
        else:
            self.ui.paramFrame.show()
            self.setVisiblePreviewFile = True

    def hidePreviewFileConf(self):
        """ """
        if self.setVisiblePreviewFile == True:
            self.ui.paramFrame.hide()
            self.setVisiblePreviewFile = False
        else:
            self.ui.paramFrame.show()
            self.setVisiblePreviewFile = True

if ACTIVE_GRAPH == True: 
    class Graph(QwtPlot):
        def __init__(self, latticeSizeX= None, totalEmigrationRate = None, Shape = None, nameGraph = None):
            QwtPlot.__init__(self)
    
            xmax = int(latticeSizeX)
            M = float(totalEmigrationRate)
            
            Xlist = range(xmax+1)
            Xlist2 = Xlist[1:]
            Ylist = [0] 
            if nameGraph in ['Pareto', '0', '2', '3', '6', '7', '9']:
                n = float(Shape)
                for x in Xlist2:
                    y = M/(2* abs(x)**n)
                    Ylist.append(y)
    
            if nameGraph in ['Geometric']:
                n = float(Shape)
                for x in Xlist2:
                    y = (M/2)* n ** (abs(x)-1)
                    Ylist.append(y)
                    
            if nameGraph in ['Stepping Stone']:
                Ylist.append(M/2)
                for x in Xlist2:
                    Ylist.append(0)
            Xlist = Xlist[1:]
            Ylist = Ylist[1:]
            
            self.setBackgroundRole(QPalette.ColorRole("2"))
            self.setTitle("Dispersal Distribution")
            self.insertLegend(Qwt.QwtLegend(), Qwt.QwtPlot.RightLegend)
            self.setAxisTitle(Qwt.QwtPlot.xBottom, "Geographic Distance in Lattice Step")
            self.setAxisTitle(Qwt.QwtPlot.yLeft, "Probability")
            
            curve = QwtPlotCurve("Distribution %s" % nameGraph)
            curve.setData(Xlist, Ylist)
            curve.setSymbol(Qwt.QwtSymbol(Qwt.QwtSymbol.Star1,
                                          QBrush(),
                                          QPen(Qt.red, 0, Qt.DashLine),
                                          QSize(15,15)))
    
            curve.attach(self)
            self.replot()


def _test():
    """Fonction qui permet de lancer les docTests"""
    import doctest
    doctest.testmod()


def main():
    # Lancement des doctests
    _test()

    # sous linux, on appelle gconf pour voir les icones dans les menus et boutons
    if "linux" in sys.platform:
        cmd_args_list = ["gconftool-2", "--type", "boolean", "--set", "/desktop/gnome/interface/buttons_have_icons", "true"]
        Popen(cmd_args_list) 
        cmd_args_list = ["gconftool-2", "--type", "boolean", "--set", "/desktop/gnome/interface/menus_have_icons", "true"]
        Popen(cmd_args_list) 

    #print sys.argv+["-cmd", "-gnome-terminal"]
    nargv = sys.argv

    # instanciation des objets principaux
    app = QApplication(nargv)
    myapp = IBDSim(app)

    myapp.showMaximized()
    # les .app sous macos nécessitent cela pour que l'appli s'affiche en FG
    if "darwin" in sys.platform and ".app/" in sys.argv[0]:
        myapp.raise_()

    # lancement de la boucle Qt
    sys.exit(app.exec_())



if __name__ == '__main__':
    _test()
    main()
#    app = QApplication(sys.argv)
#    myapp = IBDSim(app)
#    myapp.showMaximized()
#   if "darwin" in sys.platform and ".app/" in sys.argv[0]:
#     	myapp.raise_()
#    sys.exit(app.exec_())

